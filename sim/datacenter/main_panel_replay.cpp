// -*- c-basic-offset: 4; indent-tabs-mode: nil -*-
// Standalone MoE dispatch replay over PanelTopology: reads per-layer flow
// files ("src dst bytes layer", 1-based ids), and per layer runs a dispatch
// phase followed by its combine (exact transpose), with a barrier between
// phases -- the GLASS MoEEP semantics, but over the panel fabrics and the
// frontend's routing policies + leasing OCS instead of a star topology.
// Transport is TcpSrc with CC disabled (cwnd = flowsize), mirroring the
// ASTRA frontend's -nocc mode. Deterministic; no randomness.
//
// Usage:
//   htsim_panel_replay -flowdir DIR -layers N -panel hybrid -planes 2
//     -linkGiBps 200 -latencyNs 1000 -policy directpref -ocs -reconfNs 10
//     [-q BYTES] [-maxwin BYTES] [-nocombine]
#include "config.h"
#include "eventlist.h"
#include "logfile.h"
#include "clock.h"
#include "tcp.h"
#include "mtcp.h"
#include "panel_topology.h"
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <set>
#include <utility>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

struct Flow { uint32_t s, d; uint64_t bytes; };

static EventList* g_ev = NULL;
static PanelTopology* g_top = NULL;
static Logfile* g_logfile = NULL;
static TcpRtxTimerScanner* g_rtx = NULL;

// config
static int g_policy = 2;              // 0 static, 1 adaptive, 2 directpref
static bool g_ocs = false;
static simtime_picosec g_reconf = 0;
static uint64_t g_maxwin = 2097152;
static double g_pref_factor = 1.10;

// leasing OCS ledgers (mirror of HTSimProtoTcp)
static vector<vector<simtime_picosec>> up_free, down_free;
static vector<vector<int>> up_peer, down_peer;
static uint64_t g_reconfigs = 0, g_reuses = 0;
static simtime_picosec g_wait_total = 0;

// phase state
static vector<vector<Flow>> g_layers;         // dispatch flows per layer
static bool g_combine = true;
static size_t g_phase = 0;                    // phase index: layer*2 (+1 combine)
static size_t g_outstanding = 0;
static int g_next_tag = 1;
static set<int> g_done_tags;
static vector<simtime_picosec> g_phase_end;
static uint64_t g_direct_bytes = 0, g_plane_bytes = 0;

static void launch_phase();
static void start_flow(uint32_t s, uint32_t d, uint64_t bytes, int forced_plane);

// ---- compiled-plan mode (SPECTRA et al.) --------------------------------
struct PlanCfg { std::vector<Flow> circ; };
static std::vector<std::vector<PlanCfg>> plan_cfgs;   // [plane][seq]
static std::vector<size_t> plan_cur;
static simtime_picosec plan_reconf = 0;
static bool plan_mode = false, plan_transposed = false;
static uint64_t plan_reconfigs = 0;

static bool plan_matching_changed(const PlanCfg& a, const PlanCfg& b) {
    std::set<std::pair<uint32_t,uint32_t>> sa, sb;
    for (size_t i = 0; i < a.circ.size(); i++) sa.insert({a.circ[i].s, a.circ[i].d});
    for (size_t i = 0; i < b.circ.size(); i++) sb.insert({b.circ[i].s, b.circ[i].d});
    return sa != sb;
}

class PlanDriver : public EventSource {
  public:
    PlanDriver(EventList& ev, int plane)
        : EventSource(ev, "plandrv"), _plane(plane) {}
    void schedule_at(simtime_picosec t) { eventlist().sourceIsPending(*this, t); }
    virtual void doNextEvent();
  private:
    int _plane;
};
static std::vector<PlanDriver*> plan_drivers;

static void plan_install(int plane) {
    size_t ci = plan_cur[plane];
    if (ci >= plan_cfgs[plane].size()) return;
    PlanCfg& c = plan_cfgs[plane][ci];
    double Bpns = 0; simtime_picosec drain = g_ev->now();
    for (size_t i = 0; i < c.circ.size(); i++) {
        uint32_t s = plan_transposed ? c.circ[i].d : c.circ[i].s;
        uint32_t d = plan_transposed ? c.circ[i].s : c.circ[i].d;
        start_flow(s, d, c.circ[i].bytes, plane);
    }
    // drain when the largest circuit finishes serializing on its uplink
    uint64_t mx = 0;
    for (size_t i = 0; i < c.circ.size(); i++) if (c.circ[i].bytes > mx) mx = c.circ[i].bytes;
    Bpns = 200.0 * 1073741824.0 / 1e9;   // plane rate, B/ns (GiB/s convention)
    drain = g_ev->now() + (simtime_picosec)((double)mx / Bpns * 1000.0)
          + timeFromNs(3.0 * 1500.0 / Bpns);
    plan_drivers[plane]->schedule_at(drain);
}

void PlanDriver::doNextEvent() {
    size_t ci = plan_cur[_plane];
    if (ci + 1 >= plan_cfgs[_plane].size()) return;
    bool ch = plan_matching_changed(plan_cfgs[_plane][ci], plan_cfgs[_plane][ci + 1]);
    plan_cur[_plane]++;
    if (ch) { plan_reconfigs++; }
    simtime_picosec when = eventlist().now() + (ch ? plan_reconf : 0);
    if (when <= eventlist().now()) { plan_install(_plane); }
    else {
        // dark period: install after T_r
        simtime_picosec t = when;
        // reuse the driver: schedule install via a zero-length config trick
        struct Inst : public EventSource {
            int pl; Inst(EventList& ev, int p) : EventSource(ev, "inst"), pl(p) {}
            virtual void doNextEvent() { plan_install(pl); }
        };
        Inst* iv = new Inst(eventlist(), _plane);
        eventlist().sourceIsPending(*iv, t);
    }
}

static void plan_launch_all() {
    for (size_t p = 0; p < plan_cfgs.size(); p++) {
        plan_cur[p] = 0;
        if (!plan_cfgs[p].empty()) plan_install((int)p);
    }
}
static void report_and_exit();

// --- compiled mode (QTP-style): per phase, offload the heaviest pairs to the
// planes until plane serialization balances DOR torus load, edge-color the
// offloaded pairs into matching epochs, and gate optical flows by epoch with
// drain-based advancement (matchings persist; T_r on epoch change).
static bool g_compiled = false;
static int g_forced_plane = -2;        // -2 = policy; -1 = direct; >=0 plane
struct EpochFlow { uint32_t s, d; uint64_t bytes; };
struct Flow_fwd;   // (Flow defined earlier)
static std::vector<std::vector<Flow>> g_plane_epochs[2];   // per plane chains
static size_t g_plane_cur[2] = {0, 0};
static size_t g_plane_outstanding[2] = {0, 0};
static std::map<int, int> g_epoch_tag_plane;
static int torus_hops_replay(uint32_t s, uint32_t d) {
    int ex = 8, hops = 0;
    int xs = s % ex, ys = s / ex, xd = d % ex, yd = d / ex;
    int dx = abs(xd - xs); dx = std::min(dx, ex - dx);
    int dy = abs(yd - ys); dy = std::min(dy, ex - dy);
    return dx + dy;
}
static void launch_epoch_on(int pl);

static bool g_flow_is_epoch = false;   // set while launching epoch flows

static void flow_done_cb(int /*src*/, int /*dst*/, int /*size*/, int tag) {
    if (g_done_tags.count(tag)) return;
    g_done_tags.insert(tag);
    if (g_compiled) {
        auto it = g_epoch_tag_plane.find(tag);
        if (it != g_epoch_tag_plane.end()) {
            int pl = it->second;
            g_epoch_tag_plane.erase(it);
            if (--g_plane_outstanding[pl] == 0) {
                g_plane_cur[pl]++; launch_epoch_on(pl);
            }
        }
    }
    if (--g_outstanding == 0) {
        if (plan_mode) return;      // plan main loop owns phase transitions
        g_phase_end.push_back(g_ev->now());
        g_phase++;
        launch_phase();
    }
}

static void start_flow(uint32_t s, uint32_t d, uint64_t bytes);
static void start_flow_forced(uint32_t s, uint32_t d, uint64_t bytes, int plane);

static void start_flow(uint32_t s, uint32_t d, uint64_t bytes) {
    start_flow(s, d, bytes, -2);
}
static void start_flow(uint32_t s, uint32_t d, uint64_t bytes, int forced_plane) {
    vector<PanelTopology::Candidate>* cands = g_top->get_candidates(s, d);
    simtime_picosec now = g_ev->now();
    if (forced_plane >= 0) {
        // plan mode: ride the installed circuit on this plane, no policy.
        for (size_t ci = 0; ci < cands->size(); ci++) {
            PanelTopology::Candidate& cd = (*cands)[ci];
            if (cd.is_plane && cd.plane == forced_plane) {
                g_plane_bytes += bytes;
                TcpSrc* src = new TcpSrc(NULL, NULL, *g_ev);
                MultipathTcpSrc* mtcp = new MultipathTcpSrc(UNCOUPLED, *g_ev, NULL);
                mtcp->setName("rpm_" + ntoa(g_next_tag));
                g_logfile->writeName(*mtcp);
                mtcp->addSubflow(src);
                src->_debug_srcid = (int)s; src->_debug_dstid = (int)d;
                src->set_flowsize(bytes);
                uint64_t win = bytes + 2 * Packet::data_packet_size();
                if (win > g_maxwin) win = g_maxwin;
                src->set_cwnd(win); src->set_ssthresh(win);
                TcpSink* snk = new TcpSink();
                snk->_debug_srcid = (int)s; snk->_debug_dstid = (int)d;
                snk->astrasim_flow_finish_recv_cb = &flow_done_cb;
                src->setName("rp_" + ntoa(s) + "_" + ntoa(d) + "_" + ntoa(g_next_tag));
                g_logfile->writeName(*src);
                snk->setName("rps_" + ntoa(s) + "_" + ntoa(d) + "_" + ntoa(g_next_tag));
                g_logfile->writeName(*snk);
                g_rtx->registerTcp(*src);
                Route* ro = new Route(*(cd.route)); ro->push_back(snk);
                Route* ri = new Route(); ri->push_back(src);
                int tag = g_next_tag++;
                // picosecond stagger: break exact event-time ties that can
                // wedge the event loop in same-timestamp storms
                src->connect(*ro, *ri, *snk, g_ev->now() + (simtime_picosec)(tag % 997));
                src->setFlowId(tag); snk->setFlowId(tag);
                break;
            }
        }
        for (size_t k = 0; k < cands->size(); k++) delete (*cands)[k].route;
        delete cands;
        return;
    }
    int best = -1; double best_cost = 0; bool best_reuse = false;
    simtime_picosec best_start = now;
    int direct_idx = -1; double direct_cost = 0;
    if (g_forced_plane != -2) {
        for (size_t ci = 0; ci < cands->size(); ci++) {
            PanelTopology::Candidate& cd = (*cands)[ci];
            if (g_forced_plane < 0 && !cd.is_plane) { best = (int)ci; break; }
            if (g_forced_plane >= 0 && cd.is_plane && cd.plane == g_forced_plane) {
                best = (int)ci; break;
            }
        }
        assert(best >= 0);
    } else {
    for (size_t ci = 0; ci < cands->size(); ci++) {
        PanelTopology::Candidate& cd = (*cands)[ci];
        double cost; bool reuse = false; simtime_picosec tstart = now;
        (void)0;
        if (!cd.is_plane) {
            double ser = 0;
            for (size_t k = 0; k < cd.hop_queues.size(); k++) {
                LedgerQueue* q = cd.hop_queues[k];
                double Bpns = (double)q->link_bitrate() / 8.0 / 1e9;
                double queued = (g_policy != 0) ? (double)q->reserved_bytes() : 0.0;
                double t = (queued + (double)bytes) / Bpns;
                if (t > ser) ser = t;
            }
            cost = (double)cd.latency_sum / 1000.0 + ser;
            direct_cost = cost; direct_idx = (int)ci;
        } else if (g_ocs) {
            int pl = cd.plane;
            simtime_picosec t0 = now;
            if (up_free[pl][s] > t0) t0 = up_free[pl][s];
            if (down_free[pl][d] > t0) t0 = down_free[pl][d];
            reuse = (up_peer[pl][s] == (int)d && down_peer[pl][d] == (int)s);
            tstart = t0 + (reuse ? 0 : g_reconf);
            double Bpns = (double)cd.hop_queues[0]->link_bitrate() / 8.0 / 1e9;
            if (g_policy == 0) {
                cost = (double)g_reconf / 1000.0
                     + (double)cd.latency_sum / 1000.0 + (double)bytes / Bpns;
            } else {
                cost = (double)(tstart - now) / 1000.0
                     + (double)cd.latency_sum / 1000.0 + (double)bytes / Bpns;
            }
        } else {
            double ser = 0;
            for (size_t k = 0; k < cd.hop_queues.size(); k++) {
                LedgerQueue* q = cd.hop_queues[k];
                double Bpns = (double)q->link_bitrate() / 8.0 / 1e9;
                double queued = (g_policy != 0) ? (double)q->reserved_bytes() : 0.0;
                double t = (queued + (double)bytes) / Bpns;
                if (t > ser) ser = t;
            }
            cost = (double)cd.latency_sum / 1000.0 + ser;
        }
        if (best < 0 || cost < best_cost) {
            best = (int)ci; best_cost = cost; best_reuse = reuse; best_start = tstart;
        }
    }
    if (g_policy == 2 && direct_idx >= 0 && best != direct_idx &&
        direct_cost <= g_pref_factor * best_cost) {
        best = direct_idx;
    }
    }
    PanelTopology::Candidate* choice = &(*cands)[best];
    (void)choice;
    simtime_picosec flow_delay = 0;
    if (choice->is_plane && g_ocs) {
        int pl = choice->plane;
        double Bpns = (double)choice->hop_queues[0]->link_bitrate() / 8.0 / 1e9;
        simtime_picosec occ = (simtime_picosec)(((double)bytes / Bpns) * 1000.0)
                              + timeFromNs(3.0 * 1500.0 / Bpns);
        simtime_picosec rel = best_start + occ;
        up_free[pl][s] = rel; down_free[pl][d] = rel;
        up_peer[pl][s] = (int)d; down_peer[pl][d] = (int)s;
        if (best_reuse) g_reuses++; else g_reconfigs++;
        g_wait_total += (best_start - now);
        flow_delay = best_start - now;
        g_plane_bytes += bytes;
    } else {
        for (size_t k = 0; k < choice->hop_queues.size(); k++)
            choice->hop_queues[k]->reserve_bytes(bytes);
        if (choice->is_plane) g_plane_bytes += bytes; else g_direct_bytes += bytes;
    }

    TcpSrc* src = new TcpSrc(NULL, NULL, *g_ev);
    MultipathTcpSrc* mtcp = new MultipathTcpSrc(UNCOUPLED, *g_ev, NULL);
    mtcp->setName("rpm_" + ntoa(g_next_tag));
    g_logfile->writeName(*mtcp);
    mtcp->addSubflow(src);
    src->_debug_srcid = (int)s; src->_debug_dstid = (int)d;
    src->set_flowsize(bytes);
    uint64_t win = bytes + 2 * Packet::data_packet_size();
    if (win > g_maxwin) win = g_maxwin;
    src->set_cwnd(win); src->set_ssthresh(win);
    TcpSink* snk = new TcpSink();
    snk->_debug_srcid = (int)s; snk->_debug_dstid = (int)d;
    snk->astrasim_flow_finish_recv_cb = &flow_done_cb;
    src->setName("rp_" + ntoa(s) + "_" + ntoa(d) + "_" + ntoa(g_next_tag));
    g_logfile->writeName(*src);
    snk->setName("rps_" + ntoa(s) + "_" + ntoa(d) + "_" + ntoa(g_next_tag));
    g_logfile->writeName(*snk);
    g_rtx->registerTcp(*src);
    Route* routeout = new Route(*(choice->route));
    routeout->push_back(snk);
    Route* routein = new Route();
    routein->push_back(src);
    src->connect(*routeout, *routein, *snk, g_ev->now() + flow_delay);
    src->setFlowId(g_next_tag); snk->setFlowId(g_next_tag);
    if (g_flow_is_epoch) g_epoch_tag_plane[g_next_tag] = g_forced_plane;
    g_next_tag++;
    for (size_t k = 0; k < cands->size(); k++) delete (*cands)[k].route;
    delete cands;
}

static void launch_phase() {
    size_t nphases = g_layers.size() * (g_combine ? 2 : 1);
    if (g_phase >= nphases) { report_and_exit(); }
    size_t layer = g_combine ? g_phase / 2 : g_phase;
    bool combine = g_combine && (g_phase % 2 == 1);
    vector<Flow>& fl = g_layers[layer];
    g_outstanding = fl.size();
    if (g_outstanding == 0) { g_phase_end.push_back(g_ev->now()); g_phase++; launch_phase(); return; }
    if (!g_compiled) {
        for (size_t i = 0; i < fl.size(); i++) {
            if (combine) start_flow(fl[i].d, fl[i].s, fl[i].bytes);
            else         start_flow(fl[i].s, fl[i].d, fl[i].bytes);
        }
        return;
    }
    // compiled: offset-grouped matchings. Each destination offset o is a
    // perfect matching i -> (i+o) mod N; offload the heaviest hop-weighted
    // offsets to the planes until plane serialization balances DOR torus
    // load; each plane runs its own epoch chain (one offset per epoch),
    // matchings persist for the whole epoch.
    const uint32_t NN = 64;
    std::vector<long double> offW(NN, 0.0L);
    std::vector<std::vector<Flow>> offFlows(NN);
    long double torus_hopbytes = 0;
    for (size_t i = 0; i < fl.size(); i++) {
        uint32_t s = combine ? fl[i].d : fl[i].s;
        uint32_t d = combine ? fl[i].s : fl[i].d;
        uint32_t o = (d + NN - s) % NN;
        int h = torus_hops_replay(s, d);
        offW[o] += (long double)h * fl[i].bytes;
        Flow f; f.s = s; f.d = d; f.bytes = fl[i].bytes;
        offFlows[o].push_back(f);
        torus_hopbytes += (long double)h * fl[i].bytes;
    }
    std::vector<uint32_t> order;
    for (uint32_t o = 1; o < NN; o++) if (!offFlows[o].empty()) order.push_back(o);
    std::sort(order.begin(), order.end(),
              [&](uint32_t a, uint32_t b) { return offW[a] > offW[b]; });
    g_plane_epochs[0].clear(); g_plane_epochs[1].clear();
    long double plane_sel = 0;
    size_t k = 0;
    for (; k < order.size(); k++) {
        uint64_t ob = 0;
        for (size_t j = 0; j < offFlows[order[k]].size(); j++)
            ob += offFlows[order[k]][j].bytes;
        long double torus_load = torus_hopbytes / (4.0L * NN);
        long double plane_load = (plane_sel + ob) / 2.0L / NN;
        if (plane_load >= torus_load) break;
        plane_sel += ob;
        torus_hopbytes -= offW[order[k]];
        g_plane_epochs[k % 2].push_back(offFlows[order[k]]);
    }
    for (size_t i = k; i < order.size(); i++)
        for (size_t j = 0; j < offFlows[order[i]].size(); j++) {
            Flow& f = offFlows[order[i]][j];
            start_flow_forced(f.s, f.d, f.bytes, -1);
        }
    for (int pl = 0; pl < 2; pl++) { g_plane_cur[pl] = 0; launch_epoch_on(pl); }
}

static void launch_epoch_on(int pl) {
    if (g_plane_cur[pl] >= g_plane_epochs[pl].size()) return;
    std::vector<Flow>& ef = g_plane_epochs[pl][g_plane_cur[pl]];
    g_plane_outstanding[pl] = ef.size();
    for (size_t i = 0; i < ef.size(); i++)
        start_flow_forced(ef[i].s, ef[i].d, ef[i].bytes, pl);
}


static void start_flow_forced(uint32_t s, uint32_t d, uint64_t bytes, int plane) {
    g_forced_plane = plane;
    g_flow_is_epoch = (plane >= 0);
    start_flow(s, d, bytes);
    g_flow_is_epoch = false;
    g_forced_plane = -2;
}

int main(int argc, char** argv) {
    string flowdir, flowlist, plan_file; int nlayers = 1, nodes = 64, planes = 2;
    double link_gibps = 200, plane_gibps = -1;
    simtime_picosec lat = timeFromNs(1000);
    mem_b qsize = 90000 * 1500;
    string panel = "hybrid", policy = "directpref";
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-flowdir")) flowdir = argv[++i];
        else if (!strcmp(argv[i], "-flowlist")) flowlist = argv[++i];
        else if (!strcmp(argv[i], "-layers")) nlayers = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-nodes")) nodes = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-panel")) panel = argv[++i];
        else if (!strcmp(argv[i], "-planes")) planes = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-linkGiBps")) link_gibps = atof(argv[++i]);
        else if (!strcmp(argv[i], "-planeGiBps")) plane_gibps = atof(argv[++i]);
        else if (!strcmp(argv[i], "-latencyNs")) lat = timeFromNs(atof(argv[++i]));
        else if (!strcmp(argv[i], "-policy")) policy = argv[++i];
        else if (!strcmp(argv[i], "-ocs")) g_ocs = true;
        else if (!strcmp(argv[i], "-reconfNs")) g_reconf = timeFromNs(atof(argv[++i]));
        else if (!strcmp(argv[i], "-q")) qsize = (mem_b)atol(argv[++i]) * 1500;
        else if (!strcmp(argv[i], "-maxwin")) g_maxwin = atol(argv[++i]);
        else if (!strcmp(argv[i], "-nocombine")) g_combine = false;
        else if (!strcmp(argv[i], "-plan")) { plan_mode = true; plan_file = argv[++i]; }
        else if (!strcmp(argv[i], "-planreconfNs")) plan_reconf = timeFromNs(atof(argv[++i]));
        else if (!strcmp(argv[i], "-compiled")) g_compiled = true;
        else { fprintf(stderr, "unknown arg %s\n", argv[i]); return 1; }
    }
    if (plane_gibps < 0) plane_gibps = link_gibps;
    g_policy = (policy == "static") ? 0 : (policy == "adaptive") ? 1 : 2;

    EventList eventlist; g_ev = &eventlist;
    eventlist.setEndtime(timeFromSec(60));
    Clock c(timeFromSec(5 / 100.), eventlist);
    Logfile lf("replay_logout.dat", eventlist); g_logfile = &lf;
    TcpRtxTimerScanner rtx(timeFromMs(250), eventlist); g_rtx = &rtx;

    PanelTopology::Base base =
        (panel == "torus3d") ? PanelTopology::Base::Torus3D :
        (panel == "mesh3d")  ? PanelTopology::Base::Mesh3D :
        (panel == "fullswitch") ? PanelTopology::Base::None :
        PanelTopology::Base::Torus2D;
    int p = (panel == "fullswitch") ? 6 :
            (panel == "torus3d" || panel == "mesh3d") ? 0 : planes;
    g_top = new PanelTopology(nodes, base, p, link_gibps, lat,
                              plane_gibps, lat, qsize, &lf, &eventlist);
    up_free.assign(p, vector<simtime_picosec>(nodes, 0));
    down_free.assign(p, vector<simtime_picosec>(nodes, 0));
    up_peer.assign(p, vector<int>(nodes, -1));
    down_peer.assign(p, vector<int>(nodes, -1));

    vector<string> phase_files;
    if (!flowlist.empty()) {
        // manifest mode: each line is a flow-file path = one dispatch phase
        // (plus its combine), in execution order -- e.g. iteration x layer
        // chains from gen_iter_replay.py. OCS ledgers persist across phases.
        ifstream mf(flowlist.c_str());
        string line;
        while (getline(mf, line)) if (!line.empty()) phase_files.push_back(line);
    } else {
        for (int L = 0; L < nlayers; L++) {
            char path[512];
            snprintf(path, sizeof(path), "%s/dispatch_layer%d.htsim", flowdir.c_str(), L);
            phase_files.push_back(path);
        }
    }
    for (size_t pi = 0; pi < phase_files.size(); pi++) {
        ifstream f(phase_files[pi].c_str());
        vector<Flow> fl;
        if (f.good()) {
            string line;
            while (getline(f, line)) {
                if (line.empty()) continue;
                istringstream is(line);
                long s, d; unsigned long long b; int lid;
                if (!(is >> s >> d >> b >> lid)) continue;
                Flow fw; fw.s = (uint32_t)(s - 1); fw.d = (uint32_t)(d - 1); fw.bytes = b;
                fl.push_back(fw);
            }
        }
        g_layers.push_back(fl);
    }
    size_t total_flows = 0;
    for (size_t i = 0; i < g_layers.size(); i++) total_flows += g_layers[i].size();
    printf("REPLAY start: layers=%zu flows/layer0=%zu total=%zu panel=%s policy=%s ocs=%d\n",
           g_layers.size(), g_layers.empty() ? 0 : g_layers[0].size(),
           total_flows, panel.c_str(), policy.c_str(), (int)g_ocs);

    if (plan_mode) {
        // parse flat plan program
        ifstream pf(plan_file.c_str());
        string tok; int nplanes_plan = 0;
        pf >> tok >> nplanes_plan;
        plan_cfgs.assign(nplanes_plan, std::vector<PlanCfg>());
        plan_cur.assign(nplanes_plan, 0);
        int curp = -1;
        string line;
        while (pf >> tok) {
            if (tok == "P") { pf >> curp; plan_cfgs[curp].push_back(PlanCfg()); }
            else if (tok == "C") {
                long a, b2; unsigned long long by; pf >> a >> b2 >> by;
                Flow fw; fw.s = (uint32_t)a; fw.d = (uint32_t)b2; fw.bytes = by;
                plan_cfgs[curp].back().circ.push_back(fw);
            }
        }
        size_t tot = 0;
        for (size_t p = 0; p < plan_cfgs.size(); p++) {
            plan_drivers.push_back(new PlanDriver(eventlist, (int)p));
            for (size_t c2 = 0; c2 < plan_cfgs[p].size(); c2++) tot += plan_cfgs[p][c2].circ.size();
        }
        printf("PLAN mode: planes=%zu total_flows=%zu (x2 with combine)\n",
               plan_cfgs.size(), tot);
        g_outstanding = tot;
        g_phase = 0;
        plan_transposed = false;
        plan_launch_all();
        uint64_t evn = 0;
        while (eventlist.doNextEvent()) {
            if (eventlist.now() > timeFromSec(2)) {
                fprintf(stderr, "SAFETY: sim exceeded 2s, terminating (out=%zu)\n", g_outstanding);
                break;
            }
            if (++evn % 2000000 == 0)
                fprintf(stderr, "HB ev=%lluM now_ns=%.0f out=%zu ph=%zu\n",
                        (unsigned long long)(evn / 1000000),
                        timeAsNs(eventlist.now()), g_outstanding, g_phase_end.size());
            if (g_outstanding == 0 && g_combine && !plan_transposed) {
                g_phase_end.push_back(eventlist.now());
                plan_transposed = true;
                g_outstanding = tot;
                plan_launch_all();
            } else if (g_outstanding == 0 && (plan_transposed || !g_combine)
                       && g_phase_end.size() < (g_combine ? 2u : 1u)) {
                g_phase_end.push_back(eventlist.now());
            }
        }
        printf("PLAN_RESULT makespan_ns=%.0f reconfigs=%llu plane_bytes=%llu rtx=%llu\n",
               g_phase_end.empty() ? -1.0 : timeAsNs(g_phase_end.back()),
               (unsigned long long)plan_reconfigs,
               (unsigned long long)g_plane_bytes,
               (unsigned long long)TcpSrc::_global_rtx_count);
        return 0;
    }
    launch_phase();
    while (eventlist.doNextEvent()) {}
    report_and_exit();
    return 0;
}

static void report_and_exit() {
    simtime_picosec prev = 0;
    for (size_t i = 0; i < g_phase_end.size(); i++) {
        printf("PHASE %zu (%s L%zu) end_ns=%.0f dur_ns=%.0f\n", i,
               (g_combine && i % 2) ? "combine" : "dispatch",
               g_combine ? i / 2 : i,
               timeAsNs(g_phase_end[i]), timeAsNs(g_phase_end[i] - prev));
        prev = g_phase_end[i];
    }
    printf("REPLAY_RESULT makespan_ns=%.0f phases=%zu direct_bytes=%llu plane_bytes=%llu "
           "reconfigs=%llu reuses=%llu wait_ns=%.0f rtx=%llu\n",
           g_phase_end.empty() ? -1.0 : timeAsNs(g_phase_end.back()),
           g_phase_end.size(),
           (unsigned long long)g_direct_bytes, (unsigned long long)g_plane_bytes,
           (unsigned long long)g_reconfigs, (unsigned long long)g_reuses,
           timeAsNs(g_wait_total), (unsigned long long)TcpSrc::_global_rtx_count);
    fflush(stdout);
    exit(0);
}
