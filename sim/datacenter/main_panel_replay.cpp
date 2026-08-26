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
#include "panel_topology.h"
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <set>
#include <vector>
#include <string>

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
static void report_and_exit();

static void flow_done_cb(int /*src*/, int /*dst*/, int /*size*/, int tag) {
    if (g_done_tags.count(tag)) return;
    g_done_tags.insert(tag);
    if (--g_outstanding == 0) {
        g_phase_end.push_back(g_ev->now());
        g_phase++;
        launch_phase();
    }
}

static void start_flow(uint32_t s, uint32_t d, uint64_t bytes) {
    vector<PanelTopology::Candidate>* cands = g_top->get_candidates(s, d);
    simtime_picosec now = g_ev->now();
    int best = -1; double best_cost = 0; bool best_reuse = false;
    simtime_picosec best_start = now;
    int direct_idx = -1; double direct_cost = 0;
    for (size_t ci = 0; ci < cands->size(); ci++) {
        PanelTopology::Candidate& cd = (*cands)[ci];
        double cost; bool reuse = false; simtime_picosec tstart = now;
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
    PanelTopology::Candidate* choice = &(*cands)[best];
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
    for (size_t i = 0; i < fl.size(); i++) {
        if (combine) start_flow(fl[i].d, fl[i].s, fl[i].bytes);
        else         start_flow(fl[i].s, fl[i].d, fl[i].bytes);
    }
}

int main(int argc, char** argv) {
    string flowdir; int nlayers = 1, nodes = 64, planes = 2;
    double link_gibps = 200, plane_gibps = -1;
    simtime_picosec lat = timeFromNs(1000);
    mem_b qsize = 90000 * 1500;
    string panel = "hybrid", policy = "directpref";
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-flowdir")) flowdir = argv[++i];
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

    for (int L = 0; L < nlayers; L++) {
        char path[512];
        snprintf(path, sizeof(path), "%s/dispatch_layer%d.htsim", flowdir.c_str(), L);
        ifstream f(path);
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
