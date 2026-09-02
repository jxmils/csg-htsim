// -*- c-basic-offset: 4; indent-tabs-mode: nil -*-
#include "panel_topology.h"
#include <fstream>
#include <iostream>
#include "main.h"
#include <sstream>
#include <cassert>
#include <cmath>

using namespace std;

static string ntos(int v) { stringstream s; s << v; return s.str(); }

PanelTopology::PanelTopology(uint32_t npus, Base base, int planes,
                             double base_gibps, simtime_picosec base_latency,
                             double plane_gibps, simtime_picosec plane_latency,
                             mem_b queuesize, Logfile* logfile, EventList* ev,
                             const std::vector<int>& extents, bool quiet,
                             const std::string& graphfile)
    : _n(npus), _base(base), _planes(planes), _queuesize(queuesize),
      _logfile(logfile), _ev(ev), _quiet(quiet) {
    switch (base) {
        case Base::None:    _dims = 0; _wrap = false; break;
        case Base::Custom:  _dims = 0; _wrap = false; _graphfile = graphfile; break;
        case Base::Ring1D:  _dims = 1; _wrap = true;  break;
        case Base::Mesh2D:  _dims = 2; _wrap = false; break;
        case Base::Torus2D: _dims = 2; _wrap = true;  break;
        case Base::Mesh3D:  _dims = 3; _wrap = false; break;
        case Base::Torus3D: _dims = 3; _wrap = true;  break;
        case Base::RingRows: _dims = 2; _wrap = true; break;
    }
    if (_dims > 0) {
        if (!extents.empty() && (int)extents.size() != _dims) {
            cerr << "PanelTopology: extents size mismatch" << endl; exit(1);
        }
        if (!extents.empty()) _extents = extents;
        else _extents.assign(_dims,
                             (int)lround(pow((double)npus, 1.0 / _dims)));
        int check = 1;
        for (int d = 0; d < _dims; d++) check *= _extents[d];
        if ((uint32_t)check != npus) {
            cerr << "PanelTopology: npus " << npus << " is not a perfect "
                 << _dims << "-dim grid" << endl;
            exit(1);
        }
        build_base(base_gibps, base_latency);
    }
    if (_base == Base::Custom) build_custom(base_gibps, base_latency);
    if (_planes > 0) build_planes(plane_gibps, plane_latency);
    assert(_dims > 0 || _planes > 0 || _base == Base::Custom);
}

LedgerQueue* PanelTopology::make_queue(double gibps, const string& name) {
    QueueLoggerSampling* ql = NULL;
    if (!_quiet) {
        ql = new QueueLoggerSampling(timeFromUs((uint32_t)1000), *_ev);
        _logfile->addLogger(*ql);
    }
    LedgerQueue* q = new LedgerQueue(speedFromGiBps(gibps), _queuesize, *_ev, ql);
    q->setName(name);
    _logfile->writeName(*q);
    return q;
}

Pipe* PanelTopology::make_pipe(simtime_picosec lat, const string& name) {
    Pipe* p = new Pipe(lat, *_ev);
    p->setName(name);
    _logfile->writeName(*p);
    return p;
}

int PanelTopology::coord(uint32_t id, int dim) const {
    uint32_t v = id;
    for (int d = 0; d < dim; d++) v /= _extents[d];
    return (int)(v % _extents[dim]);
}

uint32_t PanelTopology::id_of(const vector<int>& c) const {
    uint32_t id = 0;
    for (int d = _dims - 1; d >= 0; d--) id = id * _extents[d] + c[d];
    return id;
}

void PanelTopology::build_base(double gibps, simtime_picosec lat) {
    // port index: 2*dim + 0 for +direction, 2*dim + 1 for -direction
    _dir_q.assign(_n, vector<LedgerQueue*>(2 * _dims, (LedgerQueue*)NULL));
    _dir_p.assign(_n, vector<Pipe*>(2 * _dims, (Pipe*)NULL));
    for (uint32_t node = 0; node < _n; node++) {
        for (int d = 0; d < _dims; d++) {
            if (_base == Base::RingRows && d > 0) continue;  // no column links
            int c = coord(node, d);
            bool plus_exists = (c + 1 < _extents[d]) || _wrap;
            bool minus_exists = (c > 0) || _wrap;
            if (_extents[d] <= 2) {
                // extent 2: +1 and -1 reach the same neighbour; wrap adds no
                // second edge. Only the + port exists (matches Mesh2D's
                // "wraparound && width > 2" guard in the analytical study).
                minus_exists = false;
                plus_exists = (c + 1 < _extents[d]) || _wrap;
            }
            if (plus_exists) {
                _dir_q[node][2 * d] = make_queue(gibps,
                    "PQ" + ntos(node) + "_d" + ntos(d) + "p");
                _dir_p[node][2 * d] = make_pipe(lat,
                    "PP" + ntos(node) + "_d" + ntos(d) + "p");
            }
            if (minus_exists) {
                _dir_q[node][2 * d + 1] = make_queue(gibps,
                    "PQ" + ntos(node) + "_d" + ntos(d) + "m");
                _dir_p[node][2 * d + 1] = make_pipe(lat,
                    "PP" + ntos(node) + "_d" + ntos(d) + "m");
            }
        }
    }
}

void PanelTopology::build_planes(double gibps, simtime_picosec lat) {
    _up_q.assign(_planes, vector<LedgerQueue*>(_n, (LedgerQueue*)NULL));
    _up_p.assign(_planes, vector<Pipe*>(_n, (Pipe*)NULL));
    _down_q.assign(_planes, vector<LedgerQueue*>(_n, (LedgerQueue*)NULL));
    _down_p.assign(_planes, vector<Pipe*>(_n, (Pipe*)NULL));
    for (int pl = 0; pl < _planes; pl++) {
        for (uint32_t node = 0; node < _n; node++) {
            _up_q[pl][node] = make_queue(gibps, "UQ" + ntos(pl) + "_" + ntos(node));
            _up_p[pl][node] = make_pipe(lat, "UP" + ntos(pl) + "_" + ntos(node));
            _down_q[pl][node] = make_queue(gibps, "DQ" + ntos(pl) + "_" + ntos(node));
            _down_p[pl][node] = make_pipe(lat, "DP" + ntos(pl) + "_" + ntos(node));
        }
    }
}

int PanelTopology::step_towards(int& cur, int target, bool tie_backward,
                                int extent) const {
    int port;
    if (!_wrap) {
        if (target > cur) { cur++; port = 0; }
        else { cur--; port = 1; }
    } else {
        int fwd = ((target - cur) % extent + extent) % extent;
        int bwd = extent - fwd;
        bool go_fwd;
        if (fwd < bwd) go_fwd = true;
        else if (bwd < fwd) go_fwd = false;
        else go_fwd = !tie_backward;   // antipodal: split by source parity
        if (go_fwd) { cur = (cur + 1) % extent; port = 0; }
        else { cur = (cur - 1 + extent) % extent; port = 1; }
    }
    return port;
}

void PanelTopology::build_custom(double gibps, simtime_picosec lat) {
    std::ifstream f(_graphfile.c_str());
    if (!f.good()) { std::cerr << "custom graph file open failed: " << _graphfile << std::endl; abort(); }
    // Two passes: the graph may carry non-endpoint devices (HammingMesh row and
    // column switches) whose ids run past the endpoint count, so size to the
    // highest id present rather than to _n.
    std::vector<std::pair<long,long>> edges;
    // Optional third column: per-link rate in GiB/s ("E a b 4800"). Links
    // without it take -linkGiBps, so files without the column are unchanged.
    std::vector<double> ebw;
    long maxid = (long)_n - 1;
    std::string tok;
    while (f >> tok) {
        if (tok == "E") {
            long a, b; f >> a >> b;
            std::string rest; std::getline(f, rest);
            double w = 0.0;
            { std::istringstream rs(rest); rs >> w; if (rs.fail()) w = 0.0; }
            edges.push_back(std::make_pair(a, b));
            ebw.push_back(w);
            if (a > maxid) maxid = a;
            if (b > maxid) maxid = b;
        } else { std::string rest; std::getline(f, rest); }
    }
    _ndev = (uint32_t)(maxid + 1);
    _adj.assign(_ndev, std::vector<uint32_t>());
    std::vector<std::vector<double>> adjbw(_ndev);
    size_t n_explicit = 0;
    for (size_t i = 0; i < edges.size(); i++) {
        _adj[edges[i].first].push_back((uint32_t)edges[i].second);
        adjbw[edges[i].first].push_back(ebw[i]);
        if (ebw[i] > 0.0) n_explicit++;
    }
    std::cerr << "custom graph: " << edges.size() << " directed links, "
              << n_explicit << " with explicit per-link GiB/s (rest at -linkGiBps "
              << gibps << ")" << std::endl;
    _dir_q.assign(_ndev, std::vector<LedgerQueue*>());
    _dir_p.assign(_ndev, std::vector<Pipe*>());
    for (uint32_t u = 0; u < _ndev; u++) {
        for (size_t p = 0; p < _adj[u].size(); p++) {
            char nm[64];
            snprintf(nm, sizeof(nm), "cq_%u_%zu", u, p);
            double g = adjbw[u][p] > 0.0 ? adjbw[u][p] : gibps;
            _dir_q[u].push_back(make_queue(g, nm));
            snprintf(nm, sizeof(nm), "cp_%u_%zu", u, p);
            _dir_p[u].push_back(make_pipe(lat, nm));
        }
    }
    // per-source BFS next-hop port table
    _nh.assign(_ndev, std::vector<int>(_ndev, -1));
    for (uint32_t s = 0; s < _ndev; s++) {
        std::vector<int> par(_ndev, -1), parport(_ndev, -1);
        std::vector<uint32_t> q; q.push_back(s); par[s] = (int)s;
        for (size_t qi = 0; qi < q.size(); qi++) {
            uint32_t u = q[qi];
            for (size_t p = 0; p < _adj[u].size(); p++) {
                uint32_t v = _adj[u][p];
                if (par[v] < 0) { par[v] = (int)u; parport[v] = (int)p; q.push_back(v); }
            }
        }
        for (uint32_t d = 0; d < _ndev; d++) {
            if (d == s || par[d] < 0) continue;
            uint32_t cur = d;
            while ((uint32_t)par[cur] != s) cur = (uint32_t)par[cur];
            _nh[s][d] = parport[cur];
        }
    }
}

PanelTopology::Candidate PanelTopology::direct_candidate(uint32_t src, uint32_t dest) {
    Candidate cand;
    cand.route = new Route();
    cand.is_plane = false; cand.plane = -1;
    cand.hops = 0; cand.latency_sum = 0;

    if (_base == Base::Custom) {
        uint32_t cur = src;
        while (cur != dest) {
            int port = _nh[cur][dest];
            assert(port >= 0);
            LedgerQueue* q = _dir_q[cur][port];
            Pipe* p = _dir_p[cur][port];
            cand.route->push_back(q); cand.route->push_back(p);
            cand.hop_queues.push_back(q);
            cand.latency_sum += p->delay();
            cand.hops++;
            cur = _adj[cur][port];
        }
        check_non_null(cand.route);
        return cand;
    }
    vector<int> c(_dims), t(_dims);
    for (int d = 0; d < _dims; d++) { c[d] = coord(src, d); t[d] = coord(dest, d); }

    // dimension-order X -> Y (-> Z); parity tie-break keyed on the SOURCE
    // coordinate in that dimension (mirrors Hybrid2D/Mesh3D).
    for (int d = 0; d < _dims; d++) {
        bool tie_backward = (coord(src, d) & 1) != 0;
        while (c[d] != t[d]) {
            uint32_t at = id_of(c);
            int port = step_towards(c[d], t[d], tie_backward, _extents[d]);
            int pidx = 2 * d + port;
            if (_extents[d] == 2 && _dir_q[at][pidx] == NULL) pidx = 2 * d;  // folded pair
            LedgerQueue* q = _dir_q[at][pidx];
            Pipe* p = _dir_p[at][pidx];
            assert(q && p);
            cand.route->push_back(q);
            cand.route->push_back(p);
            cand.hop_queues.push_back(q);
            cand.latency_sum += p->delay();
            cand.hops++;
        }
    }
    check_non_null(cand.route);
    return cand;
}

PanelTopology::Candidate PanelTopology::plane_candidate(uint32_t src, uint32_t dest, int pl) {
    Candidate cand;
    cand.route = new Route();
    cand.is_plane = true; cand.plane = pl;
    cand.route->push_back(_up_q[pl][src]);
    cand.route->push_back(_up_p[pl][src]);
    cand.route->push_back(_down_q[pl][dest]);
    cand.route->push_back(_down_p[pl][dest]);
    cand.hop_queues.push_back(_up_q[pl][src]);
    cand.hop_queues.push_back(_down_q[pl][dest]);
    cand.hops = 2;
    cand.latency_sum = _up_p[pl][src]->delay() + _down_p[pl][dest]->delay();
    check_non_null(cand.route);
    return cand;
}

vector<PanelTopology::Candidate>* PanelTopology::get_candidates(uint32_t src, uint32_t dest) {
    assert(src < _n && dest < _n && src != dest);
    vector<Candidate>* out = new vector<Candidate>();
    bool direct_ok = (_dims > 0) || _base == Base::Custom;
    if (_base == Base::RingRows && coord(src, 1) != coord(dest, 1))
        direct_ok = false;   // rows are disjoint rings; cross-row is optical-only
    if (direct_ok) out->push_back(direct_candidate(src, dest));
    for (int pl = 0; pl < _planes; pl++) out->push_back(plane_candidate(src, dest, pl));
    return out;
}

vector<const Route*>* PanelTopology::get_bidir_paths(uint32_t src, uint32_t dest, bool reverse) {
    vector<Candidate>* cands = get_candidates(src, dest);
    vector<const Route*>* paths = new vector<const Route*>();
    for (size_t i = 0; i < cands->size(); i++) paths->push_back((*cands)[i].route);
    delete cands;
    return paths;
}
