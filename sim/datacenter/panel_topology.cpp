// -*- c-basic-offset: 4; indent-tabs-mode: nil -*-
#include "panel_topology.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include "main.h"
#include <sstream>
#include <cassert>
#include <cmath>
#include <set>

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
    std::vector<double> edge_gibps;
    std::vector<double> edge_latency_ns;
    struct RouteRecord {
        uint32_t source;
        uint32_t destination;
        std::vector<uint32_t> nodes;
    };
    std::vector<RouteRecord> route_records;
    long maxid = (long)_n - 1;
    std::string tok;
    while (f >> tok) {
        if (tok == "E") {
            long a, b; f >> a >> b;
            std::string rest;
            std::getline(f, rest);
            double edge_rate = 0.0;
            double edge_latency = -1.0;
            std::istringstream fields(rest);
            if (!(fields >> edge_rate)) edge_rate = 0.0;
            if (!(fields >> edge_latency)) edge_latency = -1.0;
            if (edge_rate < 0.0 || edge_latency < -1.0) {
                cerr << "PanelTopology: invalid custom edge attributes for "
                     << a << " -> " << b << endl;
                exit(1);
            }
            edges.push_back(std::make_pair(a, b));
            edge_gibps.push_back(edge_rate);
            edge_latency_ns.push_back(edge_latency);
            if (a > maxid) maxid = a;
            if (b > maxid) maxid = b;
        } else if (tok == "R") {
            long source, destination;
            f >> source >> destination;
            std::string rest;
            std::getline(f, rest);
            if (source < 0 || destination < 0 ||
                source >= (long)_n || destination >= (long)_n ||
                source == destination) {
                cerr << "PanelTopology: invalid custom route endpoints "
                     << source << " -> " << destination << endl;
                exit(1);
            }
            RouteRecord record;
            record.source = (uint32_t)source;
            record.destination = (uint32_t)destination;
            std::istringstream nodes(rest);
            long node;
            while (nodes >> node) {
                if (node < 0) {
                    cerr << "PanelTopology: negative custom route node" << endl;
                    exit(1);
                }
                record.nodes.push_back((uint32_t)node);
            }
            route_records.push_back(record);
        } else { std::string rest; std::getline(f, rest); }
    }
    _ndev = (uint32_t)(maxid + 1);
    _adj.assign(_ndev, std::vector<uint32_t>());
    std::vector<std::vector<double>> adjacency_gibps(_ndev);
    std::vector<std::vector<double>> adjacency_latency_ns(_ndev);
    size_t explicit_rates = 0;
    size_t explicit_latencies = 0;
    for (size_t i = 0; i < edges.size(); i++) {
        _adj[edges[i].first].push_back((uint32_t)edges[i].second);
        adjacency_gibps[edges[i].first].push_back(edge_gibps[i]);
        adjacency_latency_ns[edges[i].first].push_back(edge_latency_ns[i]);
        if (edge_gibps[i] > 0.0) explicit_rates++;
        if (edge_latency_ns[i] >= 0.0) explicit_latencies++;
    }
    cerr << "custom graph: " << edges.size() << " directed links, "
         << explicit_rates << " with explicit per-link GiB/s (rest at "
         << gibps << "), " << explicit_latencies
         << " with explicit per-link latency ns (rest at "
         << timeAsNs(lat) << ")" << endl;
    _dir_q.assign(_ndev, std::vector<LedgerQueue*>());
    _dir_p.assign(_ndev, std::vector<Pipe*>());
    for (uint32_t u = 0; u < _ndev; u++) {
        for (size_t p = 0; p < _adj[u].size(); p++) {
            char nm[64];
            snprintf(nm, sizeof(nm), "cq_%u_%zu", u, p);
            double edge_rate = adjacency_gibps[u][p] > 0.0
                ? adjacency_gibps[u][p] : gibps;
            _dir_q[u].push_back(make_queue(edge_rate, nm));
            snprintf(nm, sizeof(nm), "cp_%u_%zu", u, p);
            simtime_picosec edge_latency = adjacency_latency_ns[u][p] >= 0.0
                ? timeFromNs(adjacency_latency_ns[u][p]) : lat;
            _dir_p[u].push_back(make_pipe(edge_latency, nm));
        }
    }
    // Equal-cost shortest-path next-hop table.  Distances are computed from
    // each destination over the reverse graph, then every outgoing edge that
    // reduces the distance by one is retained.
    std::vector<std::vector<std::pair<uint32_t, int>>> reverse(_ndev);
    for (uint32_t u = 0; u < _ndev; u++) {
        for (size_t p = 0; p < _adj[u].size(); p++) {
            reverse[_adj[u][p]].push_back(std::make_pair(u, (int)p));
        }
    }
    _nh.assign(_ndev, std::vector<std::vector<int>>(_ndev));
    size_t candidate_sets = 0;
    size_t multipath_sets = 0;
    size_t maximum_width = 0;
    for (uint32_t d = 0; d < _ndev; d++) {
        std::vector<int> distance(_ndev, -1);
        std::vector<uint32_t> q(1, d);
        distance[d] = 0;
        for (size_t qi = 0; qi < q.size(); qi++) {
            uint32_t v = q[qi];
            for (const auto& predecessor : reverse[v]) {
                uint32_t u = predecessor.first;
                if (distance[u] < 0) {
                    distance[u] = distance[v] + 1;
                    q.push_back(u);
                }
            }
        }
        for (uint32_t u = 0; u < _ndev; u++) {
            if (distance[u] <= 0) continue;
            for (size_t p = 0; p < _adj[u].size(); p++) {
                uint32_t v = _adj[u][p];
                if (distance[v] == distance[u] - 1) {
                    _nh[u][d].push_back((int)p);
                }
            }
            if (!_nh[u][d].empty()) {
                candidate_sets++;
                maximum_width = std::max(maximum_width, _nh[u][d].size());
                if (_nh[u][d].size() > 1) multipath_sets++;
            }
        }
    }
    std::set<uint64_t> route_keys;
    for (const RouteRecord& record : route_records) {
        uint64_t key = ((uint64_t)record.source << 32) | record.destination;
        if (!route_keys.insert(key).second) {
            cerr << "PanelTopology: duplicate custom route override for "
                 << record.source << " -> " << record.destination << endl;
            exit(1);
        }
        if (record.nodes.size() < 2 ||
            record.nodes.front() != record.source ||
            record.nodes.back() != record.destination) {
            cerr << "PanelTopology: custom route endpoint mismatch for "
                 << record.source << " -> " << record.destination << endl;
            exit(1);
        }
        std::vector<int> ports;
        uint32_t current = record.source;
        for (size_t i = 1; i < record.nodes.size(); i++) {
            uint32_t next = record.nodes[i];
            if (current >= _adj.size() || next >= _adj.size()) {
                cerr << "PanelTopology: custom route node outside graph" << endl;
                exit(1);
            }
            std::vector<uint32_t>::const_iterator edge =
                std::find(_adj[current].cbegin(), _adj[current].cend(), next);
            if (edge == _adj[current].cend()) {
                cerr << "PanelTopology: custom route uses missing edge "
                     << current << " -> " << next << endl;
                exit(1);
            }
            int port = (int)std::distance(_adj[current].cbegin(), edge);
            const std::vector<int>& shortest = _nh[current][record.destination];
            if (std::find(shortest.begin(), shortest.end(), port) == shortest.end()) {
                cerr << "PanelTopology: custom route is not shortest for "
                     << record.source << " -> " << record.destination << endl;
                exit(1);
            }
            ports.push_back(port);
            current = next;
        }
        _route_overrides[key] = ports;
    }
    cout << "CUSTOM_GRAPH_ECMP"
         << " status=ENABLED"
         << " endpoints=" << _n
         << " devices=" << _ndev
         << " candidate_sets=" << candidate_sets
         << " multipath_sets=" << multipath_sets
         << " max_width=" << maximum_width
         << endl;
    cout << "CUSTOM_GRAPH_ROUTES"
         << " status=ENABLED"
         << " overrides=" << _route_overrides.size()
         << endl;
}

static uint64_t panel_pair_hash(uint32_t source, uint32_t destination,
                                uint32_t current) {
    uint64_t value = ((uint64_t)source << 32) ^ destination ^
                     ((uint64_t)current * 0x9e3779b97f4a7c15ULL);
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

PanelTopology::Candidate PanelTopology::direct_candidate(uint32_t src, uint32_t dest) {
    Candidate cand;
    cand.route = new Route();
    cand.is_plane = false; cand.plane = -1;
    cand.hops = 0; cand.latency_sum = 0;

    if (_base == Base::Custom) {
        uint32_t cur = src;
        uint64_t key = ((uint64_t)src << 32) | dest;
        std::unordered_map<uint64_t, std::vector<int>>::const_iterator override =
            _route_overrides.find(key);
        size_t override_index = 0;
        while (cur != dest) {
            const std::vector<int>& ports = _nh[cur][dest];
            assert(!ports.empty());
            int port;
            if (override != _route_overrides.end()) {
                assert(override_index < override->second.size());
                port = override->second[override_index++];
            } else {
                port = ports[panel_pair_hash(src, dest, cur) % ports.size()];
            }
            LedgerQueue* q = _dir_q[cur][port];
            Pipe* p = _dir_p[cur][port];
            cand.route->push_back(q); cand.route->push_back(p);
            cand.hop_queues.push_back(q);
            cand.latency_sum += p->delay();
            cand.hops++;
            cur = _adj[cur][port];
        }
        assert(override == _route_overrides.end() ||
               override_index == override->second.size());
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
