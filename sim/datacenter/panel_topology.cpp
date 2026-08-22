// -*- c-basic-offset: 4; indent-tabs-mode: nil -*-
#include "panel_topology.h"
#include "main.h"
#include <sstream>
#include <cassert>
#include <cmath>

using namespace std;

static string ntos(int v) { stringstream s; s << v; return s.str(); }

PanelTopology::PanelTopology(uint32_t npus, Base base, int planes,
                             double base_gibps, simtime_picosec base_latency,
                             double plane_gibps, simtime_picosec plane_latency,
                             mem_b queuesize, Logfile* logfile, EventList* ev)
    : _n(npus), _base(base), _planes(planes), _queuesize(queuesize),
      _logfile(logfile), _ev(ev) {
    switch (base) {
        case Base::None:    _dims = 0; _wrap = false; break;
        case Base::Mesh2D:  _dims = 2; _wrap = false; break;
        case Base::Torus2D: _dims = 2; _wrap = true;  break;
        case Base::Mesh3D:  _dims = 3; _wrap = false; break;
        case Base::Torus3D: _dims = 3; _wrap = true;  break;
    }
    if (_dims > 0) {
        _extent = (int)lround(pow((double)npus, 1.0 / _dims));
        int check = 1;
        for (int d = 0; d < _dims; d++) check *= _extent;
        if ((uint32_t)check != npus) {
            cerr << "PanelTopology: npus " << npus << " is not a perfect "
                 << _dims << "-dim grid" << endl;
            exit(1);
        }
        build_base(base_gibps, base_latency);
    }
    if (_planes > 0) build_planes(plane_gibps, plane_latency);
    assert(_dims > 0 || _planes > 0);
}

LedgerQueue* PanelTopology::make_queue(double gibps, const string& name) {
    QueueLoggerSampling* ql = new QueueLoggerSampling(timeFromUs((uint32_t)1000), *_ev);
    _logfile->addLogger(*ql);
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
    for (int d = 0; d < dim; d++) v /= _extent;
    return (int)(v % _extent);
}

uint32_t PanelTopology::id_of(const vector<int>& c) const {
    uint32_t id = 0;
    for (int d = _dims - 1; d >= 0; d--) id = id * _extent + c[d];
    return id;
}

void PanelTopology::build_base(double gibps, simtime_picosec lat) {
    // port index: 2*dim + 0 for +direction, 2*dim + 1 for -direction
    _dir_q.assign(_n, vector<LedgerQueue*>(2 * _dims, (LedgerQueue*)NULL));
    _dir_p.assign(_n, vector<Pipe*>(2 * _dims, (Pipe*)NULL));
    for (uint32_t node = 0; node < _n; node++) {
        for (int d = 0; d < _dims; d++) {
            int c = coord(node, d);
            bool plus_exists = (c + 1 < _extent) || _wrap;
            bool minus_exists = (c > 0) || _wrap;
            if (_extent <= 2) {
                // extent 2: +1 and -1 reach the same neighbour; wrap adds no
                // second edge. Only the + port exists (matches Mesh2D's
                // "wraparound && width > 2" guard in the analytical study).
                minus_exists = false;
                plus_exists = (c + 1 < _extent) || _wrap;
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

int PanelTopology::step_towards(int& cur, int target, bool tie_backward) const {
    int port;
    if (!_wrap) {
        if (target > cur) { cur++; port = 0; }
        else { cur--; port = 1; }
    } else {
        int fwd = ((target - cur) % _extent + _extent) % _extent;
        int bwd = _extent - fwd;
        bool go_fwd;
        if (fwd < bwd) go_fwd = true;
        else if (bwd < fwd) go_fwd = false;
        else go_fwd = !tie_backward;   // antipodal: split by source parity
        if (go_fwd) { cur = (cur + 1) % _extent; port = 0; }
        else { cur = (cur - 1 + _extent) % _extent; port = 1; }
    }
    return port;
}

PanelTopology::Candidate PanelTopology::direct_candidate(uint32_t src, uint32_t dest) {
    Candidate cand;
    cand.route = new Route();
    cand.is_plane = false; cand.plane = -1;
    cand.hops = 0; cand.latency_sum = 0;

    vector<int> c(_dims), t(_dims);
    for (int d = 0; d < _dims; d++) { c[d] = coord(src, d); t[d] = coord(dest, d); }

    // dimension-order X -> Y (-> Z); parity tie-break keyed on the SOURCE
    // coordinate in that dimension (mirrors Hybrid2D/Mesh3D).
    for (int d = 0; d < _dims; d++) {
        bool tie_backward = (coord(src, d) & 1) != 0;
        while (c[d] != t[d]) {
            uint32_t at = id_of(c);
            int port = step_towards(c[d], t[d], tie_backward);
            int pidx = 2 * d + port;
            if (_extent == 2 && _dir_q[at][pidx] == NULL) pidx = 2 * d;  // folded pair
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
    if (_dims > 0) out->push_back(direct_candidate(src, dest));
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
