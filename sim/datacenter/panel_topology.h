// -*- c-basic-offset: 4; indent-tabs-mode: nil -*-
// PanelTopology: direct-connect grid fabrics (2D/3D mesh and torus) with an
// optional set of switch planes attached to every node -- the packet-level
// counterpart of the analytical study's Mesh3D / Torus3D / Hybrid2D /
// MultiPlaneSwitch classes. One class covers all four fixed-port-budget
// configurations:
//
//   Torus3D      base=Torus3D, planes=0      (6D+0F)
//   Mesh3D       base=Mesh3D,  planes=0      (avg 4.5 ports, sparse control)
//   Hybrid       base=Torus2D, planes=2      (4D+2F)
//   FullSwitch   base=None,    planes=6      (0D+6F)
//
// Routing semantics deliberately mirror the analytical implementation
// (Hybrid2D.cpp / Mesh3D.cpp): dimension-order X->Y(->Z) with a
// source-coordinate-parity tie-break on equidistant torus wraps; switch-plane
// routes are exactly two hops (uplink, downlink). Route selection among the
// candidates does NOT happen here -- get_bidir_paths returns the full
// candidate set (direct first, then plane 0..P-1) and the ASTRA frontend
// applies the routing policy at flow-injection time.
#ifndef PANEL_TOPOLOGY_H
#define PANEL_TOPOLOGY_H

#include "topology.h"
#include "queue.h"
#include "pipe.h"
#include "eventlist.h"
#include "logfile.h"
#include "loggers.h"
#include <vector>

// Bandwidth convention: the analytical backend treats "GB/s" as GiB/s
// (bw_GBps_to_Bpns multiplies by 2^30/1e9). To reproduce its operating points
// exactly, panel links are specified in the same unit.
inline linkspeed_bps speedFromGiBps(double gibps) {
    return (linkspeed_bps)(gibps * 8.0 * 1073741824.0);
}

// A FIFO queue that additionally tracks bytes RESERVED by flows routed through
// it but not yet begun serializing -- the packet-level mirror of the
// analytical Link::outstanding_bytes (reserve at injection, release when
// serialization starts). The Adaptive cost term reads reserved_bytes().
class LedgerQueue : public Queue {
  public:
    LedgerQueue(linkspeed_bps bitrate, mem_b maxsize, EventList& eventlist,
                QueueLogger* logger)
        : Queue(bitrate, maxsize, eventlist, logger), _reserved(0) {}
    void reserve_bytes(uint64_t b) { _reserved += b; }
    uint64_t reserved_bytes() const { return _reserved; }
    linkspeed_bps link_bitrate() const { return _bitrate; }
    virtual void beginService() {
        if (!_enqueued.empty()) {
            uint64_t sz = (uint64_t)_enqueued.back()->size();
            _reserved = (_reserved > sz) ? (_reserved - sz) : 0;
        }
        Queue::beginService();
    }
  private:
    uint64_t _reserved;
};

class PanelTopology : public Topology {
  public:
    enum class Base { None, Ring1D, Mesh2D, Torus2D, Mesh3D, Torus3D, Custom, RingRows };

    // One route candidate, with the metadata the frontend's cost function and
    // telemetry need. hop_queues are the LedgerQueues along the path in order.
    struct Candidate {
        Route* route;                     // ends at last pipe; caller appends sink
        bool is_plane;                    // false = direct fabric
        int plane;                        // valid when is_plane
        int hops;                         // link count
        simtime_picosec latency_sum;      // sum of pipe latencies
        std::vector<LedgerQueue*> hop_queues;
    };

    // extents: per-dimension sizes (e.g. {4,8,8}); empty = uniform grid.
    PanelTopology(uint32_t npus, Base base, int planes,
                  double base_gibps, simtime_picosec base_latency,
                  double plane_gibps, simtime_picosec plane_latency,
                  mem_b queuesize, Logfile* logfile, EventList* ev,
                  const std::vector<int>& extents = std::vector<int>(),
                  bool quiet = false,
                  const std::string& graphfile = "");

    // Topology interface. Returns the candidate Routes (copies of route
    // structure, same queue objects): direct first (if a base fabric exists),
    // then one per plane.
    virtual vector<const Route*>* get_bidir_paths(uint32_t src, uint32_t dest,
                                                  bool reverse);
    virtual vector<uint32_t>* get_neighbours(uint32_t) { return NULL; }
    virtual uint32_t no_of_nodes() const { return _n; }

    // Rich candidate set for the frontend router. Caller owns the vector and
    // the Route objects inside (queues/pipes are topology-owned).
    std::vector<Candidate>* get_candidates(uint32_t src, uint32_t dest);

    int planes() const { return _planes; }
    // Custom base only: -ecmp selects among all equal-cost next hops by
    // dest mod k at each hop (merlin fat-tree deterministic rule). Off = old
    // single-path BFS table, unchanged.
    static bool ecmp;
    bool has_base() const { return _base != Base::None; }

  private:
    uint32_t _n;          // endpoints (ranks visible to ASTRA)
    uint32_t _ndev = 0;   // physical devices incl. switches (Custom)
    Base _base;
    int _dims;            // 2 or 3 (0 when Base::None)
    std::vector<int> _extents;   // per-dimension extents
    bool _wrap;
    int _planes;
    mem_b _queuesize;
    Logfile* _logfile;
    EventList* _ev;
    bool _quiet;   // suppress per-queue sampling loggers (large-N runs)

    // Direct fabric: per node, per port (2*dim + sign; sign 0 = +, 1 = -).
    std::vector<std::vector<LedgerQueue*>> _dir_q;
    std::vector<std::vector<Pipe*>> _dir_p;
    // Planes: [plane][node] uplink and downlink.
    std::vector<std::vector<LedgerQueue*>> _up_q;
    std::vector<std::vector<Pipe*>> _up_p;
    std::vector<std::vector<LedgerQueue*>> _down_q;
    std::vector<std::vector<Pipe*>> _down_p;

    // Custom base: arbitrary directed graph from file ("E src dst" lines),
    // BFS shortest-path next-hop routing.
    std::vector<std::vector<uint32_t>> _adj;      // [node] -> out-neighbors
    std::vector<std::vector<int>> _nh;            // [src][dst] -> out-port (-1 none)
    std::vector<std::vector<std::vector<int>>> _nhs; // [src][dst] -> all shortest-path out-ports (-ecmp)
    std::string _graphfile;
    void build_custom(double gibps, simtime_picosec lat);

    int coord(uint32_t id, int dim) const;
    uint32_t id_of(const std::vector<int>& c) const;
    // one dimension-order step; returns port index used and updates cur
    int step_towards(int& cur, int target, bool tie_backward, int extent) const;
    LedgerQueue* make_queue(double gibps, const std::string& name);
    Pipe* make_pipe(simtime_picosec lat, const std::string& name);
    void build_base(double gibps, simtime_picosec lat);
    void build_planes(double gibps, simtime_picosec lat);
    Candidate direct_candidate(uint32_t src, uint32_t dest);
    Candidate plane_candidate(uint32_t src, uint32_t dest, int plane);
};

#endif
