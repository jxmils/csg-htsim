#include "panel_topology.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <unistd.h>

using std::string;
using std::vector;

namespace {

void require_ns(const char* label, simtime_picosec actual, double expected_ns) {
    const double actual_ns = timeAsNs(actual);
    if (actual_ns != expected_ns) {
        std::cerr << label << ": expected " << expected_ns << " ns, got "
                  << actual_ns << " ns" << std::endl;
        std::exit(1);
    }
    std::cout << "PASS " << label << "=" << actual_ns << "ns" << std::endl;
}

simtime_picosec candidate_latency(
    PanelTopology& topology, uint32_t src, uint32_t dst, bool plane) {
    vector<PanelTopology::Candidate>* candidates =
        topology.get_candidates(src, dst);
    simtime_picosec result = 0;
    bool found = false;
    for (size_t i = 0; i < candidates->size(); ++i) {
        if ((*candidates)[i].is_plane == plane && !found) {
            result = (*candidates)[i].latency_sum;
            found = true;
        }
        delete (*candidates)[i].route;
    }
    delete candidates;
    if (!found) {
        std::cerr << "missing " << (plane ? "plane" : "direct")
                  << " candidate" << std::endl;
        std::exit(1);
    }
    return result;
}

}  // namespace

int main() {
    EventList eventlist;
    const string prefix = "/tmp/htsim-panel-latency-" +
        std::to_string(static_cast<long long>(getpid()));
    const string graph_path = prefix + ".edges";
    Logfile logfile("/dev/null", eventlist);
    const mem_b queue_bytes = 1024 * 1024;

    PanelTopology one_hop(
        4, PanelTopology::Base::Mesh2D, 0,
        200.0, timeFromNs(25.0), 200.0, timeFromNs(25.0),
        queue_bytes, &logfile, &eventlist, vector<int>{2, 2}, true);
    require_ns("direct_one_hop",
               candidate_latency(one_hop, 0, 1, false), 25.0);

    PanelTopology four_hops(
        5, PanelTopology::Base::Mesh2D, 0,
        200.0, timeFromNs(25.0), 200.0, timeFromNs(25.0),
        queue_bytes, &logfile, &eventlist, vector<int>{5, 1}, true);
    require_ns("direct_four_hops",
               candidate_latency(four_hops, 0, 4, false), 100.0);

    PanelTopology ocs(
        2, PanelTopology::Base::None, 1,
        200.0, timeFromNs(25.0), 200.0,
        panelPlaneLegLatencyFromWholePathNs(55.0),
        queue_bytes, &logfile, &eventlist, vector<int>(), true);
    require_ns("ocs_whole_path",
               candidate_latency(ocs, 0, 1, true), 55.0);

    PanelTopology legacy_ocs(
        2, PanelTopology::Base::None, 1,
        200.0, timeFromNs(1000.0), 200.0, timeFromNs(1000.0),
        queue_bytes, &logfile, &eventlist, vector<int>(), true);
    require_ns("legacy_ocs_two_edges",
               candidate_latency(legacy_ocs, 0, 1, true), 2000.0);

    {
        std::ofstream graph(graph_path.c_str());
        graph << "E 0 2 200 75\n"
              << "E 2 1 200 325\n"
              << "E 1 2 200 75\n"
              << "E 2 0 200 325\n";
    }
    PanelTopology packet_switch(
        2, PanelTopology::Base::Custom, 0,
        200.0, timeFromNs(1000.0), 200.0, timeFromNs(1000.0),
        queue_bytes, &logfile, &eventlist, vector<int>(), true, graph_path);
    require_ns("packet_switch_whole_path",
               candidate_latency(packet_switch, 0, 1, false), 400.0);

    LedgerQueue serialization_queue(
        speedFromGiBps(200.0), queue_bytes, eventlist, NULL);
    require_ns("serialization_1500B_at_200GiBps",
               serialization_queue.serializationTime(1500), 6.985);

    std::remove(graph_path.c_str());
    std::cout << "PANEL_LATENCY_TEST status=PASS" << std::endl;
    return 0;
}
