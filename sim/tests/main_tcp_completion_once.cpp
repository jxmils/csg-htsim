// Verify that duplicate terminal packets cannot emit duplicate Astra callbacks.

#include <iostream>
#include <memory>
#include <set>
#include <vector>

#include "eventlist.h"
#include "network.h"
#include "tcp.h"
#include "tcppacket.h"

namespace {

int send_callbacks = 0;
int recv_callbacks = 0;
std::set<int> send_tags;
std::set<int> recv_tags;

void finish_send(int, int, int, int tag) {
    ++send_callbacks;
    send_tags.insert(tag);
}

void finish_recv(int, int, int, int tag) {
    ++recv_callbacks;
    recv_tags.insert(tag);
}

}  // namespace

int main() {
    constexpr int flow_count = 3;
    EventList eventlist;
    eventlist.setEndtime(timeFromSec(2));
    Packet::set_packet_size(1000);

    std::vector<std::unique_ptr<TcpSrc>> sources;
    std::vector<std::unique_ptr<TcpSink>> sinks;
    std::vector<std::unique_ptr<Route>> outgoing;
    std::vector<std::unique_ptr<Route>> returning;
    std::vector<std::unique_ptr<PacketFlow>> packet_flows;

    for (int index = 0; index < flow_count; ++index) {
        const int runtime_flow_id = 5121 + index;
        sources.emplace_back(new TcpSrc(nullptr, nullptr, eventlist));
        sinks.emplace_back(new TcpSink());
        outgoing.emplace_back(new Route());
        returning.emplace_back(new Route());
        packet_flows.emplace_back(new PacketFlow(nullptr));

        TcpSrc& source = *sources.back();
        TcpSink& sink = *sinks.back();
        source.setFlowId(runtime_flow_id);
        source.set_flowsize(0);
        source._established = true;
        source._highest_sent = source._flow_size;
        source._last_acked = source._flow_size;
        source._debug_srcid = 32;
        source._debug_dstid = 33;
        source.astrasim_flow_finish_send_cb = finish_send;
        sink._debug_srcid = 32;
        sink._debug_dstid = 33;
        sink.astrasim_flow_finish_recv_cb = finish_recv;
        outgoing.back()->push_back(&sink);
        returning.back()->push_back(&source);
        source.connect(*outgoing.back(), *returning.back(), sink,
                       timeFromSec(1));
        packet_flows.back()->set_flowid(runtime_flow_id);

        // The second packet models a retransmitted terminal segment. It
        // generates another cumulative final ACK, but neither terminal
        // callback may be emitted a second time.
        for (int copy = 0; copy < 2; ++copy) {
            TcpPacket* packet = TcpPacket::newpkt(
                *packet_flows.back(), *outgoing.back(), 1,
                static_cast<int>(source._flow_size));
            packet->set_ts(0);
            sink.receivePacket(*packet);
        }
    }

    const bool pass =
        send_callbacks == flow_count && recv_callbacks == flow_count &&
        send_tags.size() == flow_count && recv_tags.size() == flow_count;
    std::cout << "TCP_COMPLETION_AUDIT"
              << " runtime_flows=" << flow_count
              << " terminal_packets=" << 2 * flow_count
              << " send_callbacks=" << send_callbacks
              << " recv_callbacks=" << recv_callbacks
              << " unique_send_tags=" << send_tags.size()
              << " unique_recv_tags=" << recv_tags.size()
              << " status=" << (pass ? "PASS" : "FAIL") << std::endl;
    return pass ? 0 : 1;
}
