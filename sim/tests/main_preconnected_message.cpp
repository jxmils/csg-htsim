// Verify exact-byte packetization for an already-established message transport.

#include <cstdint>
#include <iostream>

#include "eventlist.h"
#include "network.h"
#include "pipe.h"
#include "queue.h"
#include "tcp.h"

namespace {

int send_callbacks = 0;
int recv_callbacks = 0;
int send_bytes = 0;
int recv_bytes = 0;
simtime_picosec receive_time = 0;
EventList* test_eventlist = nullptr;

void finish_send(int, int, int bytes, int) {
    ++send_callbacks;
    send_bytes = bytes;
}

void finish_recv(int, int, int bytes, int) {
    ++recv_callbacks;
    recv_bytes = bytes;
    receive_time = test_eventlist->now();
}

}  // namespace

int main() {
    constexpr uint64_t message_bytes = 128;
    constexpr uint16_t maximum_packet_bytes = 8192;
    EventList eventlist;
    test_eventlist = &eventlist;
    eventlist.setEndtime(timeFromUs((uint32_t)10));
    Packet::set_packet_size(1500);

    Queue queue(speedFromGbps(400), 1 << 20, eventlist, nullptr);
    Pipe pipe(timeFromNs(80), eventlist);
    TcpSrc source(nullptr, nullptr, eventlist);
    TcpSink sink;
    Route outgoing;
    Route returning;
    outgoing.push_back(&queue);
    outgoing.push_back(&pipe);
    outgoing.push_back(&sink);
    returning.push_back(&source);

    source._debug_srcid = 0;
    source._debug_dstid = 1;
    source.astrasim_flow_finish_send_cb = finish_send;
    sink._debug_srcid = 0;
    sink._debug_dstid = 1;
    sink.astrasim_flow_finish_recv_cb = finish_recv;
    source.setFlowId(17);
    source.configure_preconnected_message(
        message_bytes, maximum_packet_bytes);
    source.set_cwnd(message_bytes + maximum_packet_bytes);
    source.set_ssthresh(message_bytes + maximum_packet_bytes);
    source.connect(outgoing, returning, sink, 0);
    sink.setFlowId(17);

    while (eventlist.doNextEvent()) {}

    const simtime_picosec expected = timeFromNs(80) +
        queue.serializationTime(message_bytes);
    const bool pass =
        source._packets_sent == message_bytes &&
        send_callbacks == 1 && recv_callbacks == 1 &&
        send_bytes == (int)message_bytes && recv_bytes == (int)message_bytes &&
        receive_time == expected;
    std::cout << "PRECONNECTED_MESSAGE_AUDIT"
              << " payload_bytes=" << message_bytes
              << " transmitted_bytes=" << source._packets_sent
              << " receive_ps=" << receive_time
              << " expected_ps=" << expected
              << " send_callbacks=" << send_callbacks
              << " recv_callbacks=" << recv_callbacks
              << " status=" << (pass ? "PASS" : "FAIL") << std::endl;
    return pass ? 0 : 1;
}
