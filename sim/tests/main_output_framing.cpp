// Deterministic output-framing fixture for the shared HTSim progress clock.

#include <iostream>

#include "clock.h"
#include "eventlist.h"

int main() {
    EventList eventlist;
    eventlist.setEndtime(timeFromSec(0.21));
    Clock clock(timeFromSec(0.05), eventlist);

    while (eventlist.doNextEvent()) {
    }

    std::cout
        << "NETWORK_ROUTE class=switch hops=2 messages=1 "
           "payload_bytes=123456789 byte_hops=246913578 "
           "propagation_ns=0 serialization_ns=0\n"
        << "OCS_PLAN_REPLAY reconfigurations=1 rounds_advanced=1 "
           "scheduled_bytes=123456789 transmitted_bytes=123456789 "
           "reconf_ns=10\n"
        << "OCS_CONFIG_DRAIN plane=0 config=0 round=0 expected_stripes=1 "
           "started_stripes=1 completed_stripes=1 status=PASS\n"
        << "OCS_PLAN_AUDIT expected_flows=1 started_flows=1 "
           "completed_flows=1 expected_stripes=1 started_stripes=1 "
           "completed_stripes=1 expected_slots=1 consumed_slots=1 "
           "fallback_lookups=0 status=PASS\n"
        << "BACKEND_CAPABILITIES backend_contract=1\n";
    return 0;
}
