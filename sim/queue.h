// -*- c-basic-offset: 4; indent-tabs-mode: nil -*-        
#ifndef QUEUE_H
#define QUEUE_H

/*
 * A simple FIFO queue
 */

//#include <list>
//#include "circular_buffer.h"
#include "config.h"
#include "eventlist.h"
#include "network.h"
#include "loggertypes.h"
#include "fairpullqueue.h"
#include "drawable.h"
#include "switch.h"
#include "circular_buffer.h"

// BaseQueue is a generic queue, but doesn't actually implement any
// queuing discipline.  Subclasses implement different queuing
// disciplines. 
class Switch;

class BaseQueue  : public EventSource, public PacketSink, public Drawable {
 public:
    BaseQueue(linkspeed_bps bitrate, EventList &eventlist, QueueLogger* logger);
    virtual void setLogger(QueueLogger* logger) {
            _logger = logger;
    }
    virtual void setName(const string& name) {
        Logged::setName(name); 
        _nodename += name;
    }
    virtual void forceName(const string& name) {
        Logged::setName(name); 
        _nodename = name;
    }

    virtual void setSwitch(Switch* s){assert(!_switch);_switch = s;}
    virtual Switch* getSwitch(){return _switch;}
    
    void setNext(PacketSink* next_sink) {
            _next_sink = next_sink;
    }
    PacketSink* next() const {
            return _next_sink;
    }
    virtual const string& nodename() { return _nodename; }
    virtual mem_b queuesize() const = 0;
    virtual mem_b maxsize() const = 0;
    
    // Quantize ONCE, for the whole packet, rounding up. Rounding per byte
    // first is what produced physically impossible link rates.
    inline simtime_picosec serializationTime(mem_b bytes) const {
        if (_bitrate == 0 || bytes <= 0) return 0;
        // 128-bit product: bytes * 8e12 overflows uint64 above ~2.3 MiB.
        // Packets are far smaller, but the whole reason this function exists
        // is that a silent arithmetic assumption produced impossible link
        // rates, so do not leave another one in place.
        const unsigned __int128 num =
            (unsigned __int128)(uint64_t)bytes * 8ULL * 1000000000000ULL;
        return (simtime_picosec)((num + _bitrate - 1) / _bitrate);
    }

    inline simtime_picosec drainTime(Packet *pkt) {
            return serializationTime((mem_b)pkt->size());
    }

    inline mem_b serviceCapacity(simtime_picosec t) { 
            return (mem_b)(timeAsSec(t) * (double)_bitrate); 
    }

    virtual void log_packet_send(simtime_picosec duration);
    virtual uint16_t average_utilization();

    virtual uint64_t quantized_queuesize();
    virtual uint8_t quantized_utilization();

    static simtime_picosec _update_period;

protected:
    // Housekeeping
    PacketSink* _next_sink; // used in generic topology for linkage
    QueueLogger* _logger;
    linkspeed_bps _bitrate; 
    // Kept as a double: rounding service time per BYTE loses up to 100% of
    // it (zero above ~931.32 GiB/s). Packet serialization is computed
    // exactly by serializationTime(); this value is only used for the
    // queueing-delay estimators below.
    double _ps_per_byte;  // service time, in picoseconds per byte
    string _nodename;
    
    CircularBuffer<simtime_picosec> _busystart;
    CircularBuffer<simtime_picosec> _busyend;

    //how much time have we spent being busy in the current measurement window?
    simtime_picosec _busy;
    simtime_picosec _idle;
    simtime_picosec _window;

    simtime_picosec _last_update_qs, _last_update_utilization;
    uint8_t _last_qs, _last_utilization;

    Switch* _switch;//which switch is this queue part of?
};


// A standard FIFO packet queue of fixed size
class Queue : public BaseQueue {
 public:
    Queue(linkspeed_bps bitrate, mem_b maxsize, EventList &eventlist, 
          QueueLogger* logger);
    virtual void receivePacket(Packet& pkt);
    void doNextEvent();
    // should really be private, but loggers want to see
    mem_b _maxsize; 

    virtual mem_b queuesize() const;
    virtual mem_b maxsize() const {return _maxsize;}
    simtime_picosec serviceTime();
    int num_drops() const {return _num_drops;}
    void reset_drops() {_num_drops = 0;}

 protected:
    // Mechanism
    // start serving the item at the head of the queue
    virtual void beginService(); 

    // wrap up serving the item at the head of the queue
    virtual void completeService(); 

    mem_b _queuesize;
    CircularBuffer<Packet*> _enqueued;
    int _num_drops;
};

class HostQueue : public Queue {
    public: 
        HostQueue(linkspeed_bps bitrate, mem_b maxsize, EventList &eventlist,  QueueLogger* logger);

        void addHostSender(PacketSink* snk) {_senders.push_back(snk);};
        vector<PacketSink*> _senders;
        virtual simtime_picosec serviceTime(Packet& pkt) = 0;

};

/* implement a 3-level priority queue */
class PriorityQueue : public HostQueue {
 public:
    typedef enum {Q_LO=0, Q_MID=1, Q_HI=2, Q_NONE=3} queue_priority_t;
    PriorityQueue(linkspeed_bps bitrate, mem_b maxsize, EventList &eventlist, 
                  QueueLogger* logger);
    virtual void receivePacket(Packet& pkt);
    virtual mem_b queuesize() const;
    virtual simtime_picosec serviceTime(Packet& pkt);

 protected:
    //this is needed for lossless operation!

    // start serving the item at the head of the queue
    virtual void beginService(); 

    // wrap up serving the item at the head of the queue
    virtual void completeService(); 
    PriorityQueue::queue_priority_t getPriority(Packet& pkt);
    list <Packet*> _queue[Q_NONE];
    mem_b _queuesize[Q_NONE];
    queue_priority_t _servicing;
    int _state_send;
};

/* implement a 3-level priority queue */

class FairPriorityQueue : public HostQueue {
 public:
    typedef enum {Q_LO=0, Q_MID=1, Q_HI=2, Q_NONE=3} queue_priority_t;
    FairPriorityQueue(linkspeed_bps bitrate, mem_b maxsize, EventList &eventlist, 
                  QueueLogger* logger);
    virtual void receivePacket(Packet& pkt);
    virtual mem_b queuesize() const;
    virtual simtime_picosec serviceTime(Packet& pkt);

 protected:
    //this is needed for lossless operation!
    

    // start serving the item at the head of the queue
    virtual void beginService(); 

    // wrap up serving the item at the head of the queue
    virtual void completeService(); 
    FairPriorityQueue::queue_priority_t getPriority(Packet& pkt);
    FairPullQueue<Packet> _queue[Q_NONE];

    Packet* _sending;
    mem_b _queuesize[Q_NONE];
    queue_priority_t _servicing;
    int _state_send;
};


#endif
