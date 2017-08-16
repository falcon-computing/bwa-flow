#ifndef BWA_FLOW_MPICHANNEL_H
#define BWA_FLOW_MPICHANNEL_H
#include <boost/thread/mutex.hpp>
#include <boost/thread/lockable_adapter.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread.hpp>
#include <gtest/gtest_prod.h>
#include <unordered_set>

#include "mpi.h"

class MPILink
: public boost::basic_lockable_adapter<boost::mutex> {
 public:
  // message primitives
  bool query(MPI::Request &req);

  void send(MPI::Intercomm comm, 
      const void* buf, int count, 
      const MPI::Datatype& datatype, 
      int dest, int tag);

  void recv(MPI::Intercomm comm, 
      void* buf, int count, 
      const MPI::Datatype& datatype, 
      int source, int tag);
};

class Channel 
: public boost::basic_lockable_adapter<boost::mutex> {

 public:
  Channel(MPILink* link);

  virtual void retire() = 0; // retire from channel if is sender
  virtual void send(const char* data, int length) = 0;
  virtual void* recv(int & length) = 0;

  // separate two finish var because sometimes sender 
  // and receiver share the same channel object
  bool sendFinished() { return is_send_finished_; }
  bool recvFinished() { return is_recv_finished_; }
  
 protected:
  FRIEND_TEST(ChannelTests, ChannelSetup); 

  enum Msg {
    req    = 0,
    length = 1,
    data   = 2,
    MsgMax = data
  };

  int getTag(Msg m);
  int getId() { return id_; }

  // static channel id
  int id_;
  static int counter_;

  MPILink* link_;

  int rank_;
  int nproc_;

  bool is_send_finished_;
  bool is_recv_finished_;
};

class SourceChannel : public Channel {

 public:
  /*
   * - send_to_source: sometimes the sender process will also
   *   receive data, but sometimes it does not. In case it
   *   does not, the sender process need to be excluded from
   *   the receivers list, otherwise the recv() in retire will
   *   block
   */
  SourceChannel(MPILink* link, 
      int source_rank = 0, 
      bool send_to_source = true);

  void retire();
  void send(const char* data, int length);
  void* recv(int & length);
  
 private:
  int source_rank_;
  std::unordered_set<int> active_receiver_;
};

class SinkChannel : public Channel {

 public:
  /*
   * - recv_from_sink: similar to send_to_source,
   *   prevent recv() from hang
   */
  SinkChannel(MPILink* link, 
      int sink_rank = 0,
      bool recv_from_sink = true);

  void retire();
  void send(const char* data, int length);
  void* recv(int & length);
  
 private:
  int sink_rank_;
  std::unordered_set<int> active_senders_;
};


#endif
