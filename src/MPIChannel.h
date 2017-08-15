#ifndef BWA_FLOW_MPICHANNEL_H
#define BWA_FLOW_MPICHANNEL_H
#include <boost/thread/mutex.hpp>
#include <boost/thread/lockable_adapter.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread.hpp>
#include <unordered_set>

#include "mpi.h"

class Channel 
: public boost::basic_lockable_adapter<boost::mutex> {

 public:
  Channel();

  //void retire(); // retire from channel if is sender
  virtual void finish() = 0;
  virtual void send(const char* data, int length) = 0;
  virtual void recv(char* data, int & length) = 0;
  //bool isFinish(); // if all senders of this channel have retired

  int getId() { return id_; }
  bool isFinished() { return is_finished_; }
  
 protected:
  // message primitives
  bool query(MPI::Request &req);

  void do_send(MPI::Intercomm comm, 
      const void* buf, int count, 
      const MPI::Datatype& datatype, 
      int dest, int tag);

  void do_recv(MPI::Intercomm comm, 
      void* buf, int count, 
      const MPI::Datatype& datatype, 
      int source, int tag);

  enum Msg {
    req    = 0,
    length = 1,
    data   = 2
  };
  int getTag(Msg m);

  // static channel id
  int id_;
  static int counter_;

  int rank_;
  int nproc_;

  bool is_finished_;
};

class SourceChannel : public Channel {

 public:
  SourceChannel(int source_rank = 0);

  void finish();
  void send(const char* data, int length);
  void recv(char* data, int & length);
  
 private:
  int source_rank_;
};

class SinkChannel : public Channel {

 public:
  SinkChannel(int sink_rank = 0);

  void finish();
  void send(const char* data, int length);
  void recv(char* data, int & length);
  
 private:
  int sink_rank_;
  std::unordered_set<int> active_senders_;
};


#endif
