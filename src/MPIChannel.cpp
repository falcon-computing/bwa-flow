#include <boost/asio.hpp>
#include <boost/function_types/result_type.hpp>
#include <boost/make_shared.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread.hpp>
#include <fstream>
#include <glog/logging.h>
#include <list>
#include <queue>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "mpi.h"
#include "MPIChannel.h"

int Channel::counter_ = 0;

int Channel::getTag(Msg m) {
  int tag = std::numeric_limits<Msg>::max()*id_ + (int)m;
  return tag;
}

bool Channel::query(MPI::Request &req) {
  boost::lock_guard<Channel> guard(*this);
  return req.Test();
}

void Channel::do_send(MPI::Intercomm comm,
    const void* buf,
    int count, const MPI::Datatype& datatype,
    int dest, int tag
) {
  MPI::Request req;
  {
    boost::lock_guard<Channel> guard(*this);
    req = comm.Isend(buf, count,
        datatype, dest, tag);
  }
  while (!query(req)) {
    boost::this_thread::sleep_for(boost::chrono::milliseconds(1));
  }
}

void Channel::do_recv(MPI::Intercomm comm, void* buf,
    int count, const MPI::Datatype& datatype,
    int source, int tag
) {
  MPI::Request req;
  {
    boost::lock_guard<Channel> guard(*this);
    req = comm.Irecv(buf, count,
        datatype, source, tag);
  }
  while (!query(req)) {
    boost::this_thread::sleep_for(boost::chrono::milliseconds(1));
  }
}

Channel::Channel(): is_finished_(false) {
  
  // accumulate the class id
  id_ = counter_;
  counter_++;

  rank_  = MPI::COMM_WORLD.Get_rank();
  nproc_ = MPI::COMM_WORLD.Get_size();

#if 0 // removed and re-enable later
  /* Register/broadcast senders:
   * all procs enter, and communicate who are the senders of this Channel
   */
  std::vector<int> sender_list;
  std::vector<int> receiver_list;
  
  for (int p = 0; p < nproc_; p++) {
    // no matter if p == rank_ at this point, cause is_p_sender
    // would be overwritten if p != rank
    int is_p_sender = is_sender ? 1 : 0;

    MPI::COMM_WORLD.Bcast(&is_p_sender, 1, MPI::INT, p);

    if (is_p_sender) {
      sender_list.push_back(p);
    }
    else {
      receiver_list.push_back(p);
    }
  }

  // create mpi group and communicators 
  MPI::Group world_group = MPI::COMM_WORLD.Get_group();
  MPI::Group send_group = world_group.Incl(sender_list.size(), &sender_list[0]);
  MPI::Group recv_group = world_group.Incl(receiver_list.size(), &receiver_list[0]);

  send_comm_ = MPI::COMM_WORLD.Create(send_group);
  recv_comm_ = MPI::COMM_WORLD.Create(recv_group);

  num_send_ = send_comm_.Get_size();
  num_recv_ = recv_comm_.Get_size();

  if (is_sender) {
    g_rank_  = send_comm_.Get_rank();
  }
  else {
    g_rank_  = recv_comm_.Get_rank();
  }
#endif
}

SourceChannel::SourceChannel(int source_rank): source_rank_(source_rank) {
  ;
}

// dispatch data to receivers of this channel
void SourceChannel::send(const char* data, int length) {
  // pid of the requester
  int req_id = -1;

  do_recv(MPI::COMM_WORLD, &req_id, 1, 
      MPI::INT, MPI::ANY_SOURCE, getTag(Msg::req));
  do_send(MPI::COMM_WORLD, &length, 1, 
      MPI::INT, req_id, getTag(Msg::length));
  do_send(MPI::COMM_WORLD, data, length, 
      MPI::CHAR, req_id, getTag(Msg::data));
}

void SourceChannel::finish() {
  if (rank_ != source_rank_) {
    throw std::runtime_error("unexpected caller of finish()");
  }

  for (int p = 0; p < nproc_; p++) {
    int proc_id = 0;
    int length = 0;
    do_recv(MPI::COMM_WORLD, &proc_id, 1, 
        MPI::INT, MPI::ANY_SOURCE, getTag(Msg::req));

    // Send finish signal to child-process
    do_send(MPI::COMM_WORLD, &length, 1, 
        MPI::INT, p, getTag(Msg::length));
  }
  //is_finished_ = true;
}

void SourceChannel::recv(char* data, int & length) {

  do_send(MPI::COMM_WORLD, &rank_, 1, 
      MPI::INT, source_rank_, getTag(Msg::req));
  do_recv(MPI::COMM_WORLD, &length, 1, 
      MPI::INT, source_rank_, getTag(Msg::length));
  if (length > 0) {
    do_recv(MPI::COMM_WORLD, data, length, 
        MPI::CHAR, source_rank_, getTag(Msg::data));
  }
  else {
    is_finished_ = true;
  }
}

SinkChannel::SinkChannel(int sink_rank): sink_rank_(sink_rank) {
  for (int p = 0; p < nproc_; p++) {
    active_senders_.insert(p);
  } 
}

void SinkChannel::send(const char* data, int length) {
  do_send(MPI::COMM_WORLD, &rank_, 1, 
      MPI::INT, sink_rank_, getTag(Msg::req));
  do_send(MPI::COMM_WORLD, &length, 1, 
      MPI::INT, sink_rank_, getTag(Msg::length));
  do_send(MPI::COMM_WORLD, data, length, 
      MPI::CHAR, sink_rank_, getTag(Msg::data));
}

void SinkChannel::recv(char* data, int & length) {
  if (active_senders_.empty()) return;
  
  // non-blocking query for tasks
  int req_id = 0;
  do_recv(MPI::COMM_WORLD, &req_id, 1, MPI::INT, MPI::ANY_SOURCE, getTag(Msg::req));
  do_recv(MPI::COMM_WORLD, &length, 1, MPI::INT, req_id, getTag(Msg::length));
  while (length == 0) {
    active_senders_.erase(req_id); 
    if (active_senders_.empty()) break;
    do_recv(MPI::COMM_WORLD, &req_id, 1, MPI::INT, MPI::ANY_SOURCE, getTag(Msg::req));
    do_recv(MPI::COMM_WORLD, &length, 1, MPI::INT, req_id, getTag(Msg::length));
  }
  if (active_senders_.empty()) {
    is_finished_ = true;
  }
  else {
    do_recv(MPI::COMM_WORLD, data, length, MPI::CHAR, req_id, getTag(Msg::data));
  }
}

void SinkChannel::finish() {
  int length = 0;
  do_send(MPI::COMM_WORLD, &rank_, 1, MPI::INT, sink_rank_, getTag(Msg::req));
  do_send(MPI::COMM_WORLD, &length, 1, MPI::INT, sink_rank_, getTag(Msg::length));
  // 
  // is_finished_ = true;
}
