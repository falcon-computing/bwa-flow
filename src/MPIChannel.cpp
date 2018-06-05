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
#include "allocation_wrapper.h"

#include "mpi.h"
#include "MPIChannel.h"

bool MPILink::query(MPI::Request &req) {
  boost::lock_guard<MPILink> guard(*this);
  return req.Test();
}

void MPILink::send(MPI::Intercomm comm,
    const void* buf,
    int count, const MPI::Datatype& datatype,
    int dest, int tag
) {
  MPI::Request req;
  {
    boost::lock_guard<MPILink> guard(*this);
    req = comm.Isend(buf, count,
        datatype, dest, tag);
    DLOG_IF(INFO, VLOG_IS_ON(4)) << "Send to " << dest 
               << " from " << MPI::COMM_WORLD.Get_rank()
               << " tag: " << tag;
  }
  while (!query(req)) {
    boost::this_thread::sleep_for(boost::chrono::milliseconds(1));
  }
}

void MPILink::recv(MPI::Intercomm comm, void* buf,
    int count, const MPI::Datatype& datatype,
    int source, int tag
) {
  MPI::Request req;
  {
    boost::lock_guard<MPILink> guard(*this);
    req = comm.Irecv(buf, count,
        datatype, source, tag);

    DLOG_IF(INFO, VLOG_IS_ON(4)) << "Recv from " << source 
               << " to " << MPI::COMM_WORLD.Get_rank()
               << " tag: " << tag;
  }
  while (!query(req)) {
    boost::this_thread::sleep_for(boost::chrono::milliseconds(1));
  }
}

int Channel::counter_ = 0;

Channel::Channel(MPILink* link): 
  link_(link),
  is_send_finished_(false),
  is_recv_finished_(false) 
{
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

int Channel::getTag(Msg m) {
  int tag = ((int)Msg::MsgMax+1)*id_ + (int)m;
  return tag;
}

SourceChannel::SourceChannel(MPILink* link, 
    int source_rank,
    bool send_to_source): 
  Channel(link), 
  source_rank_(source_rank) 
{
  for (int p = 0; p < nproc_; p++) {
    if (send_to_source || p != source_rank) {
      active_receiver_.insert(p);
    }
  }
}

// dispatch data to receivers of this channel
void SourceChannel::send(const char* data, int length) {
  // pid of the requester
  int req_id = -1;

  link_->recv(MPI::COMM_WORLD, &req_id, 1, 
      MPI::INT, MPI::ANY_SOURCE, getTag(Msg::req));
  link_->send(MPI::COMM_WORLD, &length, 1, 
      MPI::INT, req_id, getTag(Msg::length));
  link_->send(MPI::COMM_WORLD, data, length, 
      MPI::CHAR, req_id, getTag(Msg::data));
}

void SourceChannel::retire() {
  if (rank_ != source_rank_) {
    throw std::runtime_error("unexpected caller of finish()");
  }

  for (auto p : active_receiver_) {
  //for (int p = 0; p < nproc_; p++) {
    int proc_id = 0;
    int length = 0;
    link_->recv(MPI::COMM_WORLD, &proc_id, 1, 
        MPI::INT, MPI::ANY_SOURCE, getTag(Msg::req));

    // Send finish signal to child-process
    link_->send(MPI::COMM_WORLD, &length, 1, 
        MPI::INT, p, getTag(Msg::length));
  }
  is_send_finished_ = true;
}

void* SourceChannel::recv(int & length) {

  link_->send(MPI::COMM_WORLD, &rank_, 1, 
      MPI::INT, source_rank_, getTag(Msg::req));

  //DLOG_IF(INFO, VLOG_IS_ON(2)) << __func__ << ":" << "receiving from " << source_rank_;

  link_->recv(MPI::COMM_WORLD, &length, 1, 
      MPI::INT, source_rank_, getTag(Msg::length));

  if (length > 0) {
    // allocate return buffer for data
    // allow one extra 0 for c strings
    void* data = calloc(length+1, 1);

    link_->recv(MPI::COMM_WORLD, data, length, 
        MPI::CHAR, source_rank_, getTag(Msg::data));

    return data;
  }
  else {
    is_recv_finished_ = true;
    return NULL;
  }
}

SinkChannel::SinkChannel(MPILink* link, 
    int sink_rank,
    bool recv_from_sink): 
  Channel(link),
  sink_rank_(sink_rank) 
{
  for (int p = 0; p < nproc_; p++) {
    if (recv_from_sink || p != sink_rank) {
      active_senders_.insert(p);
    }
  } 
}

void SinkChannel::send(const char* data, int length) {
  link_->send(MPI::COMM_WORLD, &rank_, 1, 
      MPI::INT, sink_rank_, getTag(Msg::req));
  DLOG_IF(INFO, VLOG_IS_ON(2)) << __func__ << ":" << "sending to " << sink_rank_;
  link_->send(MPI::COMM_WORLD, &length, 1, 
      MPI::INT, sink_rank_, getTag(Msg::length));
  link_->send(MPI::COMM_WORLD, data, length, 
      MPI::CHAR, sink_rank_, getTag(Msg::data));
}

void* SinkChannel::recv(int & length) {
  if (active_senders_.empty()) return NULL;
  
  // non-blocking query for tasks
  int req_id = 0;
  link_->recv(MPI::COMM_WORLD, &req_id, 1, MPI::INT, MPI::ANY_SOURCE, getTag(Msg::req));
  //DLOG_IF(INFO, VLOG_IS_ON(2)) << __func__ << ":" << "receiving from " << req_id;

  link_->recv(MPI::COMM_WORLD, &length, 1, MPI::INT, req_id, getTag(Msg::length));
  while (length == 0) {
    active_senders_.erase(req_id); 
    if (active_senders_.empty()) break;
    link_->recv(MPI::COMM_WORLD, &req_id, 1, MPI::INT, MPI::ANY_SOURCE, getTag(Msg::req));
    link_->recv(MPI::COMM_WORLD, &length, 1, MPI::INT, req_id, getTag(Msg::length));
  }
  if (active_senders_.empty()) {
    is_recv_finished_ = true;
    return NULL;
  }
  else {
    // allocate return buffer for data
    // allow one extra 0 for c strings
    void* data = calloc(length+1, 1);

    link_->recv(MPI::COMM_WORLD, data, length, MPI::CHAR, req_id, getTag(Msg::data));

    return data;
  }
}

void SinkChannel::retire() {
  //DLOG(INFO) << "Retiring sender " << rank_;
  int length = 0;
  link_->send(MPI::COMM_WORLD, &rank_, 1, MPI::INT, sink_rank_, getTag(Msg::req));
  link_->send(MPI::COMM_WORLD, &length, 1, MPI::INT, sink_rank_, getTag(Msg::length));

  is_send_finished_ = true;
}
