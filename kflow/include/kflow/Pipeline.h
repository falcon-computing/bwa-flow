#ifndef KFLOW_PIPELINE_H
#define KFLOW_PIPELINE_H

#include <boost/any.hpp>
#include <boost/asio.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
#include <deque>

#include "Common.h"
#include "Stage.h"

// for testing purpose
#ifndef TEST_FRIENDS_LIST
#define TEST_FRIENDS_LIST
#endif

namespace kestrelFlow {

class Pipeline {
  friend class OccupancyCounter;
  TEST_FRIENDS_LIST

 public:
  Pipeline(int num_stages, int num_threads);
  ~Pipeline();

  template <typename U, typename V,
            int IN_DEPTH,  int OUT_DEPTH> 
  bool addStage(int idx, Stage<U, V, IN_DEPTH, OUT_DEPTH> *stage);

  // branch stage[idx] of pipeline
  void branch(Pipeline& pipeline, int idx); 

  template <typename T>
  bool addConst(std::string key, T val);

  boost::any getConst(std::string key);

  void start();
  void stop();
  void wait();
  void finalize();

  // Post work to pipeline workers
  template <typename CompletionHandler>
    void post(CompletionHandler handler); 

  // Experimental feature 
  boost::shared_ptr<QueueBase> getQueue(int idx);

  boost::shared_ptr<QueueBase> getInputQueue();
  boost::shared_ptr<QueueBase> getOutputQueue();

 private:
  void schedule();

  void addSeat() {num_active_threads_.fetch_add(1);}
  void removeSeat() {num_active_threads_.fetch_sub(1);}

  boost::shared_ptr<boost::asio::io_service> ios_;
  boost::shared_ptr<boost::asio::io_service::work> ios_work_;

  boost::thread_group workers_;
  boost::thread_group scheduler_;

  int num_stages_;
  int num_threads_;
  mutable boost::atomic<int> num_active_threads_;

  std::vector<StageBase*> stages_;
  std::vector<Queue_ptr> queues_;

  std::map<std::string, boost::any> constants_;

  // Unfinished stages that require dynamic task dispatching
  std::deque<StageBase*> pending_stages_;
};

template <
  typename U, typename V, 
  int IN_DEPTH, int OUT_DEPTH
> 
bool Pipeline::addStage(int idx,
    Stage<U, V, IN_DEPTH, OUT_DEPTH> *stage)
{
  if (idx < 0 || idx >= num_stages_) {
    LOG(ERROR) << "index out of bound";
    return false;
  }
  if (stages_[idx]) {
    LOG(WARNING) << "Overwritting existing stage at idx=" << idx;
  }

  // bind input queue for stage if the stage requires it
  if (idx > 0 && IN_DEPTH > 0) {
    stage->input_queue_ = queues_[idx];
  }

  // create output queue for the stage if it requires it
  if (idx < num_stages_-1 && OUT_DEPTH <= 0) {
    LOG(ERROR) << "Intermediate stage must have an output queue";
    return false;
  }
  if (OUT_DEPTH > 0) {
    // if the output queue is already set, skip creation 
    if (!queues_[idx+1]) {
      boost::shared_ptr<QueueBase> output_queue(new Queue<V, OUT_DEPTH>);
      queues_[idx+1] = output_queue;
    }

    stage->output_queue_ = queues_[idx+1];
  }
  stages_[idx] = stage;
  if (idx > 0) {
    stages_[idx-1]->linkStage(stage);
  }
  stage->pipeline_ = this;

  return true;
}

template <typename T>
bool Pipeline::addConst(std::string key, T val) {
  if (constants_.count(key)) {
    LOG(ERROR) << key << " already exists in the constant table";
    return false;
  }
  constants_[key] = val;
  return true;
}

template <typename CompletionHandler>
void Pipeline::post(CompletionHandler handler) {
  ios_->post(handler);
}

} // namespace kestrelFlow
#endif

