#ifndef KFLOW_PIPELINE_H
#define KFLOW_PIPELINE_H

#include <boost/any.hpp>
#include <boost/asio.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
#include <deque>

#include "Common.h"
#include "Stage.h"
#include "MegaPipe.h"

// for testing purpose
#ifndef TEST_FRIENDS_LIST
#define TEST_FRIENDS_LIST
#endif

namespace kestrelFlow {

class Pipeline {
  friend class OccupancyCounter;
  friend class MegaPipe;
  TEST_FRIENDS_LIST

 public:
  Pipeline(int num_stages, int num_threads);
  ~Pipeline();

  template <typename U, typename V,
            int IN_DEPTH,  int OUT_DEPTH> 
  bool addStage(int idx, Stage<U, V, IN_DEPTH, OUT_DEPTH> *stage);

  template <typename U, typename V,
            int IN_DEPTH,  int OUT_DEPTH> 
  bool addAccxBckStage(int idx, Stage<U, V, IN_DEPTH, OUT_DEPTH> *stage, float init_priority = 1.0);

  // branch stage[idx] of pipeline
  void branch(Pipeline& pipeline, int idx); 

  // bind the output of pipeline to input of stage[idx]
  void converge(Pipeline& pipeline, int idx);

  // bind the output of stage[idx] to pipeline
  void diverge(Pipeline& pipeline, int idx);

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

  MegaPipe * megapipe_ = NULL;

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
  if (IN_DEPTH > 0) {
    if (idx > 0) {
      stage->input_queue_ = queues_[idx];
    }
    else {
      // create an input queue
      boost::shared_ptr<QueueBase> input_queue(
          new Queue<U, IN_DEPTH>);
      queues_[idx] = input_queue;
      stage->input_queue_ = input_queue;
    }
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

template <
  typename U, typename V,
  int IN_DEPTH, int OUT_DEPTH
>
bool Pipeline::addAccxBckStage(int idx,
    Stage<U, V, IN_DEPTH, OUT_DEPTH> *accx_bck_stage,
    float init_priority) {
  //if (stages_[idx] == nullptr) {
  //  return addStage(idx, accx_bck_stage);
  //}

  
  StageBase *stage = stages_[idx];  
  if (stage->accx_backend_stage_ != nullptr) {
    DLOG(WARNING) << "Overwriting the existing accelerator backend for stage " << idx;
  }
  stage->use_accx_ = false;

  if (IN_DEPTH > 0) {
    // create the internal queue for accx input
    int accx_load_queue_depth = (int)((init_priority+1)*accx_bck_stage->getMaxNumThreads());
    boost::shared_ptr<QueueBase> accx_load_queue(new Queue<U, IN_DEPTH>(accx_load_queue_depth));
    accx_bck_stage->input_queue_ = accx_load_queue;
  }
  accx_bck_stage->output_queue_ = stage->output_queue_;

  std::vector<StageBase *> shadow_downstream_stages(stage->downstream_stages_);
  for (int sid = 0; sid < shadow_downstream_stages.size(); sid++)
    accx_bck_stage->linkStage(shadow_downstream_stages[sid]);

  stage->linkStage(accx_bck_stage);
  stage->accx_backend_stage_ = accx_bck_stage;
  stage->accx_load_queue_ = accx_bck_stage->input_queue_;
  stage->accx_priority_ = init_priority;
  stage->use_accx_ = true;

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

