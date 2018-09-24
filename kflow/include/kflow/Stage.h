#ifndef KFLOW_STAGE_H
#define KFLOW_STAGE_H

#include <boost/atomic.hpp>
#include <boost/type_traits.hpp>
#include <typeinfo>

#include "Common.h"
#include "Queue.h"

// for testing purpose
#ifndef TEST_FRIENDS_LIST
#define TEST_FRIENDS_LIST
#endif

namespace kestrelFlow {

// base class for the template class Stage
class StageBase {
 friend class Pipeline;
 friend class MegaPipe;
 friend class OccupancyCounter;
 TEST_FRIENDS_LIST

 public:
  StageBase(int _num_workers=1, bool is_dyn=true);
  ~StageBase();

  StageBase* accx_backend_stage_;
  float accx_priority_ = 1.0;

  int getMaxNumThreads();
  int getNumActiveThreads();
  int incFinalizedThreads();

  virtual void start();
  virtual void stop();
  virtual void wait();

  bool isDynamic();
  int  getPriority() { return priority_; }

  bool useAccx() { return use_accx_; }
  void setUseAccx(bool flag) { use_accx_ = flag; }
  StageBase *getAccxStage() { return accx_backend_stage_; }

  boost::any getConst(std::string key);

  boost::thread *getThread(int i) { return worker_thread_ptr_[i]; }

 protected:
  void final();

  // inform the subsequent stage that there
  // is no more new records for the queue
  void finalize();

  bool isFinal();
  virtual bool inputQueueEmpty() = 0;

  // Get the number of current active worker threads
  int getNumThreads();

  // profiling counters for each worker:
  // <read_block_time, compute_time, write_block_time, total_time>
  uint64_t **perf_meters_;

  Pipeline*   pipeline_;

  // beginning and end timestamps of all stage workers
  uint64_t  start_ts_;
  uint64_t  end_ts_;

  Queue_ptr input_queue_;
  Queue_ptr output_queue_;
  Queue_ptr accx_load_queue_;

  std::vector<StageBase*> downstream_stages_;

 private:
  // Execute one task from this stage, if the input_queue is empty
  // or the output_queue is full, no work will be posted and 
  // the call will return false
  virtual bool execute() = 0;
  virtual int execute_new() = 0;

  // Function body of workers of static stage
  virtual void worker_func(int wid) = 0;

  // add a downstreams stage
  void linkStage(StageBase* next_stage);

  void addSeat() {num_active_threads_.fetch_add(1);}
  void removeSeat() {num_active_threads_.fetch_sub(1);}
  void removeSeatAndCheck() {
    num_active_threads_.fetch_sub(1);
    if (this->getNumThreads() == 0) {

    }
  }

  bool is_dynamic_;
  mutable boost::atomic<bool> use_accx_;
  int  num_workers_;
  int  num_upstream_stages_; 
  int  priority_ = 1;

  mutable boost::atomic<int>  num_active_threads_;
  mutable boost::atomic<int>  num_finalized_workers_;
  mutable boost::atomic<int>  num_finalized_upstream_stages_;
  mutable boost::atomic<bool> is_final_;
  boost::thread_group         worker_threads_;
  boost::thread             **worker_thread_ptr_;
  
  boost::mutex mtx_;
};

template <
  typename U, 
  typename V, 
  int IN_DEPTH  = 64,
  int OUT_DEPTH = 64 
>
class Stage : public StageBase {
 friend class Pipeline;
 TEST_FRIENDS_LIST

 public:
  Stage(int num_workers=1, bool is_dyn=true):
      StageBase(num_workers, is_dyn) {}

 public:
  Queue<U, IN_DEPTH>* getInputQueue() {
    if (IN_DEPTH == 0) return NULL;
    if (!input_queue_) {
      throw internalError("input queue is uninitialized");
    }
    Queue<U, IN_DEPTH>* queue = 
      dynamic_cast<Queue<U, IN_DEPTH>*>(input_queue_.get());

    if (!queue) {
      throw internalError("failed to dynamic_cast input queue");
    }

    return queue;
  }

  Queue<V, OUT_DEPTH>* getOutputQueue() {
    if (OUT_DEPTH == 0) return NULL;
    if (!output_queue_) {
      throw internalError("output queue is uninitialized");
    }
    Queue<V, OUT_DEPTH>* queue = 
      dynamic_cast<Queue<V, OUT_DEPTH>*>(output_queue_.get());

    if (!queue) {
      throw internalError("failed to dynamic_cast output queue");
    }

    return queue;
  }

  Queue<U, IN_DEPTH>* getAccxQueue() {
    if (IN_DEPTH == 0) return NULL;
    if (!accx_load_queue_) {
      throw internalError("accelerator input queue is uninitialized");
    }
    Queue<U, IN_DEPTH>* queue =
      dynamic_cast<Queue<U, IN_DEPTH>*>(accx_load_queue_.get());

    if (!queue) {
       throw internalError("failed to dynamic_cast accelerator input queue");
    }

    return queue;
  }

  bool inputQueueEmpty() {
    if (this->getInputQueue()) {
      if (!this->getInputQueue()->empty())
        return false;
      else {
        if (accx_backend_stage_ != nullptr)
          return this->getAccxQueue()->empty();
        else
          return true;
      }
    }
    else {
      return true;
    }
  }
};

} // namespace kestrelFlow
#endif
