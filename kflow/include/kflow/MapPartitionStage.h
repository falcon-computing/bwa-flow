#ifndef KFLOW_MAPPARTITIONSTAGE_H
#define KFLOW_MAPPARTITIONSTAGE_H

#include "Stage.h"

// for testing purpose
#ifndef TEST_FRIENDS_LIST
#define TEST_FRIENDS_LIST
#endif

namespace kestrelFlow {

template <
  typename U, 
  typename V, 
  int IN_DEPTH = 64,
  int OUT_DEPTH = 64
>
class MapPartitionStage : 
  public Stage<U, V, IN_DEPTH, OUT_DEPTH> 
{
 public:
  // force one worker for IO stages
  MapPartitionStage(int n=1);

 protected:
  virtual void compute(int wid) = 0;

  bool getInput(U &item);
  void pushOutput(V const & item);

 private:
  void worker_func(int wid);

  Queue<U, IN_DEPTH>*  src_queue_;
  Queue<V, OUT_DEPTH>* dst_queue_;
};

template <
  typename U,
  typename V, 
  int IN_DEPTH,
  int OUT_DEPTH
>
MapPartitionStage<U, V, IN_DEPTH, OUT_DEPTH>::MapPartitionStage(int n): 
  Stage<U, V, IN_DEPTH, OUT_DEPTH>(n) 
{}

template <
  typename U,
  typename V, 
  int IN_DEPTH,
  int OUT_DEPTH
>
bool MapPartitionStage<U, V, IN_DEPTH, OUT_DEPTH>::getInput(U &item) {
  if (!src_queue_) {
    return false; 
  }
  return src_queue_->async_pop(item);
}

template <
  typename U,
  typename V, 
  int IN_DEPTH,
  int OUT_DEPTH
>
void MapPartitionStage<U, V, IN_DEPTH, OUT_DEPTH>::pushOutput(V const & item) {
  if (!dst_queue_) {
    return; 
  }
  dst_queue_->push(item);
}

template <
  typename U,
  typename V, 
  int IN_DEPTH,
  int OUT_DEPTH
>
void MapPartitionStage<U, V, IN_DEPTH, OUT_DEPTH>::worker_func(int wid) {

#ifndef DISABLE_PROFILING
  uint64_t start_ts = getUs();
#endif

  if (!this->input_queue || !this->output_queue) {
    LOG(ERROR) << "Empty input/output queue is not allowed";
    return;
  }

  src_queue_ = dynamic_cast<Queue<U, IN_DEPTH>*>(this->input_queue);
  dst_queue_ = dynamic_cast<Queue<V, OUT_DEPTH>*>(this->output_queue);

  if (!src_queue_ || !dst_queue_) {
    LOG(ERROR) << "Wrong input/output queue type(s)";
    return;
  }

  try {
#ifndef DISABLE_PROFILING
    uint64_t start_ts = getUs();
#endif
    // call user-defined compute function
    compute(wid); 

#ifndef DISABLE_PROFILING
    // record compute time
    this->perf_meters[wid][1] += getUs() - start_ts;
#endif
  } 
  catch (boost::thread_interrupted &e) {
    VLOG(2) << "Worker thread is interrupted";
    return;
  }
  // inform the next Stage there will be no more
  // output records
  this->finalize();

#ifndef DISABLE_PROFILING
  // record total time
  this->perf_meters[wid][3] = getUs()-start_ts;
#endif
  this->end_ts = getUs();

  DLOG(INFO) << "Worker thread is terminated";
}

} // namespace kestrelFlow
#endif 
