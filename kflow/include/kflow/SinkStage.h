#ifndef SINKSTAGE_H
#define SINKSTAGE_H

#include "Stage.h"

// for testing purpose
#ifndef TEST_FRIENDS_LIST
#define TEST_FRIENDS_LIST
#endif

namespace kestrelFlow 
{

template <
  typename U, 
  int IN_DEPTH = 64 
>
class SinkStage : 
  public Stage<U, void, IN_DEPTH, 0>
{
  public:
    // force one worker for IO stages
    SinkStage();

  protected:
    virtual void compute() = 0;

    bool getInput(U &item);

  private:
    void worker_func(int wid);
    Queue<U, IN_DEPTH>* src_queue;
};

template <typename U, int IN_DEPTH>
SinkStage<U, IN_DEPTH>::SinkStage(): 
  Stage<U, void, IN_DEPTH, 0>(1) 
{}

template <typename U, int IN_DEPTH>
bool SinkStage<U, IN_DEPTH>::getInput(U &item) 
{
  if (!src_queue) {
    return false; 
  }
  return src_queue->async_pop(item);
}

template <typename U, int IN_DEPTH>
void SinkStage<U, IN_DEPTH>::worker_func(int wid) {

#ifndef DISABLE_PROFILING
  uint64_t start_ts = getUs();
#endif

  if (!this->input_queue) {
    LOG(ERROR) << "Empty input queue is not allowed";
    return;
  }

  src_queue = dynamic_cast<
    Queue<U, IN_DEPTH>*>(this->input_queue);

  if (!src_queue && this->input_queue) {
    LOG(ERROR) << "Wrong input queue type";
    return;
  }

  try {
#ifndef DISABLE_PROFILING
    uint64_t start_ts = getUs();
#endif
    // call user-defined compute function
    compute(); 

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
