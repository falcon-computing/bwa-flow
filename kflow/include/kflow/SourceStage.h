#ifndef SOURCESTAGE_H
#define SOURCESTAGE_H

#include "Stage.h"

// for testing purpose
#ifndef TEST_FRIENDS_LIST
#define TEST_FRIENDS_LIST
#endif

namespace kestrelFlow 
{

template <
  typename V, 
  int OUT_DEPTH = 64 
>
class SourceStage : public Stage<void, V, 0, OUT_DEPTH>
{
  public:
    // force one worker for IO stages
    SourceStage();

  protected:
    virtual void compute() = 0;

    void pushOutput(V const & item);

  private:
    void worker_func(int wid);

    Queue<V, OUT_DEPTH>* dst_queue;
};

template <typename V, int OUT_DEPTH>
SourceStage<V, OUT_DEPTH>::SourceStage(): 
  Stage<void, V, 0, OUT_DEPTH>(1) 
{}

template <typename V, int OUT_DEPTH>
void SourceStage<V, OUT_DEPTH>::pushOutput(V const & item) 
{
  if (!dst_queue) {
    return; 
  }
  dst_queue->push(item);
}

template <typename V, int OUT_DEPTH>
void SourceStage<V, OUT_DEPTH>::worker_func(int wid) 
{
#ifndef DISABLE_PROFILING
  uint64_t start_ts = getUs();
#endif

  if (!this->output_queue) {
    LOG(ERROR) << "Empty output queue is not allowed";
    return;
  }

  dst_queue = dynamic_cast<
    Queue<V, OUT_DEPTH>*>(this->output_queue);

  if (!dst_queue && this->output_queue) {
    LOG(ERROR) << "Wrong ouput queue type";
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
  // fix compute time
  this->perf_meters[wid][1] -= this->perf_meters[wid][2];

  // record total time
  this->perf_meters[wid][3] = getUs()-start_ts;
#endif
  this->end_ts = getUs();

  DLOG(INFO) << "Worker thread is terminated";
}

} // namespace kestrelFlow
#endif 
