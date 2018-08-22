#ifndef KFLOW_SINKSTAGE_H
#define KFLOW_SINKSTAGE_H

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
    SinkStage(int n_workers = 1, bool is_dyn = false);

  protected:
    virtual void compute(int wid) = 0;

    bool getInput(U &item);

  private:
    bool execute() {;}
    int execute_new() {;}
    void worker_func(int wid);
};

template <typename U, int IN_DEPTH>
SinkStage<U, IN_DEPTH>::SinkStage(int n_workers, bool is_dyn): 
  Stage<U, void, IN_DEPTH, 0>(n_workers, is_dyn) 
{}

template <typename U, int IN_DEPTH>
bool SinkStage<U, IN_DEPTH>::getInput(U &item) 
{
  if (!this->getInputQueue()) {
    return false; 
  }
  Queue<U, IN_DEPTH>* queue = this->getInputQueue();
  return queue->async_pop(item);
}

template <typename U, int IN_DEPTH>
void SinkStage<U, IN_DEPTH>::worker_func(int wid) {

  if (!this->getInputQueue()) {
    LOG(ERROR) << "Empty input queue is not allowed";
    return;
  }

  try {
    // call user-defined compute function
    compute(wid); 
  } 
  catch (boost::thread_interrupted &e) {
    DLOG_IF(INFO, FLAGS_v >= 2) << "Worker thread is interrupted";
    return;
  }
  // inform the next Stage there will be no more
  // output records
  this->finalize();

  DLOG(INFO) << "Worker thread is terminated";
}

} // namespace kestrelFlow
#endif 
