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
    SinkStage();

  protected:
    virtual void compute() = 0;

    bool getInput(U &item);

  private:
    bool execute() {;}
    void worker_func(int wid);
};

template <typename U, int IN_DEPTH>
SinkStage<U, IN_DEPTH>::SinkStage(): 
  Stage<U, void, IN_DEPTH, 0>(1, false) 
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
    compute(); 
  } 
  catch (boost::thread_interrupted &e) {
    VLOG(2) << "Worker thread is interrupted";
    return;
  }
  // inform the next Stage there will be no more
  // output records
  this->finalize();

  DLOG(INFO) << "Worker thread is terminated";
}

} // namespace kestrelFlow
#endif 
