#ifndef IOSTAGE_H
#define IOSTAGE_H

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
class FSourceStage : public Stage<void, V, 0, OUT_DEPTH>
{
  public:
    // force one worker for IO stages
    FSourceStage();

  protected:
    virtual void compute() = 0;

    void pushOutput(V const & item);

  private:
    void worker_func();

    Queue<V, OUT_DEPTH>* dst_queue;
};

template <
  typename U, 
  int IN_DEPTH = 64 
>
class FSinkStage : 
  public Stage<U, void, IN_DEPTH, 0>
{
  public:
    // force one worker for IO stages
    FSinkStage();

  protected:
    virtual void compute() = 0;

    bool getInput(U &item);

  private:
    void worker_func();    
    Queue<U, IN_DEPTH>* src_queue;
};

template <typename V, int OUT_DEPTH>
FSourceStage<V, OUT_DEPTH>::FSourceStage(): 
  Stage<void, V, 0, OUT_DEPTH>(1) 
{}

template <typename V, int OUT_DEPTH>
void FSourceStage<V, OUT_DEPTH>::pushOutput(V const & item) 
{
  if (!dst_queue) {
    return; 
  }
  dst_queue->push(item);
}

template <typename V, int OUT_DEPTH>
void FSourceStage<V, OUT_DEPTH>::worker_func() {

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

template <typename U, int IN_DEPTH>
FSinkStage<U, IN_DEPTH>::FSinkStage(): 
  Stage<U, void, IN_DEPTH, 0>(1) 
{}

template <typename U, int IN_DEPTH>
bool FSinkStage<U, IN_DEPTH>::getInput(U &item) {
  if (!src_queue) {
    return false; 
  }
  return src_queue->async_pop(item);
}

template <typename U, int IN_DEPTH>
void FSinkStage<U, IN_DEPTH>::worker_func() {

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
