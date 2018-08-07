#ifndef KFLOW_SOURCESTAGE_H
#define KFLOW_SOURCESTAGE_H

#include "Stage.h"

// for testing purpose
#ifndef TEST_FRIENDS_LIST
#define TEST_FRIENDS_LIST
#endif

namespace kestrelFlow {

template <
  typename V, 
  int OUT_DEPTH = 64 
>
class SourceStage : public Stage<void, V, 0, OUT_DEPTH> {
 public:
  // force one worker for IO stages
  SourceStage();

 protected:
  virtual void compute() = 0;

  void pushOutput(V const & item);

 private:
  // Only need to overwrite vtable, do nothing
  bool execute() {;}
  int execute_new() {;}
  void worker_func(int wid);
};

template <typename V, int OUT_DEPTH>
SourceStage<V, OUT_DEPTH>::SourceStage(): 
  Stage<void, V, 0, OUT_DEPTH>(1, false)
{}

template <typename V, int OUT_DEPTH>
void SourceStage<V, OUT_DEPTH>::pushOutput(V const & item) {
  if (!this->getOutputQueue()) {
    return; 
  }
  this->getOutputQueue()->push(item);
}

template <typename V, int OUT_DEPTH>
void SourceStage<V, OUT_DEPTH>::worker_func(int wid) {
  if (!this->getOutputQueue()) {
    LOG(ERROR) << "Empty output queue is not allowed";
    return;
  }

  try {
    // call user-defined compute function
    compute(); 
  } 
  catch (boost::thread_interrupted &e) {
    DLOG_IF(INFO, FLAGS_v >= 2) << "Worker thread is interrupted";
    return;
  }
  // inform the next Stage there will be no more
  // output records
  this->finalize();

  DLOG(INFO) << "SourceStage is finished";
}

} // namespace kestrelFlow
#endif 
