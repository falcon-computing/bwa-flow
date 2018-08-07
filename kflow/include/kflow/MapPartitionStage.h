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
  MapPartitionStage(int n_workers=1, bool is_dyn=true):
      Stage<U, V, IN_DEPTH, OUT_DEPTH>(n_workers, is_dyn) {;}

  bool execute();
  int execute_new() { return (execute() ? 0 : 2); }

 protected:
  virtual void compute(int wid) = 0;

  bool getInput(U &item);
  void pushOutput(V const & item);

 private:
  // Function body of worker threads
  void worker_func(int wid);

  // Function body of execute function
  void execute_func();
};

template <
  typename U, typename V, 
  int IN_DEPTH, int OUT_DEPTH
>
bool MapPartitionStage<U, V, IN_DEPTH, OUT_DEPTH>::getInput(U &item) {
  if (!this->getInputQueue()) {
    return false; 
  }
  return this->getInputQueue()->async_pop(item);
}

template <
  typename U, typename V, 
  int IN_DEPTH, int OUT_DEPTH
>
void MapPartitionStage<U, V, IN_DEPTH, OUT_DEPTH>::pushOutput(V const & item) {
  if (!this->getOutputQueue()) {
    return; 
  }
  this->getOutputQueue()->push(item);
}

template <
  typename U, typename V, 
  int IN_DEPTH, int OUT_DEPTH
>
bool MapPartitionStage<U, V, IN_DEPTH, OUT_DEPTH>::execute() {

  // Return false if input queue is empty or max num_worker_threads reached
  if (this->getInputQueue()->empty() || 
      this->getNumThreads() >= this->getMaxNumThreads()) {
    return false;
  }

  // Post work from the compute function
  this->pipeline_->post(boost::bind(
        &MapPartitionStage<U, V, IN_DEPTH, OUT_DEPTH>::execute_func, this));

  return true;
}

template <
  typename U, typename V, 
  int IN_DEPTH, int OUT_DEPTH
>
void MapPartitionStage<U, V, IN_DEPTH, OUT_DEPTH>::execute_func() {

  try {
    OccupancyCounter seat(this->pipeline_, this);

    int n_workers = this->getNumThreads();

    DLOG(INFO) << "Post a work for MapPartitionStage, there are "
      << n_workers
      << " active threads in this stage";

    // call user-defined compute function
    compute(n_workers-1); 
  } 
  catch (boost::thread_interrupted &e) {
    DLOG_IF(INFO, FLAGS_v >= 2) << "Worker thread is interrupted";
    return;
  }
  catch (std::runtime_error &e) {
    LOG(ERROR) << "Error in compute(): " << e.what();  
  }
}

template <
  typename U, typename V, 
  int IN_DEPTH, int OUT_DEPTH
>
void MapPartitionStage<U, V, IN_DEPTH, OUT_DEPTH>::worker_func(int wid) {

  if (this->isDynamic()) {
    DLOG(WARNING) << "Dynamic stages are not supposed to "
      << "start worker threads, exiting.";
    return;
  }
    
  if (!this->getInputQueue() || !this->getOutputQueue()) {
    DLOG(ERROR) << "Empty input/output queue is not allowed";
    return;
  }

  try {
    // call user-defined compute function
    compute(wid); 
  } 
  catch (boost::thread_interrupted &e) {
    DLOG_IF(INFO, FLAGS_v >= 2) << "Worker thread is interrupted";
  }
  // inform the next Stage there will be no more
  // output records
  this->finalize();

  DLOG(INFO) << "MapPartition worker thread is terminated";
}

} // namespace kestrelFlow
#endif 
