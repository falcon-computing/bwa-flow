#ifndef KFLOW_MAPSTAGE_H
#define KFLOW_MAPSTAGE_H

#include "Stage.h"
#include "OccupancyCounter.h"

// for testing purpose
#ifndef TEST_FRIENDS_LIST
#define TEST_FRIENDS_LIST
#endif

namespace kestrelFlow {

template <
  typename U, 
  typename V, 
  int IN_DEPTH  = 64,
  int OUT_DEPTH = 64 
>
class MapStage : public Stage<U, V, IN_DEPTH, OUT_DEPTH> {
 public:
  MapStage(int n_workers=1, bool is_dyn=true):
    Stage<U, V, IN_DEPTH, OUT_DEPTH>(n_workers, is_dyn) {;}

  bool execute();

 protected:
  virtual V compute(U const & input) = 0;

 private:
  void deleteIfPtr(U &obj, boost::true_type) {delete obj;}
  void deleteIfPtr(U &obj, boost::false_type) {}

  // Function body of worker threads
  void worker_func(int wid);

  // Function body of execute function
  void execute_func(U input);
};

template <
  typename U, typename V, 
  int IN_DEPTH, int OUT_DEPTH
>
bool MapStage<U, V, IN_DEPTH, OUT_DEPTH>::execute() {

  // Return false if input queue is empty or max num_worker_threads reached
  if (this->getNumThreads() >= this->getMaxNumThreads() || 
      this->getOutputQueue()->almost_full()) {
    return false;
  }

  // Try to get one input from the input queue
  U input;
  bool ready = this->getInputQueue()->async_pop(input);

  if (!ready) {
    // return false if input queue is empty
    return false;
  }

  // Post work from the compute function
  this->pipeline_->post(boost::bind(
        &MapStage<U, V, IN_DEPTH, OUT_DEPTH>::execute_func, this, input));

  return true;
}

template <
  typename U, typename V, 
  int IN_DEPTH, int OUT_DEPTH
>
void MapStage<U, V, IN_DEPTH, OUT_DEPTH>::execute_func(U input) {

  try {
    OccupancyCounter seat(this->pipeline_, this);

    DLOG(INFO) << "Post a work for MapStage, there are "
      << this->getNumThreads()
      << " active threads in this stage";

    // call user-defined compute function
    V output = compute(input); 

    // write result to output_queue
    if (this->getOutputQueue()) {
      uint64_t start_ts = getUs();
      if (sizeof(V) != sizeof(int))
        this->getOutputQueue()->push(output);
      uint64_t end_ts = getUs();
      if (end_ts - start_ts >= 2000) {
        DLOG(WARNING) << "Output queue is full for " << end_ts - start_ts
                     << " us, blocking progress";
      }
    }

    // free input if it is a pointer
    deleteIfPtr(input, boost::is_pointer<U>());
  } 
  catch (boost::thread_interrupted &e) {
    DLOG_IF(INFO, FLAGS_v >= 2) << "Worker thread is interrupted";
  }
  catch (std::runtime_error &e) {
    LOG(ERROR) << "Error in compute(): " << e.what();  
  }
}

template <
  typename U, typename V, 
  int IN_DEPTH, int OUT_DEPTH
>
void MapStage<U, V, IN_DEPTH, OUT_DEPTH>::worker_func(int wid) {

  if (this->isDynamic()) {
    LOG(WARNING) << "Dynamic stages are not supposed to "
      << "start worker threads, exiting.";
    return;
  }
    
  if (!this->getInputQueue() || !this->getOutputQueue()) {
    throw std::runtime_error("Empty input/output queue is not allowed");
  }

  while (true) {
    try {
      // first read input from input queue
      U input;
      bool ready = this->getInputQueue()->async_pop(input);

      while (!this->isFinal() && !ready) {
        boost::this_thread::sleep_for(boost::chrono::microseconds(100));
        ready = this->getInputQueue()->async_pop(input);
      }
      if (!ready) { 
        // this means isFinal() is true and input queue is empty
        break; 
      }
      
      // call user-defined compute function
      V output = compute(input); 

      // write result to output_queue
      this->getOutputQueue()->push(output);

      // free input if it is a pointer
      deleteIfPtr(input, boost::is_pointer<U>());
    } 
    catch (boost::thread_interrupted &e) {
      DLOG_IF(INFO, FLAGS_v >= 2) << "Worker thread is interrupted";
      break;
    }
    catch (std::runtime_error &e) {
      LOG(INFO) << "compute() failed from " << e.what();
      return;
    }
  }
  // inform the next Stage there will be no more
  // output records
  this->finalize();

  DLOG(INFO) << "Map Worker thread is terminated";
}
} // namespace kestrelFlow
#endif 
