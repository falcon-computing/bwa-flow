#ifndef KFLOW_MAPSTAGE_H
#define KFLOW_MAPSTAGE_H

#include "Stage.h"
#include "Pipeline.h"
#include "MegaPipe.h"
#include "OccupancyCounter.h"

#include <deque>

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
  int execute_new();

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
int MapStage<U, V, IN_DEPTH, OUT_DEPTH>::execute_new() {
  // Return STATUS-2 if output queue is half full
  if (this->getOutputQueue() && this->getOutputQueue()->almost_full()) {
    return 2;
  }

  // Before checking the input queue, merge accx input queue if accx is down (due to time/error)
  if ( !((StageBase*)this)->useAccx() &&
       ((StageBase*)this)->accx_backend_stage_ != NULL ) {
    while (this->getAccxQueue()->get_size() > 0) {
      U temp;
      bool isValid = this->getAccxQueue()->async_pop(temp);
      if (isValid) this->getInputQueue()->push(temp);
    }
  }
  // Try to get one input from the input queue
  if (this->inputQueueEmpty()) return 1;
  U input;
  bool ready = this->getInputQueue()->async_pop(input);
  if (!ready) {
    // return STATUS-1 if input queue is empty
    return 1;
  }

  // Try to get token
  if (((StageBase*)this)->useAccx()) {
    int accx_queue_size = this->getAccxQueue()->get_size();
    int accx_queue_capacity = this->getAccxQueue()->get_capacity();
    if ( accx_queue_size < accx_queue_capacity &&
         accx_queue_size <= ((StageBase*)this)->accx_backend_stage_->getNumActiveThreads() * ((StageBase*)this)->accx_priority_ ) {
      this->getAccxQueue()->push(input);
      return 0;
    }
  }

  this->execute_func(input);

  return 0;
}


template <
  typename U, typename V, 
  int IN_DEPTH, int OUT_DEPTH
>
void MapStage<U, V, IN_DEPTH, OUT_DEPTH>::execute_func(U input) {

  try {
    /* disabled pipeline-level threading control */
    OccupancyCounter seat(this->pipeline_, this); 

    //DLOG(INFO) << "Instantiate a work for MapStage, there are "
    //           << this->getNumThreads()
    //           << " active threads in this stage";

    // wait to get a cpu token
    while (!this->pipeline_->megapipe_->acqThrd()) {
      boost::this_thread::sleep_for(boost::chrono::microseconds(100));
    }
    // call user-defined compute function
    V output = compute(input); 
    // release a cpu token
    this->pipeline_->megapipe_->relThrd();

    // write result to output_queue
    if (this->getOutputQueue()) {
      uint64_t start_ts = getUs();
      if (sizeof(V) != sizeof(int))
        this->getOutputQueue()->push(output);
      uint64_t end_ts = getUs();
      if (end_ts - start_ts >= 500) {
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
    
  if (!this->getInputQueue() || (OUT_DEPTH > 0 && !this->getOutputQueue()) ) {
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
      if (OUT_DEPTH > 0)
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
