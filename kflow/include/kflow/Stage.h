#ifndef KFLOW_STAGE_H
#define KFLOW_STAGE_H

#include <boost/atomic.hpp>
#include <boost/type_traits.hpp>
#include <typeinfo>

#include "Common.h"
#include "Queue.h"

// for testing purpose
#ifndef TEST_FRIENDS_LIST
#define TEST_FRIENDS_LIST
#endif

namespace kestrelFlow {

// base class for the template class Stage
class StageBase {
 friend class Pipeline;
 friend class OccupancyCounter;
 TEST_FRIENDS_LIST

 public:
  StageBase(int _num_workers=1, bool is_dyn=true);
  ~StageBase();

  virtual void start();
  virtual void stop();
  virtual void wait();

  bool isDynamic();

  boost::any getConst(std::string key);

 protected:
  void final();

  // inform the subsequent stage that there
  // is no more new records for the queue
  void finalize();

  bool isFinal();

  // Get the number of current active worker threads
  int getNumThreads();

  int getMaxNumThreads();

  // profiling counters for each worker:
  // <read_block_time, compute_time, write_block_time, total_time>
  uint64_t **perf_meters_;

  Pipeline*   pipeline_;

  // beginning and end timestamps of all stage workers
  uint64_t  start_ts_;
  uint64_t  end_ts_;

 private:
  // Execute one task from this stage, if the input_queue is empty
  // or the output_queue is full, no work will be posted and 
  // the call will return false
  virtual bool execute() = 0;

  // Function body of workers of static stage
  virtual void worker_func(int wid) = 0;

  void addSeat() {num_active_threads_.fetch_add(1);}
  void removeSeat() {num_active_threads_.fetch_sub(1);}

  bool        is_dynamic_;
  int         num_workers_;
  StageBase*  next_stage_;

  mutable boost::atomic<int>  num_active_threads_;
  mutable boost::atomic<int>  num_finalized_;
  mutable boost::atomic<bool> is_final_;
  boost::thread_group         worker_threads_;
};

template <
  typename U, 
  typename V, 
  int IN_DEPTH  = 64,
  int OUT_DEPTH = 64 
>
class Stage : public StageBase {
 friend class Pipeline;
 TEST_FRIENDS_LIST

 public:
  Stage(int num_workers=1, bool is_dyn=true):
      StageBase(num_workers, is_dyn) {}

 protected:
  Queue<U, IN_DEPTH>*  input_queue_;
  Queue<V, OUT_DEPTH>* output_queue_;
};

} // namespace kestrelFlow
#endif
