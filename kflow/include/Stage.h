#ifndef STAGE_H
#define STAGE_H

#include <boost/atomic.hpp>
#include <boost/type_traits.hpp>
#include <typeinfo>

#include "Common.h"
#include "Queue.h"

// for testing purpose
#ifndef TEST_FRIENDS_LIST
#define TEST_FRIENDS_LIST
#endif

namespace kestrelFlow 
{

// base class for the template class Stage
class StageBase 
{
  friend class Pipeline;
  TEST_FRIENDS_LIST

  public:
    StageBase(int _num_workers=1);

    void start();
    void stop();
    void wait();

  protected:
    virtual void worker_func() = 0;

    void final();

    // inform the subsequent stage that there
    // is no more new records for the queue
    void finalize();

    bool isFinal();

  private:
    int                 num_workers;
    boost::atomic<int>  num_finalized;
    boost::atomic<bool> is_final;

    StageBase* next_stage;

    boost::thread_group worker_threads;
};

template <
  typename U, 
  typename V, 
  int IN_DEPTH  = 64,
  int OUT_DEPTH = 64 
>
class Stage : public StageBase 
{
  friend class Pipeline;
  TEST_FRIENDS_LIST

  public:
    Stage(int _num_workers=1): StageBase(_num_workers) {}

  protected:
    void deleteIfPtr(U &obj, boost::true_type) {delete obj;}
    void deleteIfPtr(U &obj, boost::false_type) {}

    Queue<U, IN_DEPTH>*  input_queue;
    Queue<V, OUT_DEPTH>* output_queue;
};

StageBase::StageBase(int _num_workers): 
  num_workers(_num_workers),
  num_finalized(0),
  is_final(false), 
  next_stage(NULL) 
{
  if (_num_workers<1) {
    throw paramError("Invalid parameters");
  }
}

void StageBase::start() {
  if (worker_threads.size()) {
    worker_threads.interrupt_all();
    worker_threads.join_all();
  }
  for (int i=0; i<num_workers; i++) 
  {
    worker_threads.create_thread(
        boost::bind(&StageBase::worker_func, this));
  }
}

void StageBase::stop() {
  if (worker_threads.size()) {
    worker_threads.interrupt_all();
    worker_threads.join_all();
  }
}

void StageBase::wait() {
  if (worker_threads.size()) {
    worker_threads.join_all();
  }
}

void StageBase::final() {
  is_final.store(true);   
}

void StageBase::finalize() {
  num_finalized.fetch_add(1);
  if (next_stage && 
      num_finalized.load() == num_workers) 
  {
    next_stage->final();
    DLOG(INFO) << "Stage is finalized";
  }
}

bool StageBase::isFinal() {
  return is_final.load();
}


} // namespace kestrelFlow
#endif
