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

    boost::any getConst(std::string key);

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
    Pipeline*  pipeline;

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

} // namespace kestrelFlow
#endif
