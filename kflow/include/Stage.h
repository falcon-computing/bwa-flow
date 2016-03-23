#ifndef STAGE_H
#define STAGE_H

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
class StageBase {
  public:
    virtual void start() = 0;
    virtual void stop() = 0;
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
    Stage(int _num_workers=1);

    void start();
    void stop();

  protected:
    // produce a single output record
    virtual V compute(U const & in) = 0;

  private:
    void worker_func();

    void deleteIfPtr(U &obj, boost::true_type);
    void deleteIfPtr(U &obj, boost::false_type);

    int num_workers;

    boost::thread_group worker_threads;

    Queue<U, IN_DEPTH>*  input_queue;
    Queue<V, OUT_DEPTH>* output_queue;
};

template <
  typename U, typename V, 
  int IN_DEPTH, int OUT_DEPTH
>
Stage<U, V, IN_DEPTH, OUT_DEPTH>::Stage(int _num_workers): 
  num_workers(_num_workers)
{
  if (_num_workers<1) {
    throw paramError("Invalid parameters");
  }
}

template <
  typename U, typename V, 
  int IN_DEPTH, int OUT_DEPTH
>
void Stage<U, V, IN_DEPTH, OUT_DEPTH>::start() {

  worker_threads.interrupt_all();
  worker_threads.join_all();

  for (int i=0; i<num_workers; i++) 
  {
    worker_threads.create_thread(
        boost::bind(&Stage<U, V, IN_DEPTH, OUT_DEPTH>::worker_func, this));
  }  
}

template <
  typename U, typename V, 
  int IN_DEPTH, int OUT_DEPTH
>
void Stage<U, V, IN_DEPTH, OUT_DEPTH>::stop() {
  worker_threads.interrupt_all();
  worker_threads.join_all();
}

template <
  typename U, typename V, 
  int IN_DEPTH, int OUT_DEPTH
>
void Stage<U, V, IN_DEPTH, OUT_DEPTH>::worker_func() {
  
  while (true) {
    try {
      // first read input from input queue
      U input;
      if (input_queue) {
        input_queue->pop(input);
      }

      // call user-defined compute function
      V output = compute(input); 

      // write result to output_queue
      if (output_queue) {
        output_queue->push(output);
      }

      // free input if it is a pointer
      deleteIfPtr(input, boost::is_pointer<U>());
    } 
    catch (boost::thread_interrupted &e) {
      DLOG(INFO) << "Worker thread is interrupted";
      break;
    }
  }
  DLOG(INFO) << "Worker terminated";
}

template <
  typename U, typename V, 
  int IN_DEPTH, int OUT_DEPTH
>
void Stage<U, V, IN_DEPTH, OUT_DEPTH>::deleteIfPtr(U &obj, boost::true_type) 
{
  delete obj;
}

template <
  typename U, typename V, 
  int IN_DEPTH, int OUT_DEPTH
>
void Stage<U, V, IN_DEPTH, OUT_DEPTH>::deleteIfPtr(U &obj, boost::false_type) 
{
  ;
}

} // namespace kestrelFlow
#endif
