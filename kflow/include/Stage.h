#ifndef STAGE_H
#define STAGE_H

#include "Common.h"
#include "Queue.h"

// for testing purpose
#ifndef TEST_FRIENDS_LIST
#define TEST_FRIENDS_LIST
#endif

namespace kestrelFlow 
{

class StageBase {
  public:
    virtual void start() = 0;
    virtual void stop() = 0;

    // produce a single output record
    virtual void compute() = 0;
};

template <
  typename U, int IN_DEPTH,
  typename V, int OUT_DEPTH
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
    // blocking call to get one element from input_queue
    U readInput();

    // blocking call to put one element to output_queue
    void writeOutput(V item);

  private:
    void worker_func();

    int num_workers;

    boost::thread_group worker_threads;

    Queue<U, IN_DEPTH>*  input_queue;
    Queue<V, OUT_DEPTH>* output_queue;
};

template <
  typename U, int IN_DEPTH,
  typename V, int OUT_DEPTH
>
Stage<U, IN_DEPTH, V, OUT_DEPTH>::Stage(int _num_workers): 
  num_workers(_num_workers)
{
  if (_num_workers<1) {
    throw paramError("Invalid parameters");
  }
}

template <
  typename U, int IN_DEPTH,
  typename V, int OUT_DEPTH
>
void Stage<U, IN_DEPTH, V, OUT_DEPTH>::start() {

  worker_threads.interrupt_all();
  worker_threads.join_all();

  for (int i=0; i<num_workers; i++) 
  {
    worker_threads.create_thread(
        boost::bind(&Stage<U, IN_DEPTH, V, OUT_DEPTH>::worker_func, this));
  }  
}

template <
  typename U, int IN_DEPTH,
  typename V, int OUT_DEPTH
>
void Stage<U, IN_DEPTH, V, OUT_DEPTH>::stop() {
  worker_threads.interrupt_all();
  worker_threads.join_all();
}

template <
  typename U, int IN_DEPTH,
  typename V, int OUT_DEPTH
>
void Stage<U, IN_DEPTH, V, OUT_DEPTH>::worker_func() {
  
  while (true) {
    try {
      compute(); 
    } 
    catch (boost::thread_interrupted &e) {
      DLOG(INFO) << "Worker thread is interrupted";
      break;
    }
  }
  DLOG(INFO) << "Worker terminated";
}

template <
  typename U, int IN_DEPTH,
  typename V, int OUT_DEPTH
>
U Stage<U, IN_DEPTH, V, OUT_DEPTH>::readInput()
{
  U ret;
  if (input_queue) {
    input_queue->pop(ret);
  }
  return ret;
}

template <
  typename U, int IN_DEPTH,
  typename V, int OUT_DEPTH
>
void Stage<U, IN_DEPTH, V, OUT_DEPTH>::writeOutput(V item)
{
  if (output_queue) {
    output_queue->push(item);
  }
}
}
#endif
