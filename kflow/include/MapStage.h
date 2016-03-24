#ifndef MAPSTAGE_H
#define MAPSTAGE_H

#include "Stage.h"

// for testing purpose
#ifndef TEST_FRIENDS_LIST
#define TEST_FRIENDS_LIST
#endif

namespace kestrelFlow 
{

template <
  typename U, 
  typename V, 
  int IN_DEPTH  = 64,
  int OUT_DEPTH = 64 
>
class MapStage : public Stage<U, V, IN_DEPTH, OUT_DEPTH>
{
  public:
    MapStage(int _num_workers=1):
      Stage<U, V, IN_DEPTH, OUT_DEPTH>(_num_workers) {}

  protected:
    virtual V compute(U const & input) = 0;

  private:
    void worker_func() {
      if (!this->input_queue) {
        LOG(ERROR) << "Empty input queue for MapStage is not allowed";
        return;
      }

      Queue<U, IN_DEPTH>* input_queue = dynamic_cast<
        Queue<U, IN_DEPTH>*>(this->input_queue);

      Queue<V, OUT_DEPTH>* output_queue = dynamic_cast<
        Queue<V, OUT_DEPTH>*>(this->output_queue);

      if (!input_queue) {
        LOG(ERROR) << "Mismatched input queue type";
        return;
      }
      if (!output_queue && this->output_queue) {
        LOG(ERROR) << "Mismatched ouput queue type";
        return;
      }

      while (true) 
      {
        try {
          // first read input from input queue
          U input;
          bool ready = input_queue->async_pop(input);

          while (!this->isFinal() && !ready) {
            boost::this_thread::sleep_for(boost::chrono::microseconds(100));
            ready = input_queue->async_pop(input);
          }
          if (!ready) { 
            // this means isFinal() is true and input queue is empty
            break; 
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
          VLOG(2) << "Worker thread is interrupted";
          break;
        }
      }
      // inform the next Stage there will be no more
      // output records
      this->finalize();

      DLOG(INFO) << "Worker thread is terminated";
    }
};

} // namespace kestrelFlow
#endif 
