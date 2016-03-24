#ifndef PIPELINE_H
#define PIPELINE_H

#include "Common.h"
#include "Stage.h"

// for testing purpose
#ifndef TEST_FRIENDS_LIST
#define TEST_FRIENDS_LIST
#endif

namespace kestrelFlow 
{
class Pipeline 
{
  TEST_FRIENDS_LIST

  public:
    Pipeline(int _num_stages);

    template <typename U, typename V,
              int IN_DEPTH,  int OUT_DEPTH> 
    void addStage(int idx, Stage<U, V, IN_DEPTH, OUT_DEPTH> *stage);
    
    void start();
    void stop();
    void wait();
    void finalize();

    QueueBase* getInputQueue();
    QueueBase* getOutputQueue();

  private:
    int num_stages;
    std::vector<StageBase*> stages;
    std::vector<boost::shared_ptr<QueueBase> > queues;
};

Pipeline::Pipeline(int _num_stages): 
  num_stages(_num_stages),
  stages(_num_stages, NULL),
  queues(_num_stages+1, NULL_QUEUE_PTR) {}

void Pipeline::finalize() {
  stages[0]->final();
}

template <
  typename U, typename V, 
  int IN_DEPTH, int OUT_DEPTH
> 
void Pipeline::addStage(int idx,
    Stage<U, V, IN_DEPTH, OUT_DEPTH> *stage)
{
  if (idx < 0 || idx >= num_stages) {
    throw paramError("idx out of bound");
  }
  if (stages[idx]) {
    LOG(WARNING) << "Overwritting existing stage at idx=" << idx;
  }

  // bind input queue for stage if it exists
  if (idx == 0 && IN_DEPTH > 0) {
    boost::shared_ptr<QueueBase> queue(new Queue<U, IN_DEPTH>);
    queues[idx] = queue;
  }
  Queue<U, IN_DEPTH>* input_queue = 
    dynamic_cast<Queue<U, IN_DEPTH>*>(queues[idx].get());

  stage->input_queue = input_queue;

  // create output queue for the stage
  if (idx < num_stages-1 && OUT_DEPTH <= 0) {
    LOG(ERROR) << "Intermediate stage must have an output queue";
    return;
  }
  if (OUT_DEPTH > 0) {
    boost::shared_ptr<QueueBase> output_queue(new Queue<V, OUT_DEPTH>);
    queues[idx+1] = output_queue;

    stage->output_queue = dynamic_cast<
      Queue<V, OUT_DEPTH>*>(output_queue.get());
  }
  stages[idx] = stage;
  if (idx > 0) {
    stages[idx-1]->next_stage = stage;
  }
}

void Pipeline::start() 
{
  bool initialized = true;
  for (int i=0; i<num_stages; i++) {
    if (!stages[i]) {
      LOG(ERROR) << "Stage " << i << " is uninitialized";
      initialized = false;
    }
  }
  if (!initialized) {
    LOG(ERROR) << "Cannot start the pipeline due to previous errors";
  }
  for (int i=0; i<num_stages; i++) {
    stages[i]->start();
    VLOG(1) << "Start workers for stage " << i;
  }
}

void Pipeline::stop() 
{
  for (int i=0; i<num_stages; i++) {
    if (!stages[i]) {
      LOG(WARNING) << "Stage " << i << " is uninitialized";
    }
    else {
      stages[i]->stop();
    }
  }
}

void Pipeline::wait() 
{
  if (stages[num_stages-1]) {
    stages[num_stages-1]->wait();
  }
}

QueueBase* Pipeline::getInputQueue() {
  return queues[0].get();;
}

QueueBase* Pipeline::getOutputQueue() {
  return queues[num_stages].get();
}

} // namespace kestrelFlow
#endif

