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
    bool addStage(int idx, Stage<U, V, IN_DEPTH, OUT_DEPTH> *stage);
    
    bool addConst(std::string key, RecordBase* val);

    RecordBase* getConst(std::string key);

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
    std::map<std::string, RecordBase*> constants;
};

template <
  typename U, typename V, 
  int IN_DEPTH, int OUT_DEPTH
> 
bool Pipeline::addStage(int idx,
    Stage<U, V, IN_DEPTH, OUT_DEPTH> *stage)
{
  if (idx < 0 || idx >= num_stages) {
    LOG(ERROR) << "index out of bound";
    return false;
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
    return false;
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
  stage->pipeline = this;

  return true;
}
} // namespace kestrelFlow
#endif

