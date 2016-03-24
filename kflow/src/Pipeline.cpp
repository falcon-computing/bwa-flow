#include "Pipeline.h"

namespace kestrelFlow 
{

Pipeline::Pipeline(int _num_stages): 
  num_stages(_num_stages),
  stages(_num_stages, NULL),
  queues(_num_stages+1, NULL_QUEUE_PTR) {}

void Pipeline::finalize() {
  stages[0]->final();
}

boost::any Pipeline::getConst(std::string key) {
  if (!constants.count(key)) {
    return 0;
  }
  else {
    return constants[key];
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

} // namepspace kestrelFlow
