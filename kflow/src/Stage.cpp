#include "Pipeline.h"
#include "Stage.h"

namespace kestrelFlow 
{
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

boost::any StageBase::getConst(std::string key) {
  return pipeline->getConst(key);
}

} // namepsace kestrelFlow
