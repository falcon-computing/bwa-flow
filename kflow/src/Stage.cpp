#include <iomanip>
#include "Pipeline.h"
#include "Stage.h"

namespace kestrelFlow {

StageBase::StageBase(int num_workers, bool is_dyn): 
  is_dynamic_(is_dyn),
  num_workers_(num_workers),
  pipeline_(NULL),
  num_upstream_stages_(0),
  num_active_threads_(0),
  num_finalized_workers_(0),
  num_finalized_upstream_stages_(0),
  is_final_(false),
  use_accx_(false),
  accx_backend_stage_(NULL)
{
 // if (num_workers<1) {
 //   throw paramError("Invalid parameters");
 // }

  if (is_dynamic_) {
    DLOG(INFO) << "Created stage of maximum " << num_workers << " workers";
  }
  else {
    DLOG(INFO) << "Created stage of " << num_workers << " workers";
  }

  // Deprecated
  perf_meters_ = new uint64_t*[num_workers];
  for (int i = 0; i < num_workers; i++) {
    perf_meters_[i] = new uint64_t[4]();
  }
  worker_thread_ptr_ = new boost::thread*[num_workers];
}

StageBase::~StageBase() {
  for (int i = 0; i < num_workers_; i++) {
    delete [] perf_meters_[i];
  }
  delete [] perf_meters_;
  delete [] worker_thread_ptr_;
}

void StageBase::start() {
  if (is_dynamic_) return;

  if (worker_threads_.size()) {
    worker_threads_.interrupt_all();
    worker_threads_.join_all();
  }
  // set the start timer 
  start_ts_ = getUs();

  for (int i = 0; i < num_workers_; i++) {
    worker_thread_ptr_[i] = worker_threads_.create_thread(
        boost::bind(&StageBase::worker_func, this, i));
  }
}

void StageBase::stop() {
  if (is_dynamic_) return;
  if (worker_threads_.size()) {
    worker_threads_.interrupt_all();
    worker_threads_.join_all();
  }
}

void StageBase::wait() {
  if (is_dynamic_) return;
  if (worker_threads_.size()) {
    worker_threads_.join_all();
  }
}

void StageBase::final() {
  if (!is_final_.load()) {
    num_finalized_upstream_stages_.fetch_add(1) ;
    DLOG(INFO) << "Current finalized upstream stages: " << num_finalized_upstream_stages_;
    if (num_finalized_upstream_stages_.load() >= num_upstream_stages_) {
      is_final_.store(true);   
      DLOG(INFO) << "Stage is finalized";
    }
  }
}

void StageBase::finalize() {
  if (downstream_stages_.empty()) return;
  if (is_dynamic_) { 
    // if dynamic this function will be called by pipeline
    for (int i = 0; i < downstream_stages_.size(); i++) {
      downstream_stages_[i]->final();
    }
  }
  else {
    // otherwise this will be called by compute()
    if (incFinalizedThreads() == num_workers_) {
      for (int i = 0; i < downstream_stages_.size(); i++) {
        downstream_stages_[i]->final();
      }
    }
  }
}

bool StageBase::isDynamic() {
  return is_dynamic_;
}

bool StageBase::isFinal() {
  return is_final_.load();
}

boost::any StageBase::getConst(std::string key) {
  return pipeline_->getConst(key);
}

int StageBase::getNumThreads() {
  if (is_dynamic_) {
    return num_active_threads_.load();
  }
  else {
    return num_workers_;
  }
}

int StageBase::getNumActiveThreads() {
  return num_active_threads_.load();
}

int StageBase::incFinalizedThreads() {
  boost::lock_guard<boost::mutex> guard(mtx_);
  num_finalized_workers_.fetch_add(1); 
  int val = num_finalized_workers_.load();
  return val;
}

int StageBase::getMaxNumThreads() {
  return num_workers_;
}

void StageBase::linkStage(StageBase* next_stage) {
  // never link a stage to itself
  if (this == next_stage) return;

  downstream_stages_.push_back(next_stage);
  next_stage->num_upstream_stages_++;
  if (use_accx_ && accx_backend_stage_ != nullptr) {
    accx_backend_stage_->linkStage(next_stage);
  }
}

} // namepsace kestrelFlow
