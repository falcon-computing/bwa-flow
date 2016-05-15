#include <iomanip>
#include "Pipeline.h"
#include "Stage.h"

namespace kestrelFlow {

StageBase::StageBase(int num_workers, bool is_dyn): 
  is_dynamic_(is_dyn),
  num_workers_(num_workers),
  next_stage_(NULL),
  pipeline_(NULL),
  num_active_threads_(0),
  num_finalized_(0),
  is_final_(false)
{
  if (num_workers<1) {
    throw paramError("Invalid parameters");
  }

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
}

StageBase::~StageBase() {
  for (int i = 0; i < num_workers_; i++) {
    delete [] perf_meters_[i];
  }
  delete [] perf_meters_;
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
    worker_threads_.create_thread(
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
  is_final_.store(true);   
  DLOG(INFO) << "Stage is finalized";
}

void StageBase::finalize() {
  if (!next_stage_) return;
  if (is_dynamic_) {
    next_stage_->final();
  }
  else {
    num_finalized_.fetch_add(1);
    if (num_finalized_.load() == num_workers_) {
      next_stage_->final();
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

int StageBase::getMaxNumThreads() {
  return num_workers_;
}

/*
std::string StageBase::printPerf() {
  using namespace std;
  stringstream ss;
  ss << "Stage total time: " << (double)(end_ts - start_ts)/1e3 << "ms\n";
  ss << "===============================================================\n";
  ss << "| worker | load time | compute time | store time | total time |\n";
  for (int i=0; i<num_workers; i++) {
    ss <<  "| " << setw(4)  << setfill(' ') << i << "  ";
    ss << " | " << setw(7)  << setfill(' ') << perf_meters[i][0] << "us";
    ss << " | " << setw(10) << setfill(' ') << perf_meters[i][1] << "us";
    ss << " | " << setw(8)  << setfill(' ') << perf_meters[i][2] << "us";
    ss << " | " << setw(8)  << setfill(' ') << perf_meters[i][3] << "us";
    ss << " |\n";
  }
  ss << "===============================================================\n";
  return ss.str();
}
*/

} // namepsace kestrelFlow
