#include <iomanip>
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
  perf_meters = new uint64_t*[_num_workers];
  for (int i=0; i<_num_workers; i++) {
    perf_meters[i] = new uint64_t[4]();
  }
}

StageBase::~StageBase() {
  for (int i=0; i<num_workers; i++) {
    delete [] perf_meters[i];
  }
  delete [] perf_meters;
}

void StageBase::start() {
  if (worker_threads.size()) {
    worker_threads.interrupt_all();
    worker_threads.join_all();
  }
  // set the start timer 
  start_ts = getUs();

  for (int i=0; i<num_workers; i++) 
  {
    worker_threads.create_thread(
        boost::bind(&StageBase::worker_func, this, i));
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

} // namepsace kestrelFlow
