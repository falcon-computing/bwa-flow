#include "Pipeline.h"
#include "Stage.h"
#include "MegaPipe.h"

namespace kestrelFlow {

MegaPipe::MegaPipe( int num_threads, int num_accelerators ):
  num_threads_(num_threads),
  num_accelerators_(num_accelerators),
  num_avail_thrds_(num_threads),
  num_avail_accs_(num_accelerators)
{
}

MegaPipe::~MegaPipe()
{
  if (!dworkers_.size()) {
    dworkers_.interrupt_all();
    dworkers_.join_all();
  }
}

bool MegaPipe::addPipeline( Pipeline *pipeline, int priority )
{
  for (int sid = 0; sid < pipeline->stages_.size(); sid++) {
    StageBase *stage = pipeline->stages_[sid];
    stage->priority_ = priority;
  }
  
  bool is_added = false;
  for (int i = 0; i < pipelines_.size(); i++) {
    if (priority > priorities_[i]) {
      pipelines_.insert(pipelines_.begin()+i, pipeline);
      priorities_.insert(priorities_.begin()+i, priority);
      is_added = true;
      break;
    }
  }
  if (!is_added) {
    pipelines_.push_back(pipeline);
    priorities_.push_back(priority);
    is_added = true;
  }
  
  pipeline->megapipe_ = this;

  return is_added;
}

/* Start all pipelines.
 *   Start the delicated worker group for each static stage.
 *   Create a common worker group for all dynamic stages.
 */
void MegaPipe::start()
{
  // static part
  for (int pid = 0; pid < pipelines_.size(); pid++) {
    Pipeline *pipeline = pipelines_[pid];
    for (int sid = 0; sid < pipeline->stages_.size(); sid++) {
      StageBase *stage = pipeline->stages_[sid];
      if (!stage->isDynamic()) stage->start();
      if (stage->useAccx())
        stage->accx_backend_stage_->start();
    }
  }

  // dynamic part
  if (dworkers_.size()) {
    dworkers_.interrupt_all();
    dworkers_.join_all();
  }
  boost::this_thread::sleep_for(boost::chrono::milliseconds(1000));

  for (int pid = 0; pid < pipelines_.size(); pid++)
    for (int t = 0; t < pipelines_[pid]->num_threads_; t++)
      dworkers_.create_thread(boost::bind(&MegaPipe::dworker_func, this, pid));
  
  boost::this_thread::sleep_for(boost::chrono::milliseconds(1000));
}

/* Stop all pipelines.
 *   Stop the delicated worker group for each static stage.
 *   Stop the common worker group for all dynamic stages.
 */
void MegaPipe::stop()
{
  // dynamic part
  if (dworkers_.size()) {
    dworkers_.interrupt_all();
    dworkers_.join_all();
  }
}

/* Wait all pipelines to finish.
 *   Wait all static stages.
 *   Wait the common worker group for all dynamic stages.
 */
void MegaPipe::wait()
{
  for (int pid = 0; pid < pipelines_.size(); pid++) {
    Pipeline *pipeline = pipelines_[pid];
    // Only need to wait for the last stage
    if (!pipeline->stages_[pipeline->num_stages_-1]->isDynamic()) {
      pipeline->stages_[pipeline->num_stages_-1]->wait();
    }
  }
  // dynamic part
  if (dworkers_.size()) {
    dworkers_.join_all();
  }
}

/* Finalize all pipelines */
void MegaPipe::finalize()
{
  for (int pid = 0; pid < pipelines_.size(); pid++) {
    DLOG(INFO) << "Finalize pipeline " << pid;
    pipelines_[pid]->finalize();
  }
}

/* Dynamic worker procedure. */
void MegaPipe::dworker_func(int pipeline_id)
{
  std::deque<StageBase *> pipeline_structure;
  Pipeline *pipeline = pipelines_[pipeline_id];

  for (int sid = 0; sid < pipeline->stages_.size(); sid++) {
    StageBase *stage = pipeline->stages_[sid];
    if (!stage->isDynamic()) continue;
    pipeline_structure.push_back(stage);
  }

  while (pipeline_structure.size() != 0) {
    int rval = -1;
    for (int sid = pipeline_structure.size()-1; sid >= 0; sid--) {
      StageBase *stage = pipeline_structure[sid];
      rval = stage->execute_new();
      if (rval == 0) {
        break;
      }
      else if (rval == 1) {
        boost::this_thread::sleep_for(boost::chrono::milliseconds(5));
        if (sid == 0) {
          if (stage->isFinal()) {
            if (stage->getNumThreads()==0 && stage->inputQueueEmpty()) {
              if (stage->incFinalizedThreads() == pipeline->num_threads_)
                stage->finalize();
              pipeline_structure.pop_front();
            }
          }
          else {
            sid++;
          }
        }
      }
      else if (rval >= 2) {
        boost::this_thread::sleep_for(boost::chrono::milliseconds(5));
      }
    }
  }
}

bool MegaPipe::acqAccx()
{
  if ( num_avail_accs_.fetch_sub(1) <= 0 ) {
    num_avail_accs_++;
    return false;
  }
  return true;
}

void MegaPipe::relAccx()
{
  num_avail_accs_++;
}

bool MegaPipe::acqThrd()
{
  if ( num_avail_thrds_.fetch_sub(1) <= 0 ) {
    num_avail_thrds_++;
    return false;
  }
  return true;
}

void MegaPipe::relThrd()
{
  num_avail_thrds_++;
}

void MegaPipe::getThrd()
{
  num_avail_thrds_--;
}

}
