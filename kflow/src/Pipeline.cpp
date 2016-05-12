#include "Pipeline.h"

namespace kestrelFlow {

Pipeline::Pipeline(int num_stages, int num_threads): 
  num_stages_(num_stages),
  num_threads_(num_threads),
  stages_(num_stages, NULL),
  queues_(num_stages+1, NULL_QUEUE_PTR),
  num_active_threads_(0) 
{
  boost::shared_ptr<boost::asio::io_service> ios(new boost::asio::io_service);
  ios_ = ios;

  // Start io service processing loop
  boost::shared_ptr<boost::asio::io_service::work> work(
      new boost::asio::io_service::work(*ios));

  ios_work_ = work;

  for (int t = 0; t < num_threads; t++) {
    workers_.create_thread(
        boost::bind(&boost::asio::io_service::run, ios.get()));
  }
  DLOG(INFO) << "Started " << workers_.size() << " Pipeline workers.";
}

Pipeline::~Pipeline() {
  // Allow io_service::run to exit
  ios_work_.reset();
  workers_.join_all();
}

void Pipeline::finalize() {
  stages_[0]->final();
  DLOG(INFO) << "Finalized first stage";
}

boost::any Pipeline::getConst(std::string key) {
  if (!constants_.count(key)) {
    return 0;
  }
  else {
    return constants_[key];
  }
}

void Pipeline::start() {

  bool initialized = true;
  for (int i = 0; i < num_stages_; i++) {
    if (!stages_[i]) {
      LOG(ERROR) << "Stage " << i << " is uninitialized";
      initialized = false;
    }
  }
  if (!initialized) {
    LOG(ERROR) << "Cannot start the pipeline due to previous errors";
  }

  for (int i = 0; i < num_stages_; i++) {
    StageBase* stage = stages_[i];

    // Start only stages that require dynamic task dispatch
    if (!stage->isDynamic()) {
      stage->start();
      DLOG(INFO) << "Start workers for stage " << i;
    }
    else {
      pending_stages_.push_back(stage);
    }
  }
  // Start scheduler
  scheduler_.create_thread(boost::bind(&Pipeline::schedule, this));
}

void Pipeline::stop() {

  // Stop unschedulable stages
  for (int i = 0; i < num_stages_; i++) {
    if (!stages_[i]) {
      LOG(WARNING) << "Stage " << i << " is uninitialized";
    }
    else {
      stages_[i]->stop();
    }
  }
  // Stop scheduler
  scheduler_.interrupt_all();
}

void Pipeline::wait() {
  // Only need to wait for the last stage
  if (!stages_[num_stages_-1]->isDynamic()) {
    stages_[num_stages_-1]->wait();
  }
  else {
    scheduler_.join_all(); 

    // Allow io_service::run to exit
    ios_work_.reset();
    workers_.join_all();
  }
}

void Pipeline::schedule() {
  
  DLOG(INFO) << "Started a scheduler";
  try {
    while (!pending_stages_.empty()) {

      bool dispatched = false;
      // Check if there is idle thread
      if (num_active_threads_.load() < num_threads_) {
        // Always try to start execute from later stages
        // if the corresponding input queue is not empty
        for (std::deque<StageBase*>::reverse_iterator 
             iter  = pending_stages_.rbegin();
             iter != pending_stages_.rend();
             iter ++) {
          if ((*iter)->execute()) {
            dispatched = true;
            break; 
          }
        }
      }
      if (!dispatched) {
        StageBase* stage = pending_stages_.front();
        if (stage->isFinal() && stage->getNumThreads()==0) {
          stage->finalize();
          pending_stages_.pop_front();

          DLOG(INFO) << "Front-of-queue stage is finished";
        }
        else {
          boost::this_thread::sleep_for(boost::chrono::milliseconds(10));
        }
      }
    }
    DLOG(INFO) << "Scheduler is finished";
  }
  catch (boost::thread_interrupted &e) {
    // TODO: stop io_service and join all threads
    VLOG(2) << "Scheduler is interrupted";
  }
}

QueueBase* Pipeline::getInputQueue() {
  return queues_[0].get();
}

QueueBase* Pipeline::getOutputQueue() {
  return queues_[num_stages_].get();
}

} // namepspace kestrelFlow
