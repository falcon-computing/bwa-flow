#ifndef BWA_FLOW_FPGAPIPELINE_H
#define BWA_FLOW_FPGAPIPELINE_H

#include "bwa_wrapper.h"
#include "Pipeline.h"
#include "SWTask.h"
#include "SMemTask.h"

#include "kflow/MapStage.h"
#include "kflow/MapPartitionStage.h"

template<typename T, typename R>
class TimeoutMonitor
{
 public:
  TimeoutMonitor(int n, kestrelFlow::StageBase *stage) : n_(n), timecard(n), status(n), tasklist(2*n), records(n), stage_(stage) {;}
  ~TimeoutMonitor() { t_.join(); };

  std::vector<boost::atomic<uint64_t>>   timecard;
  std::vector<boost::atomic<int>>        status;
  std::vector<T*>                        tasklist;
  std::vector<R>                         records;

  kestrelFlow::StageBase* stage_;

  boost::thread t_;
  int n_;

  void timeout_func() {
    for (int i = 0; i < n_; i++) {
      status[i].store(0);
    }
    while (true) {
      boost::this_thread::sleep_for(boost::chrono::seconds(60));
      bool is_finished = true;
      for (int i = 0; i < n_; i++) {
        if (status[i].load() == 1) {
          uint64_t cur_time = getUs()/1000;
          if (cur_time > 60000 + timecard[i].load()) {
            if (status[i].load() == 1) {
              stage_->getThread(i)->interrupt(); 
              DLOG(WARNING) << "Monitor interrupted one thread due to timeout";
              status[i].store(-1);
            }
          }
        }
        if (status[i].load() != -1) {
          is_finished = false;
        }
      } 
      if (is_finished) break;
    }
  }

  void kick_off() {
    //t_ = boost::thread(boost::bind(this, TimeoutMonitor<T, R>::timeout_func));
    DLOG(INFO) << "Started one timeout monitor";
    t_ = boost::thread(boost::bind(&TimeoutMonitor::timeout_func, this));
  }
};

class ChainsToRegionsFPGA
: public kestrelFlow::MapPartitionStage<
      ChainsRecord, RegionsRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>
{
 public:
  ChainsToRegionsFPGA(int n=1): kestrelFlow::MapPartitionStage<
      ChainsRecord, RegionsRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>(n, false),
      timeout_(n, static_cast<StageBase*>(this))
  {
    n_active_ = n;
    if (n > 0) {
      timeout_.kick_off();
    }
  }

  void compute(int wid);
  void compute_func(int wid);
  void processOutput(SWTask* task, uint64_t &deq_ts, uint64_t &read_ts, uint64_t &post_ts);

  ChainsToRegions * cpu_stage;

 private:
  TimeoutMonitor<SWTask, ChainsRecord> timeout_;
  boost::mutex mtx_;
  int n_active_;
};

class SeqsToChainsFPGA
: public kestrelFlow::MapPartitionStage<
      SeqsRecord, ChainsRecord, INPUT_DEPTH, COMPUTE_DEPTH> {
 public:
  SeqsToChainsFPGA(int n=1) : kestrelFlow::MapPartitionStage<
      SeqsRecord, ChainsRecord, INPUT_DEPTH, COMPUTE_DEPTH>(n, false),
      timeout_(n, static_cast<kestrelFlow::StageBase*>(this))
  {
    n_active_ = n;
    if (n > 0) {
      timeout_.kick_off();
    }
  }

  //ChainsRecord compute(SeqsRecord const & record);
  void compute(int wid);
  void compute_func(int wid);

  SeqsToChains * cpu_stage;

 private:
  TimeoutMonitor<SMemTask, SeqsRecord> timeout_;
  boost::mutex mtx_;
  int n_active_;
};

#endif

