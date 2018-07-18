#ifndef KFLOW_MEGAPIPE_H
#define KFLOW_MEGAPIPE_H

#include <boost/any.hpp>
#include <boost/asio.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
#include <deque>

#include "Common.h"
#include "Pipeline.h"

// for testing purpose
#ifndef TEST_FRIENDS_LIST
#define TEST_FRIENDS_LIST
#endif

namespace kestrelFlow {

#if 0
class StageHandler {
  friend class MegaPipe;
  TEST_FRIENDS_LIST

 public:
  StageHandler(StageBase *cur_stage, 
               std::vector<StageBase *> &downstream_stages_list,
               int num_downnstream_stages,
               bool is_dynamic = true,
               bool use_accx = false,
               bool is_final = false )
  : stage_(cur_stage),
    num_downnstream_stages_(num_downnstream_stages),
    is_dynamic_(is_dynamic),
    use_accx_(use_accx),
    is_final_(is_final)
  {
    num_downnstream_stages_ = std::min(num_downnstream_stages, downstream_stages_list.size() );
    for (int i = 0 ; i < num_downnstream_stages_; i++ )
      downstream_stages_list_.push_back(downstream_stages_list[i]);
  }

 private:
  StageBase *stage_;
  std::vector<StageBase *> downstream_stages_list_;
  int num_downnstream_stages_;
  bool is_dynamic_;
  bool use_accx_;
  bool is_final_;
}
#endif

class MegaPipe {
  friend class OccupancyCounter;
  TEST_FRIENDS_LIST

 public:
  MegaPipe(int num_threads, int num_accelerators = 0);
  ~MegaPipe();

  bool addPipeline(Pipeline *pipeline, int priority = 1);

  void start();
  void stop();
  void wait();
  void finalize();

  bool acqAccx();
  void relAccx();

  bool acqThrd();
  void relThrd();
  void getThrd();

 private:
  void dworker_func(int pipeline_id);

  boost::thread_group dworkers_;

  std::vector<Pipeline *> pipelines_; // priority descending
  std::vector<int>        priorities_;

  int num_threads_;
  int num_accelerators_;
  mutable boost::atomic<int> num_avail_thrds_;
  mutable boost::atomic<int> num_avail_accs_;

};

} // namespace kestrelFlow
#endif

