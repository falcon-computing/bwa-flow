#ifndef KFLOW_OCCUPANCY_H
#define KFLOW_OCCUPANCY_H

#include <boost/atomic.hpp>
#include <boost/type_traits.hpp>
#include <typeinfo>

#include "Common.h"
#include "Pipeline.h"
#include "Stage.h"

// for testing purpose
#ifndef TEST_FRIENDS_LIST
#define TEST_FRIENDS_LIST
#endif

namespace kestrelFlow {

class OccupancyCounter {
 public:
  OccupancyCounter(Pipeline* p_host, StageBase* s_host):
    p_host_(p_host), s_host_(s_host) 
  {
    // p_host_->addSeat(); 
    s_host_->addSeat(); 
  }
  ~OccupancyCounter() {
    // p_host_->removeSeat();
    s_host_->removeSeat();
  }
 private:
  Pipeline*  p_host_;
  StageBase* s_host_;
};

} // namespace kestrelFlow
#endif
