#include "Common.h"

namespace kestrelFlow {

  uint64_t getUs() {
    struct timespec tr;
    clock_gettime(CLOCK_REALTIME, &tr);

    return (uint64_t)tr.tv_sec*1e6 + tr.tv_nsec/1e3;
  }

  uint64_t getMs() {
    struct timespec tr;
    clock_gettime(CLOCK_REALTIME, &tr);

    return (uint64_t)tr.tv_sec*1e3 + tr.tv_nsec/1e6;
  }
} // namespace kestrelFlow


