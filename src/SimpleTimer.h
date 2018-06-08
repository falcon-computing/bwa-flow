#ifndef SIMPLE_TIMER_H
#define SIMPLE_TIMER_H

#include <sys/time.h>
#include <time.h>
#include <unordered_map>

class SimpleTimer {
  public:
    SimpleTimer(std::string func = "", 
        bool debug = true,
        int verbose = 1): 
      func_(func),
      debug_(debug),
      verbose_(verbose) 
    {
      start(-1);
    }

    ~SimpleTimer() {
      if (!func_.empty()) {
        stop(-1);
        if (debug_) {
          DLOG_IF(INFO, VLOG_IS_ON(verbose_)) << 
            func_ << " takes " << get_ms(-1) << " ms";
        }
        else {
          VLOG(verbose_) << func_ << " takes " << get_ms(-1) << " ms";
        }
      }
    }

    uint64_t start(int idx = 0) {
      if (!elapsed_ns_.count(idx)) {
        elapsed_ns_[idx] = 0;
        laps_[idx] = 0;
      }
      if (!timer_on_.count(idx) || !timer_on_[idx]) {
        start_ns_[idx] = getNs();
        timer_on_[idx] = true;
      }
      return start_ns_[idx];
    }

    uint64_t stop(int idx = 0) {
      if (timer_on_.count(idx) && timer_on_[idx]) {
        elapsed_ns_[idx] += getNs() - start_ns_[idx];
        timer_on_[idx] = false;
        laps_[idx] += 1;
        return elapsed_ns_[idx];
      }
      else {
        timer_on_[idx] = false;
        laps_[idx] = 0;
        return 0;
      }
    }

    void reset(int idx) {
      start_ns_[idx] = getNs();
      elapsed_ns_[idx] = 0;
      timer_on_[idx] = false;
      laps_[idx] = 0;
    }

    void reset() {
      for (auto e_time : elapsed_ns_) {
        elapsed_ns_[e_time.first] = 0;
        start_ns_[e_time.first] = getNs();
        timer_on_[e_time.first] = false;
        laps_[e_time.first] = 0;
      }
    }

    uint64_t get_ns(int idx = 0) {
      if (!timer_on_.count(idx) || timer_on_[idx] || !elapsed_ns_.count(idx)) {
        return 0;
      }
      else return elapsed_ns_[idx];
    }

    double get_ms(int idx = 0) {
      return (double)get_ns(idx)/1e6;
    }

    int get_laps(int idx = 0) {
      if (!laps_.count(idx)) laps_[idx] = 0;
      return laps_[idx];
    }

  private:
    uint64_t getNs() {
      struct timespec tr;
      clock_gettime(CLOCK_REALTIME, &tr);
      return (uint64_t)tr.tv_sec*1e9 + tr.tv_nsec;
    }

    std::string func_;
    bool        debug_;
    int         verbose_;

    std::unordered_map<int, uint64_t> elapsed_ns_;
    std::unordered_map<int, uint64_t> start_ns_;
    std::unordered_map<int, bool>     timer_on_;
    std::unordered_map<int, int>      laps_;
};

#endif
