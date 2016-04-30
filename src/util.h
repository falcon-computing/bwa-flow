#ifndef BWAFLOW_UTIL_H
#define BWAFLOW_UTIL_H

#include <stdexcept>
#include <string>
#include <sys/syscall.h>
#include <sys/time.h>
#include <syscall.h>
#include <time.h>

#include "bwa/bwamem.h"

class resultsError : public std::runtime_error {
 public:
  explicit resultsError(const std::string& what_arg):
    std::runtime_error(what_arg) {;}
};

void regionsCompare(
    mem_alnreg_v *alnreg_base,
    mem_alnreg_v *alnreg_test,
    int num_seqs);

void smemCompare(
    smem_aux_t **smem_base,
    smem_aux_t **smem_test,
    int num
);

inline uint64_t getUs() {
  struct timespec tr;
  clock_gettime(CLOCK_REALTIME, &tr);

  return (uint64_t)tr.tv_sec*1e6 + tr.tv_nsec/1e3;
}

#endif
