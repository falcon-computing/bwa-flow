#ifndef BWAFLOW_UTIL_H
#define BWAFLOW_UTIL_H

#include <stdexcept>
#include <string>

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

#endif
