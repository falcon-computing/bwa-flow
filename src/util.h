#ifndef BWAFLOW_UTIL_H
#define BWAFLOW_UTIL_H

#include <stdexcept>
#include <string>
#include <sys/syscall.h>
#include <sys/time.h>
#include <syscall.h>
#include <time.h>

#include "bwa/bwamem.h"
#include "bwa_wrapper.h"
#include "Pipeline.h"  

class resultsError : public std::runtime_error {
 public:
  explicit resultsError(const std::string& what_arg):
    std::runtime_error(what_arg) {;}
};

class fpgaResultsError: public std::runtime_error {
 public:
  explicit fpgaResultsError(const std::string& what_arg):
    std::runtime_error(what_arg) {;}
};

class fpgaHangError : public std::runtime_error {
 public:
  explicit fpgaHangError(const std::string& what_arg):
    std::runtime_error(what_arg) {;}
};

inline uint64_t getUs() {
  struct timespec tr;
  clock_gettime(CLOCK_REALTIME, &tr);

  return (uint64_t)tr.tv_sec*1e6 + tr.tv_nsec/1e3;
}

inline uint64_t getNs() {
  struct timespec tr;
  clock_gettime(CLOCK_REALTIME, &tr);

  return (uint64_t)tr.tv_sec*1e9 + tr.tv_nsec;
}

void serialize(std::stringstream &ss, bseq1_t& seq);
void serialize(std::stringstream &ss, mem_chain_v& chains);
void serialize(std::stringstream &ss, mem_chain_t& chain);
void serialize(std::stringstream &ss, mem_seed_t& seed);
void serialize(std::stringstream &ss, mem_alnreg_v& alnregs);
void serialize(std::stringstream &ss, mem_alnreg_t& alnreg);

void deserialize(std::stringstream &ss, bseq1_t& seq);
void deserialize(std::stringstream &ss, mem_chain_v& chains);
void deserialize(std::stringstream &ss, mem_chain_t& chain);
void deserialize(std::stringstream &ss, mem_alnreg_t& alnreg);
void deserialize(std::stringstream &ss, mem_alnreg_v& alnregs);

std::string serialize(SeqsRecord& record);
std::string serialize(ChainsRecord& record);
std::string serialize(RegionsRecord& record);

void deserialize(const char* data, size_t length, SeqsRecord& record);
void deserialize(const char* data, size_t length, ChainsRecord& record);
void deserialize(const char* data, size_t length, RegionsRecord& record);

#endif
