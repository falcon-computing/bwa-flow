#ifndef SWREAD_H
#define SWREAD_H

#include "bwa_wrapper.h"
#include "SWTask.h"

class SWRead {
 public:

  SWRead(int idx,
      ktp_aux_t* aux,
      const bseq1_t* seq, 
      const mem_chain_v* chains,
      mem_alnreg_v* alnregs);

  ~SWRead();

  enum TaskStatus {
    Successful = 0,
    Pending = 1,
    Finished = 2
  };

  // Get one task from chains of current read
  // return values:
  // - 0: successful
  // - 1: task pending for extension
  // - 2: no more extension to do for this read
  enum TaskStatus nextTask(ExtParam* &task);

  // Called after seed extension is done
  void finish();

  int index() { return read_idx_; }
 
 private:

  static inline ExtParam* getTask(
      mem_opt_t *opt,
      const mem_seed_t *seed, 
      const uint8_t *query,
      int l_query,
      int64_t rmax_0, 
      int64_t rmax_1,
      uint8_t *rseq,
      int idx);

  static inline void prepareChainRef(
      const ktp_aux_t* aux,
      const bseq1_t* seq,
      const mem_chain_v* chain,
      mem_chainref_t* &ref);

  static inline int testExtension(
      mem_opt_t *opt,
      mem_seed_t& seed,
      mem_alnreg_v& alnregv, 
      int l_query);

  static inline int checkOverlap(
      int startidx,
      mem_seed_t& seed,
      mem_chain_t& chain,
      uint64_t *srt);

  static inline int cal_max_gap(
      const mem_opt_t *opt, 
      int qlen);

  bool is_pend_;
  bool is_finished_;
  
  int read_idx_;
  int chain_idx_;
  int seed_idx_;

  ktp_aux_t*         aux_;
  const bseq1_t*     seq_;
  const mem_chain_v* chains_;
  mem_chainref_t*    ref_;
  mem_alnreg_v*      alnregs_;

  mem_alnreg_t*      newreg_;
};

#endif
