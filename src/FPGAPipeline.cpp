#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <vector>
#include <queue>
#include <list>
#include <boost/asio.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>

#include "bwa_wrapper.h"
#include "config.h"
#include "FPGAAgent.h"
#include "util.h"
#include "Pipeline.h"
#include "FPGAPipeline.h"

#define MAX_BAND_TRY  2
//#define SMITHWATERMAN_SIM

#ifdef SMITHWATERMAN_SIM
// hw data structures
extern "C"{
void sw_top (int *a, int *output_a, int __inc);
}
#endif

void ChainsToRegionsFPGA::extendOnFPGA(
    FPGAAgent* agent,
    char* &kernel_buffer,
    int data_size_a,
    int data_size_b,
    int stage_cnt
    )
{
  agent->writeInput((int*)kernel_buffer, data_size_a, stage_cnt, 0);
  agent->writeInput((int*)(&kernel_buffer[data_size_a]), data_size_b, stage_cnt, 1);
  agent->start(data_size_a/4, data_size_b/4, stage_cnt);
}

void ChainsToRegionsFPGA::FPGAPostProcess(
    FPGAAgent* agent,
    short* kernel_output,
    int task_num_a,
    int task_num_b,
    mem_alnreg_t** &region_batch,
    mem_chain_t** &chain_batch,
    int stage_cnt
    )
{
  uint64_t start_ts = getUs();
  int seedcov = 0;
  agent->wait(stage_cnt);
  DLOG_IF(INFO, VLOG_IS_ON(3)) << "Wait for FPGA takes " 
    << getUs() - start_ts << " us";

  start_ts = getUs();
  agent->readOutput(kernel_output, FPGA_RET_PARAM_NUM*task_num_a*4, stage_cnt, 0);
  agent->readOutput(&kernel_output[FPGA_RET_PARAM_NUM*task_num_a*2],
      FPGA_RET_PARAM_NUM*task_num_b*4, stage_cnt, 1);
  for (int i = 0; i < task_num_a + task_num_b; i++) {
    int seed_idx = ((int)(kernel_output[1+FPGA_RET_PARAM_NUM*2*i])<<16) |
                    kernel_output[0+FPGA_RET_PARAM_NUM*2*i];
    mem_alnreg_t *newreg = region_batch[seed_idx];

    newreg->qb = kernel_output[2+FPGA_RET_PARAM_NUM*2*i]; 
    newreg->qe += kernel_output[3+FPGA_RET_PARAM_NUM*2*i];
    newreg->rb += kernel_output[4+FPGA_RET_PARAM_NUM*2*i];
    newreg->re += kernel_output[5+FPGA_RET_PARAM_NUM*2*i];
    newreg->score = kernel_output[6+FPGA_RET_PARAM_NUM*2*i]; 
    newreg->truesc = kernel_output[7+FPGA_RET_PARAM_NUM*2*i]; 
    newreg->w = kernel_output[8+FPGA_RET_PARAM_NUM*2*i];

    // compute the seedcov here
    mem_chain_t *chain = chain_batch[seed_idx];
    seedcov = 0;
    for (int k=0; k < chain->n; k++) {
      const mem_seed_t *seed = &chain->seeds[k];
      if (seed->qbeg >= newreg->qb && 
          seed->qbeg + seed->len <= newreg->qe && 
          seed->rbeg >= newreg->rb && 
          seed->rbeg + seed->len <= newreg->re){
        seedcov += seed->len; 
      }
    }
    newreg->seedcov = seedcov;
  }
  DLOG_IF(INFO, VLOG_IS_ON(3)) << "Process output takes " 
    << getUs() - start_ts << " us";
}

inline int cal_max_gap(
    const mem_opt_t *opt, 
    int qlen
) {
	int l_del = (int)((double)(qlen * opt->a - opt->o_del) / opt->e_del + 1.);
	int l_ins = (int)((double)(qlen * opt->a - opt->o_ins) / opt->e_ins + 1.);
	int l = l_del > l_ins? l_del : l_ins;
	l = l > 1? l : 1;
	return l < opt->w<<1? l : opt->w<<1;
}

inline void getChainRef(
      const ktp_aux_t* aux,
      const bseq1_t* seq,
      const mem_chain_v* chain,
      mem_chainref_t* &ref
) {
  int64_t l_pac = aux->idx->bns->l_pac;

  // input for SmithWaterman for each chain
  ref = (mem_chainref_t*)malloc(chain->n*sizeof(mem_chainref_t));

  for (int j = 0; j < chain->n; j++) {
    // Prepare the maxspan and rseq for each seed
    mem_chain_t *c = &chain->a[j]; 

    int64_t rmax[2], tmp, max = 0;
    uint8_t *rseq = 0;
    uint64_t *srt;

    if (c->n == 0) {
      continue;
    }
    // get the max possible span
    rmax[0] = l_pac<<1; rmax[1] = 0;
    for (int k = 0; k < c->n; ++k) {
      int64_t b, e;
      const mem_seed_t *t = &c->seeds[k];
      b = t->rbeg - (t->qbeg + cal_max_gap(aux->opt, t->qbeg));
      e = t->rbeg + t->len + ((seq->l_seq - t->qbeg - t->len)
          + cal_max_gap(aux->opt, seq->l_seq - t->qbeg - t->len));
      rmax[0] = rmax[0] < b? rmax[0] : b;
      rmax[1] = rmax[1] > e? rmax[1] : e;
      if (t->len > max) max = t->len;
    }
    rmax[0] = rmax[0] > 0? rmax[0] : 0;
    rmax[1] = rmax[1] < l_pac<<1? rmax[1] : l_pac<<1;
    if (rmax[0] < l_pac && l_pac < rmax[1]) { // crossing the forward-reverse boundary; then choose one side
      if (c->seeds[0].rbeg < l_pac) rmax[1] = l_pac; // this works because all seeds are guaranteed to be on the same strand
      else rmax[0] = l_pac;
    }
    // retrieve the reference sequence
    int rid;
    rseq = bns_fetch_seq(aux->idx->bns, aux->idx->pac, &rmax[0], c->seeds[0].rbeg, &rmax[1], &rid);
    assert(c->rid == rid);

    srt = (uint64_t *)malloc(c->n * 8);
    for (int l = 0; l < c->n; ++l)
      srt[l] = (uint64_t)c->seeds[l].score<<32 | l;
    ks_introsort_64(c->n, srt);

    ref[j].rmax[0] = rmax[0];
    ref[j].rmax[1] = rmax[1];
    ref[j].rseq = rseq;
    ref[j].srt = srt;
  }
}

inline void packReadData(ktp_aux_t* aux, 
    const bseq1_t* seq, 
    const mem_chain_v* chains,
    mem_alnreg_v* alnregs, 
    char* buffer, 
    int &buffer_idx, 
    int &task_num, 
    mem_alnreg_t** &region_batch,
    mem_chain_t** &chain_batch) 
{
  mem_chainref_t* ref;
  int chain_idx = 0;
  int seed_idx = chains->a[0].n - 1;
  getChainRef(aux, seq, chains, ref);

  for ( ; chain_idx < chains->n; 
          seed_idx > 0 ? seed_idx-- : seed_idx = chains->a[++chain_idx].n-1)
  {
      uint32_t sorted_idx = (uint32_t)(ref[chain_idx].srt[seed_idx]);

    // get next available seed in the current read
      mem_seed_t* seed_array = &chains->a[chain_idx].seeds[sorted_idx];
      // initialize the newreg
      mem_alnreg_t* newreg = kv_pushp(mem_alnreg_t, *alnregs);
      memset(newreg, 0, sizeof(mem_alnreg_t));

      newreg->score  = seed_array->len * aux->opt->a;
      newreg->truesc = seed_array->len * aux->opt->a;
      newreg->qb = seed_array->qbeg;
      newreg->rb = seed_array->rbeg;
      newreg->qe = seed_array->qbeg + seed_array->len;
      newreg->re = seed_array->rbeg + seed_array->len; 
      newreg->rid = chains->a[chain_idx].rid;
      newreg->seedlen0 = seed_array->len; 
      newreg->frac_rep = chains->a[chain_idx].frac_rep;
      newreg->w = aux->opt->w;
  }
  chain_idx = 0;
  seed_idx = 0;
  int counter8 = 0;
  int tmp_int = 0;
  int chain_num = 0;
  int seed_num = 0;
  int chain_num_addr = 0;
  int seed_num_addr = 0;
  int idx_end_addr = buffer_idx;
  buffer_idx += 4;

  // pack the read sequence
  *((int*)(&buffer[buffer_idx])) = seq->l_seq;
  buffer_idx += 4;
  for ( int i = 0; i < seq->l_seq; i++ ){
    counter8 += 1;
    tmp_int = tmp_int << 4 | seq->seq[i];
    if ( counter8 % 8 ==0 ) {
      *(uint32_t*) (&buffer[buffer_idx]) = tmp_int ;
      buffer_idx += 4 ;
    }
  }
  if ( counter8 % 8 !=0 ) {
    tmp_int = tmp_int << (4*(8 - counter8 % 8));
    *(uint32_t*) (&buffer[buffer_idx]) = tmp_int ;
    buffer_idx += 4 ;
  }
  counter8 = 0;
  tmp_int = 0;
  chain_num_addr = buffer_idx;
  buffer_idx += 4;
  int region_num = 0; // used to keep track of the current region

  for ( chain_idx = 0; chain_idx < chains->n; chain_idx++) {
    // Pack the maxspan and rseq
    *((int64_t*)(&buffer[buffer_idx]))= ref[chain_idx].rmax[0];
    buffer_idx += 8;
    *((int64_t*)(&buffer[buffer_idx]))= ref[chain_idx].rmax[1];
    buffer_idx += 8;
    for ( int i = 0; i < ref[chain_idx].rmax[1] - ref[chain_idx].rmax[0]; i++) {
      counter8 = counter8 + 1;
      tmp_int = tmp_int << 4 | ref[chain_idx].rseq[i];
      if ( counter8 % 8 ==0 ) {
        *(uint32_t*) (&buffer[buffer_idx]) = tmp_int ;
        buffer_idx += 4 ;
      }
    } 
    if ( counter8 % 8 !=0 ) {
      tmp_int = tmp_int << (4*(8 - counter8 % 8));
      *(uint32_t*) (&buffer[buffer_idx]) = tmp_int ;
      buffer_idx += 4 ;
    }
    counter8 = 0;
    // record the address of seed number
    seed_num_addr = buffer_idx ;
    seed_num = 0;
    buffer_idx += 4 ; 
    // pack the seed information
    for ( seed_idx = chains->a[chain_idx].n -1 ; seed_idx >= 0 ; seed_idx--) {
      uint32_t sorted_idx = (uint32_t)(ref[chain_idx].srt[seed_idx]);
      // get next available seed in the current read
      mem_seed_t* seed_array = &chains->a[chain_idx].seeds[sorted_idx];
     
      if (seed_array->qbeg > 0 ||
          seed_array->qbeg + seed_array->len != seq->l_seq)
      {
        region_batch[task_num] = &(alnregs->a[region_num]) ;
        chain_batch[task_num] = &(chains->a[chain_idx]) ;
        *((int*)(&buffer[buffer_idx])) = task_num ; 
        buffer_idx += 4;
        task_num += 1;
        seed_num += 1;
        *((int64_t*)(&buffer[buffer_idx])) = seed_array->rbeg ;
        buffer_idx += 8;
        *((int32_t*)(&buffer[buffer_idx])) = seed_array->qbeg ;
        buffer_idx += 4;
        *((int32_t*)(&buffer[buffer_idx])) = seed_array->len ;
        buffer_idx += 4;
      }
      else {
        // no need to loop again 
        alnregs->a[region_num].seedcov = seed_array->len;
      }
      region_num++;
    }
     *((int*)(&buffer[seed_num_addr])) = seed_num ;
     chain_num += 1;
  }
  *((int*)(&buffer[chain_num_addr])) = chain_num ;
  *((int*)(&buffer[idx_end_addr])) = buffer_idx/4 ;

  for (int i=0; i< chains->n; i++) {
    free(ref[i].srt);
    free(ref[i].rseq);
  }
  free(ref);
}

void ChainsToRegionsFPGA::compute(int wid) {

  int chunk_size = FLAGS_chunk_size;

  const int stage_num = 2;
  int stage_cnt = 0;

  // Create FPGAAgent
  FPGAAgent agent(opencl_env, chunk_size);

  int task_num = 0;

  // Batch of regions
  mem_alnreg_t** region_batch[stage_num];
  mem_chain_t** chain_batch[stage_num];
  int stage_task_num_a[stage_num];
  int stage_task_num_b[stage_num];
  for (int i = 0; i < stage_num; i++) {
    region_batch[i] = new mem_alnreg_t*[2*chunk_size];
    chain_batch[i] = new mem_chain_t*[2*chunk_size];
    stage_task_num_a[i] = 0;
    stage_task_num_b[i] = 0;
  } 
  // elements of the Regions record
  bseq1_t* seqs        ;
  mem_chain_v* chains  ;
  mem_alnreg_v* alnreg ;

  // kernel_buffer 
  char* kernel_buffer = new char [10000*chunk_size];
  short* kernel_buffer_out = new short [40 * chunk_size];
  int kernel_buffer_idx = 0;
  // For statistics
  uint64_t start_ts;
  uint64_t last_batch_ts;
  uint64_t last_output_ts;
  uint64_t last_read_ts;

  uint64_t process_time  = 0;
  uint64_t pending_time  = 0;
  uint64_t finish_time   = 0;
  uint64_t nextTask_time = 0;
  uint64_t batch_time    = 0;
  uint64_t output_time   = 0;

  uint64_t nextTask_num = 0;
  int      start_idx    = 0;
  int      batch_num    = 0;

  bool flag_need_reads = false;
  bool flag_more_reads = true;
  ChainsRecord record; 
  while (flag_more_reads) { 
    // get one Chains record
    bool ready = this->getInput(record);
    start_ts = getUs();
    while (!this->isFinal() && !ready) {
      boost::this_thread::sleep_for(boost::chrono::microseconds(100));
      ready = this->getInput(record);
    }
    DLOG_IF(INFO, VLOG_IS_ON(2)) << "Wait for input for FPGA takes " << getUs() - start_ts << " us";
    if(!ready) {
      flag_more_reads = false;
      break;
    }
    DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started ChainsToRegions() on FPGA ";

    last_batch_ts  = getUs();
    last_output_ts = last_batch_ts;
    last_read_ts  = last_batch_ts;
    start_idx = record.start_idx;
    batch_num = record.batch_num;
    seqs = record.seqs;
    chains = record.chains;
    alnreg = record.alnreg;

    int i = 0;
    uint64_t last_prepare_ts = getUs();
    uint64_t last_batch_ts = getUs();
    int kernel_task_num_a = 0;
    int kernel_task_num_b = 0;
    int kernel_size_a = 0;
    int kernel_size_b = 0;
    bool reach_half = false;
    while (i < batch_num) {
      if (task_num < chunk_size/2) {
        packReadData(aux, seqs+i, chains+i,
            alnreg+i, kernel_buffer, kernel_buffer_idx, task_num,
            region_batch[stage_cnt], chain_batch[stage_cnt]);
        i++;
      }
      else if (task_num >= chunk_size/2 && reach_half == false) {
        kernel_size_a = kernel_buffer_idx;
        kernel_task_num_a = task_num;
        kernel_buffer_idx = 0;
        reach_half = true;
      }
      else if (task_num < chunk_size) {
        packReadData(aux, seqs+i, chains+i,
            alnreg+i, &kernel_buffer[kernel_size_a], kernel_buffer_idx, task_num,
            region_batch[stage_cnt], chain_batch[stage_cnt]);
        i++;
      }
      else if (task_num >= chunk_size) {
        kernel_size_b = kernel_buffer_idx;
        kernel_task_num_b = task_num - kernel_task_num_a;
        DLOG_IF(INFO, VLOG_IS_ON(3)) << "Prepare data for FPGA  "<< stage_cnt << " takes " << getUs()- last_prepare_ts << " us " ;
        stage_task_num_a[stage_cnt] = kernel_task_num_a ;
        stage_task_num_b[stage_cnt] = kernel_task_num_b ;
        //extendOnFPGA(&agent, kernel_buffer, kernel_size_a, kernel_size_b, stage_cnt);
        //stage_cnt = 1 - stage_cnt;
        //if (agent.pending(stage_cnt)){
        //  FPGAPostProcess(&agent, kernel_buffer_out, stage_task_num_a[stage_cnt], 
        //      stage_task_num_b[stage_cnt],region_batch[stage_cnt], chain_batch[stage_cnt], stage_cnt );
        //}
        stage_cnt = 1 - stage_cnt;
        if (agent.pending(stage_cnt)){
          FPGAPostProcess(&agent, kernel_buffer_out, stage_task_num_a[stage_cnt], 
              stage_task_num_b[stage_cnt],region_batch[stage_cnt], chain_batch[stage_cnt], stage_cnt );
        }
        extendOnFPGA(&agent, kernel_buffer, kernel_size_a, kernel_size_b, 1 - stage_cnt);
        last_prepare_ts = getUs();
        DLOG_IF(INFO, VLOG_IS_ON(3)) << "This chunk of FPGA takes " 
          << getUs()- last_batch_ts<< " us " ;
        last_batch_ts = getUs();
        task_num = 0;
        kernel_buffer_idx = 0;
        kernel_size_a = 0;
        kernel_size_b = 0;
        kernel_task_num_a = 0;
        kernel_task_num_b = 0;
        reach_half = false;
      }
    }
    if (agent.pending(1 - stage_cnt)){
      FPGAPostProcess(&agent, kernel_buffer_out, stage_task_num_a[1 - stage_cnt], 
          stage_task_num_b[1 - stage_cnt], region_batch[1 - stage_cnt],
          chain_batch[1 - stage_cnt], 1 - stage_cnt );
    }
    // finish the remain reads even with small task number 
    extendOnFPGA(&agent, kernel_buffer, kernel_size_a, kernel_size_b, stage_cnt);
    FPGAPostProcess(&agent, kernel_buffer_out, kernel_task_num_a, kernel_task_num_b,
        region_batch[stage_cnt], chain_batch[stage_cnt], stage_cnt );
    task_num = 0;
    kernel_buffer_idx = 0;
    kernel_size_a = 0;
    kernel_size_b = 0;
    kernel_task_num_a = 0;
    kernel_task_num_b = 0;
    reach_half = false;
    i = 0;

    RegionsRecord outputRecord;
    outputRecord.start_idx = start_idx;
    outputRecord.batch_num = batch_num;
    outputRecord.seqs = seqs;
    outputRecord.chains = chains;
    outputRecord.alnreg = alnreg;

    freeChains(chains, batch_num);
    pushOutput(outputRecord);
    DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished ChainsToRegions() on FPGA for "
      << getUs() - last_output_ts << " us";
    last_output_ts = getUs();
  }
}
