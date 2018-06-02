#include <assert.h>
#include <boost/asio.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
#include <limits.h>
#include <list>
#include <math.h>
#include <deque>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <vector>

#include "bwa_wrapper.h"
#include "config.h"
#include "util.h"
#include "Pipeline.h"
#include "FPGAPipeline.h"
#include "SWTask.h"
#include "SMemTask.h"

#define MAX_BAND_TRY  2

void ChainsToRegionsFPGA::processOutput(SWTask* task) {
  // finish task and get output buffers
  task->finish();

  mem_alnreg_t** region_batch = task->region_batch;
  mem_chain_t**  chain_batch = task->chain_batch;

  int total_task_num = (task->o_size[0]+task->o_size[1])/FPGA_RET_PARAM_NUM;
  int actual_tasks = 0;

  uint64_t start_ts = getUs();
  bool wrong_results = false;
  int wrong_half = -1;
  for (int k = 0; k < 2; k++) {
    int task_num = task->o_size[k] / FPGA_RET_PARAM_NUM;
    short* kernel_output = task->o_data[k];

    //DLOG_IF(INFO, VLOG_IS_ON(3)) << "output " << k << " task num = " << task_num;
    for (int i = 0; i < task_num; i++) {
      int seed_idx = ((int)(kernel_output[1+FPGA_RET_PARAM_NUM*2*i])<<16) |
        kernel_output[0+FPGA_RET_PARAM_NUM*2*i];

      if (seed_idx > actual_tasks) {
        actual_tasks = seed_idx;
      }
      mem_alnreg_t *newreg = region_batch[seed_idx];

      if (seed_idx > total_task_num || seed_idx < 0) {
        DLOG(ERROR) << "task_num = " << i << " "
                    << "seed_idx = " << seed_idx << " ";
        wrong_results = true;
        wrong_half = k;
        goto error;
      }

      newreg->qb = kernel_output[2+FPGA_RET_PARAM_NUM*2*i]; 
      newreg->qe += kernel_output[3+FPGA_RET_PARAM_NUM*2*i];
      newreg->rb += kernel_output[4+FPGA_RET_PARAM_NUM*2*i];
      newreg->re += kernel_output[5+FPGA_RET_PARAM_NUM*2*i];
      newreg->score = kernel_output[6+FPGA_RET_PARAM_NUM*2*i]; 
      newreg->truesc = kernel_output[7+FPGA_RET_PARAM_NUM*2*i]; 
      newreg->w = kernel_output[8+FPGA_RET_PARAM_NUM*2*i];

      // compute the seedcov here
      mem_chain_t *chain = chain_batch[seed_idx];
      int seedcov = 0;
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
    // reset
    task->i_size[k] = 0;
    task->o_size[k] = 0;
  }

  if (total_task_num != actual_tasks + 1) {
    DLOG(ERROR) << "problem with seedindex";
    wrong_results = true;
  }
error: 
  if (wrong_results) {
    // dump results   
    FILE* fout = fopen("dump-input.dat", "wb+");
    int k = wrong_half;
    fwrite(task->i_data[k], sizeof(int), task->i_size[k], fout);
    fclose(fout);
    DLOG(ERROR) << "dump input data of output-" << k <<  " to dump-input.dat";
    throw std::runtime_error("wrong fpga results");
  }

  // reset
  for (int k = 0; k < 2; k++) {
    task->i_size[k] = 0;
    task->o_size[k] = 0;
  }

  DLOG_IF(INFO, VLOG_IS_ON(3)) << "Total " << actual_tasks << " tasks in output";
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
    rseq = bns_fetch_seq_fpga(aux->idx->bns, aux->idx->pac, &rmax[0], c->seeds[0].rbeg, &rmax[1], &rid);

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
    SWTask* task) 
{
  mem_alnreg_t** region_batch = task->region_batch;
  mem_chain_t** chain_batch = task->chain_batch;
  int chunk_size = FLAGS_chunk_size;
  mem_chainref_t* ref;
  int seed_idx = chains->a[0].n - 1;
  getChainRef(aux, seq, chains, ref);

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
    if (counter8 % 8 ==0 ) {
      *(uint32_t*) (&buffer[buffer_idx]) = tmp_int ;
      buffer_idx += 4 ;
    }
  }
  if (counter8 % 8 !=0 ) {
    tmp_int = tmp_int << (4*(8 - counter8 % 8));
    *(uint32_t*) (&buffer[buffer_idx]) = tmp_int ;
    buffer_idx += 4 ;
  }
  counter8 = 0;
  tmp_int = 0;
  chain_num_addr = buffer_idx;
  buffer_idx += 4;
  int region_num = 0; // used to keep track of the current region

  for (int chain_idx = 0; chain_idx < chains->n; chain_idx++) {
    // Pack the maxspan and rseq
    *((int64_t*)(&buffer[buffer_idx]))= ref[chain_idx].rmax[0];
    buffer_idx += 8;
    *((int64_t*)(&buffer[buffer_idx]))= ref[chain_idx].rmax[1];
    buffer_idx += 8;

    // record the address of seed number
    seed_num_addr = buffer_idx ;
    seed_num = 0;
    buffer_idx += 4 ; 

    // pack the seed information
    for (int seed_idx = chains->a[chain_idx].n -1 ; seed_idx >= 0 ; seed_idx--) {
      mem_seed_t* seed_array = &chains->a[chain_idx].seeds[seed_idx];
      if (seed_array->qbeg > 0 ||
          seed_array->qbeg + seed_array->len != seq->l_seq) {

        region_batch[task_num] = &(alnregs->a[region_num]) ;
        chain_batch[task_num] = &(chains->a[chain_idx]) ;
        *((int*)(&buffer[buffer_idx])) = task_num ; 
        buffer_idx += 4;
        task_num += 1;
        if (task_num >= task->max_o_size_/(2*FPGA_RET_PARAM_NUM)) {
          task->max_o_size_ = task->max_o_size_ + 1000;
          task->o_data[0] = (short*)realloc(task->o_data[0], task->max_o_size_ * sizeof(short));
          task->o_data[1] = (short*)realloc(task->o_data[1], task->max_o_size_ * sizeof(short));
        }
        seed_num += 1;
        *((int64_t*)(&buffer[buffer_idx])) = seed_array->rbeg ;
        buffer_idx += 8;
        *((int32_t*)(&buffer[buffer_idx])) = seed_array->qbeg ;
        buffer_idx += 4;
        *((int32_t*)(&buffer[buffer_idx])) = seed_array->len ;
        buffer_idx += 4;
      }
      else {
        // need to loop again 
        int seedcov = 0;
        for (int i =0; i < chains->a[chain_idx].n; i++ ) {
          mem_seed_t *t = &chains->a[chain_idx].seeds[i];
          if (t->qbeg >=0 && t->qbeg + t->len <=seq->l_seq &&
              t->rbeg >= seed_array->rbeg && t->rbeg + t->len <= seed_array->rbeg + seed_array->len) {
            seedcov += t->len;
          }
        }
        alnregs->a[region_num].seedcov = seedcov;
      }
      region_num++;
    }
     *((int*)(&buffer[seed_num_addr])) = seed_num ;
     chain_num += 1;
  }
  *((int*)(&buffer[chain_num_addr])) = chain_num ;
  *((int*)(&buffer[idx_end_addr])) = buffer_idx/4 ;
  if (buffer_idx >= task->max_i_size_ * 0.9) {
    task->max_i_size_ = (int)(task->max_i_size_ * 1.2);
    task->i_data[0] = (char*)realloc(task->i_data[0], task->max_i_size_);
  }

  free(ref);
}

void ChainsToRegionsFPGA::compute(int wid) {
  int chunk_size = FLAGS_chunk_size;

  int stage_cnt = 0;

  // Create SWTasks
  std::deque<SWTask*> task_queue;

  for (int i = 0; i < 2; i++) {
    SWTask* task = new SWTask(opencl_env, chunk_size);
    task_queue.push_back(task);
  }

  // elements of the Regions record
  bseq1_t* seqs;
  mem_chain_v* chains;
  mem_alnreg_v* alnreg;

  int kernel_buffer_idx = 0;
  int task_num = 0;

  // For statistics
  uint64_t start_ts;
  uint64_t last_output_ts;

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

    uint64_t last_prepare_ts = getUs();
    uint64_t last_batch_ts = getUs();

    last_output_ts = last_batch_ts;
    start_idx = record.start_idx;
    batch_num = record.batch_num;
    seqs = record.seqs;
    chains = record.chains;
    alnreg = record.alnreg;

    int i = 0;
    bool reach_half = false;
    bool reach_end = false;
    while (i < batch_num) {
      if (task_num < chunk_size/2) {
        SWTask* task = task_queue.front();
        packReadData(aux, seqs+i, chains+i, alnreg+i, 
            task->i_data[0], 
            kernel_buffer_idx, task_num,
            task);
        i++;
      }
      else if (task_num >= chunk_size/2 && reach_half == false) {
        SWTask* task = task_queue.front();

        task->i_size[0] = kernel_buffer_idx/sizeof(int);
        task->o_size[0] = FPGA_RET_PARAM_NUM*task_num;

        kernel_buffer_idx = 0;
        reach_half = true;
      }
      else if (task_num < chunk_size) {
        SWTask* task = task_queue.front();

        packReadData(aux, seqs+i, chains+i, alnreg+i, 
            task->i_data[1], 
            kernel_buffer_idx, task_num,
            task);

        i++;
      }
      else if (task_num >= chunk_size) {
        SWTask* task = task_queue.front();
        task->i_size[1] = kernel_buffer_idx/sizeof(int);
        task->o_size[1] = FPGA_RET_PARAM_NUM*task_num - task->o_size[0];

        DLOG_IF(INFO, VLOG_IS_ON(3)) << "Prepare data for " << 
                task_num << " tasks takes " << 
                getUs() - last_batch_ts << " us";

        uint64_t output_ts = getUs();

        // start sending task to FPGA
        task->start(task_queue[1]);

        // circulate the task
        task_queue.push_back(task);
        task_queue.pop_front();

        task = task_queue.front();
        // if task is available
        if (task->i_size[0] > 0 && task->i_size[1] > 0) { 

          processOutput(task);
          DLOG_IF(INFO, VLOG_IS_ON(3)) << "This chunk of FPGA takes " << 
                                          getUs()- last_batch_ts << " us";
        }
        last_batch_ts = getUs();

        // reset 
        task_num = 0;
        kernel_buffer_idx = 0;
        reach_half = false;
      }
      if (i == batch_num - 1) {
        reach_end = true;
      }
      else {
        reach_end = false;
      }
    }
    DLOG_IF(INFO, VLOG_IS_ON(3)) << "Starting the remaining tasks";

    // finish the remain reads even with small task number 
    if (!reach_end && task_num != 0) {
      SWTask* task = task_queue.front();
      if (task_num < chunk_size/2 || reach_half == false) {
        task->i_size[0] = kernel_buffer_idx/sizeof(int);
        task->i_size[1] = 0;
        task->o_size[0] = FPGA_RET_PARAM_NUM*task_num;
        task->o_size[1] = 0;
      }
      else {
        task->i_size[1] = kernel_buffer_idx/sizeof(int);
        task->o_size[1] = FPGA_RET_PARAM_NUM*task_num - task->o_size[0];
      }

      task->start(task_queue[1]);
    }

    for (int iter = 0; iter < 2; iter++) {
      SWTask* task = task_queue.front();
      task_queue.push_back(task);
      task_queue.pop_front();

      task = task_queue.front();
      if (task->i_size[0] > 0) { 
        processOutput(task);
      }
    }
    DLOG_IF(INFO, VLOG_IS_ON(3)) << "Finished entire batch";

    // reset for one batch
    task_num = 0;
    kernel_buffer_idx = 0;
    reach_half = false;
    i = 0;

    RegionsRecord outputRecord;
    outputRecord.start_idx = start_idx;
    outputRecord.batch_num = batch_num;
    outputRecord.seqs = seqs;
    outputRecord.alnreg = alnreg;

    freeChains(chains, batch_num);
    pushOutput(outputRecord);

    DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished ChainsToRegions() on FPGA for "
      << getUs() - last_output_ts << " us";

    last_output_ts = getUs();
  }

  // delete everything
  while (!task_queue.empty()) {
    SWTask* task = task_queue.front();

    for (int k = 0; k < 2; k++) {
      clReleaseMemObject(task->i_buf[k]);
      clReleaseMemObject(task->o_buf[k]);
      free(task->i_data[k]);
      free(task->o_data[k]);
    }
    
    delete task->region_batch;
    delete task->chain_batch;
    delete task;

    task_queue.pop_front();
  }
}

ChainsRecord ChainsPipeFPGA::compute(ChainsRecord const & record) {
  uint64_t start_ts = getUs();
  ChainsRecord output = record;
  int batch_num = record.batch_num;
  mem_alnreg_v* alnreg = new mem_alnreg_v[batch_num];
  mem_chain_v* chains = record.chains;
  for (int i=0; i<batch_num; i++) {
    int chain_idx = 0;
    int seed_idx = chains[i].a[0].n -1;
    kv_init(alnreg[i]);
    for ( ; chain_idx < chains[i].n; 
         seed_idx > 0 ? seed_idx-- : seed_idx = chains[i].a[++chain_idx].n-1) {
      mem_seed_t* seed_array = &chains[i].a[chain_idx].seeds[seed_idx];
      // initialize the newreg
      mem_alnreg_t* newreg = kv_pushp(mem_alnreg_t, alnreg[i]);
      memset(newreg, 0, sizeof(mem_alnreg_t));
      newreg->score  = seed_array->len * aux->opt->a;
      newreg->truesc = seed_array->len * aux->opt->a;
      newreg->qb = seed_array->qbeg;
      newreg->rb = seed_array->rbeg;
      newreg->qe = seed_array->qbeg + seed_array->len;
      newreg->re = seed_array->rbeg + seed_array->len; 
      newreg->rid = chains[i].a[chain_idx].rid;
      newreg->seedlen0 = seed_array->len; 
      newreg->frac_rep = chains[i].a[chain_idx].frac_rep;
      newreg->w = aux->opt->w;
    }
  }
  output.alnreg = alnreg;
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished ChainsPipeFPGA() on FPGA for "
    << getUs() - start_ts << " us";
  return output;
}

ChainsRecord SeqsToChainsFPGA::compute(SeqsRecord const & seqs_record) {
  static int seqs_record_id = -1;
  seqs_record_id++;

  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started SeqsToChains() on FPGA";

  uint64_t start_ts = getUs();
  uint64_t ref_time = 0;

  bseq1_t* srseqs = seqs_record.seqs;
  uint64_t start_idx = seqs_record.start_idx;
  int batch_num = seqs_record.batch_num;

  mem_chain_v* chains = (mem_chain_v*)malloc(batch_num*sizeof(mem_chain_v));

#if 0
  const int max_task_size = 320;
  for (int i1 = 0; i1 < (batch_num+max_task_size-1)/max_task_size; i1++) {
    int task_size = std::min( max_task_size, batch_num-i1*max_task_size );
    for (int i2 = 0; i2 < task_size; i2++) {
      int i = i1*max_task_size+i2;
      // chains[i] = seq2chain(aux, &srseqs[i]);
#if 1
      bseq1_t *seqs = &srseqs[i];
      int j;
      mem_chain_v chn;
      for (j = 0; j < seqs->l_seq; ++j) // convert to 2-bit encoding if we have not done so
        seqs->seq[j] = seqs->seq[j] < 4? seqs->seq[j] : nst_nt4_table[(int)seqs->seq[j]];
      // chn = mem_chain(aux->opt, aux->idx->bwt, aux->idx->bns, seqs->l_seq, (uint8_t*)seqs->seq, 0);         // the 0 should be reconsidered
#if 1
      smem_aux_t *smem_aux = smem_aux_init();
      mem_collect_intv_new(aux->opt, aux->idx->bwt, seqs->l_seq, (uint8_t*)seqs->seq, smem_aux);
      chn = mem_chain_postprocess(aux->opt, aux->idx->bwt, aux->idx->bns, seqs->l_seq, (uint8_t*)seqs->seq, smem_aux);
      smem_aux_destroy(smem_aux);
#endif
      chn.n = mem_chain_flt(aux->opt, chn.n, chn.a);
      mem_flt_chained_seeds(aux->opt, aux->idx->bns, aux->idx->pac, seqs->l_seq, (uint8_t*)seqs->seq, chn.n, chn.a);

      chains[i] = chn;
#endif
    }
  }
#else
  // create SMemTasks
  std::deque<SMemTask*> task_queue;
  for (int i = 0; i < 2; i++) {
    task_queue.push_back( new SMemTask(opencl_env) );
  }

  const int max_task_size = 320;
  const int max_seq_len   = 150;
  SMemTask *task;
  for (int i1 = 0; i1 < (batch_num+max_task_size-1)/max_task_size; i1++) {
    task = task_queue.front();

    // prepare the task
    uint64_t inner_ts = getUs();
    int task_size = std::min( max_task_size, batch_num-i1*max_task_size );
    DLOG_IF(INFO, VLOG_IS_ON(3)) << "Task " << i1 << ": " << task_size;
    for (int i2 = 0; i2 < task_size; i2++) {
      int i = i1*max_task_size+i2;
      bseq1_t *seqs = &srseqs[i];
      if ( seqs->l_seq > max_seq_len ) {
        LOG(ERROR) << "Do not support sequences (len: " << seqs->l_seq
                   << ") with length bigger than " << max_seq_len << " on FPGA currently";
        throw std::runtime_error("Failed to process the oversized sequences");
      }
      for (int j = 0; j < seqs->l_seq; ++j) // convert to 2-bit encoding if we have not done so
        seqs->seq[j] = seqs->seq[j] < 4? seqs->seq[j] : nst_nt4_table[(int)seqs->seq[j]];
      if (seqs->l_seq < aux->opt->min_seed_len) {
        memset( &(task->i_seq_data[i2*max_seq_len]), UCHAR_MAX, max_seq_len*sizeof(char) );
      } else {
        memcpy( &(task->i_seq_data[i2*max_seq_len]), seqs->seq, seqs->l_seq*sizeof(char) );
        memset( &(task->i_seq_data[i2*max_seq_len+seqs->l_seq]), UCHAR_MAX, std::max(max_seq_len-(int)(seqs->l_seq), 0)*sizeof(char) );
      }
    }
    task->i_seq_num = task_size;
    task->i_seq_base_idx = i1*max_task_size;
    DLOG_IF(INFO, VLOG_IS_ON(3)) << "Preparation takes " << getUs() - inner_ts << " us";

    // launch the task
    task->start( task_queue[1] );

    // circulation
    task_queue.push_back( task );
    task_queue.pop_front();
    task = task_queue.front();
     
    // finish the previous task
    if ( task->i_seq_num > 0 ) {
      DLOG_IF(INFO, VLOG_IS_ON(3)) << "Switch from Task " << i1 << " to Task " << i1-1;
      task->finish();
        
      inner_ts = getUs();
      for (int i2 = 0; i2 < task->i_seq_num; i2++) {
        // unpack the data
        int i = task->i_seq_base_idx+i2;
        smem_aux_t *smem_aux = smem_aux_init();
        bwtintv_t *temp = smem_aux->mem.a;
        smem_aux->mem.a = &(task->o_mem_data[i2*task->max_intv_alloc_]);
        smem_aux->mem.n = task->o_num_data[i2];

        // postprocess
        mem_chain_v chn;
        bseq1_t *seqs = &srseqs[i];
        if (smem_aux->mem.n > task->max_intv_alloc_) {
          smem_aux->mem.a = temp;
          mem_collect_intv_new(aux->opt, aux->idx->bwt, seqs->l_seq, (uint8_t*)seqs->seq, smem_aux);
          DLOG_IF(WARNING, VLOG_IS_ON(3)) << "Seq " << i << " in SeqsRecord " << seqs_record_id << " overflowed. Redo on CPU";
        }
        chn = mem_chain_postprocess(aux->opt, aux->idx->bwt, aux->idx->bns, seqs->l_seq, (uint8_t*)seqs->seq, smem_aux);
        smem_aux->mem.a = temp;
        smem_aux_destroy(smem_aux);
        chn.n = mem_chain_flt(aux->opt, chn.n, chn.a);
        mem_flt_chained_seeds(aux->opt, aux->idx->bns, aux->idx->pac, seqs->l_seq, (uint8_t*)seqs->seq, chn.n, chn.a);

        chains[i] = chn;
      }
      task->i_seq_num = 0;
      DLOG_IF(INFO, VLOG_IS_ON(3)) << "Post-processing takes " << getUs() - inner_ts << " us";
    }
  } // loop for batches

  // finish the remaining task
  // circulation
  task = task_queue.front();
  task_queue.push_back( task );
  task_queue.pop_front();
  task = task_queue.front();
   
  if ( task->i_seq_num > 0 ) {
    DLOG_IF(INFO, VLOG_IS_ON(3)) << "Switch to Task " << (int)( task->i_seq_base_idx/max_task_size );
    task->finish();
      
    uint64_t inner_ts = getUs();
    for (int i2 = 0; i2 < task->i_seq_num; i2++) {
      // unpack the data
      int i = task->i_seq_base_idx+i2;
      smem_aux_t *smem_aux = smem_aux_init();
      bwtintv_t *temp = smem_aux->mem.a;
      smem_aux->mem.a = &(task->o_mem_data[i2*task->max_intv_alloc_]);
      smem_aux->mem.n = task->o_num_data[i2];

      // postprocess
      mem_chain_v chn;
      bseq1_t *seqs = &srseqs[i];
      if (smem_aux->mem.n > task->max_intv_alloc_) {
        smem_aux->mem.a = temp;
        mem_collect_intv_new(aux->opt, aux->idx->bwt, seqs->l_seq, (uint8_t*)seqs->seq, smem_aux);
        DLOG_IF(WARNING, VLOG_IS_ON(3)) << "Seq " << i << " in SeqsRecord " << seqs_record_id << " overflowed. Redo on CPU";
      }
      chn = mem_chain_postprocess(aux->opt, aux->idx->bwt, aux->idx->bns, seqs->l_seq, (uint8_t*)seqs->seq, smem_aux);
      smem_aux->mem.a = temp;
      smem_aux_destroy(smem_aux);
      chn.n = mem_chain_flt(aux->opt, chn.n, chn.a);
      mem_flt_chained_seeds(aux->opt, aux->idx->bns, aux->idx->pac, seqs->l_seq, (uint8_t*)seqs->seq, chn.n, chn.a);

      chains[i] = chn;
    }
    task->i_seq_num = 0;
    DLOG_IF(INFO, VLOG_IS_ON(3)) << "Post-processing takes " << getUs() - inner_ts << " us";
  }

  // delete tasks
  while (!task_queue.empty()) {
    delete task_queue.front();
    task_queue.pop_front();
  }
#endif

  ChainsRecord ret;
  ret.start_idx    = seqs_record.start_idx;
  ret.batch_num    = batch_num;
  ret.seqs         = srseqs;
  ret.chains       = chains;
 
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished SeqsToChains() in " 
    << getUs() - start_ts << " us";

  return ret;
}

smem_aux_t *smem_aux_init()
{
	smem_aux_t *a;
	a = (smem_aux_t *)calloc(1, sizeof(smem_aux_t));
	a->tmpv[0] = (bwtintv_v *)calloc(1, sizeof(bwtintv_v));
	a->tmpv[1] = (bwtintv_v *)calloc(1, sizeof(bwtintv_v));
	return a;
}

void smem_aux_destroy(smem_aux_t *a)
{
	free(a->tmpv[0]->a); free(a->tmpv[0]);
	free(a->tmpv[1]->a); free(a->tmpv[1]);
	free(a->mem.a); free(a->mem1.a);
	free(a);
}
