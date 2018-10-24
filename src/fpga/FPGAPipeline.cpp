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
#include <errno.h>
#include <vector>
#include <iostream>
#include <fstream>

#include "bwa_wrapper.h"
#include "config.h"
#include "util.h"
#include "allocation_wrapper.h"
#include "Pipeline.h"
#include "FPGAPipeline.h"
#include "SWTask.h"
#include "SMemTask.h"


inline void processOutput(SWTask* task, uint64_t &post_ts) {
  // finish task and get output buffers
  int total_task_num = task->o_size/FPGA_RET_PARAM_NUM;
  int actual_tasks = 0;

  uint64_t post_start_ts = getUs();
  int wrong_half = -1;
  int num_redo = 0;

  // check fpga results correctness
  while (1) {
    bool wrong_results = false;
    int task_num = task->o_size / FPGA_RET_PARAM_NUM;
    short* kernel_output = task->o_data;
    for (int i = 0; i < task_num; i++) {
      int seed_idx = ((int)(kernel_output[1+FPGA_RET_PARAM_NUM*2*i])<<16) |
        kernel_output[0+FPGA_RET_PARAM_NUM*2*i];
      if (seed_idx > actual_tasks)
        actual_tasks = seed_idx;
      if (seed_idx > total_task_num || seed_idx < 0) {
        DLOG(ERROR) << "Wrong fpga results at"
                    << " task_num = " << i << " seed_idx = " << seed_idx;
        break;
      }
    }
    if ((!wrong_results) && (total_task_num != actual_tasks + 1)) {
      DLOG(ERROR) << "Wrong fpga results due to problem with seedindex";
      wrong_results = true;
    }
    if (wrong_results) {
      if (num_redo < 10) {
        DLOG(WARNING) << "Incurred wrong fpga results. Redo.";
        task->redo();
        num_redo++;
      }
      else {
#if 0
        // dump results   
        FILE* fout = fopen("dump-input.dat", "wb+");
        fwrite(task->i_data, sizeof(int), task->i_size, fout);
        fclose(fout);
        DLOG(INFO) << "dump input data of output to dump-input.dat";
#endif
        throw fpgaResultsError("wrong results in smithwaterman kernel");
      }
    }
    else {
      if (num_redo != 0) {
        DLOG(WARNING) << "Redo " << num_redo << " times.";
      }
      break;
    }
  }


  // apply fpga results
  mem_alnreg_t** region_batch = task->region_batch;
  mem_chain_t**  chain_batch = task->chain_batch;
  int task_num = task->o_size / FPGA_RET_PARAM_NUM;
  short* kernel_output = task->o_data;

  for (int i = 0; i < task_num; i++) {
    int seed_idx = ((int)(kernel_output[1+FPGA_RET_PARAM_NUM*2*i])<<16) |
      kernel_output[0+FPGA_RET_PARAM_NUM*2*i];

    if (seed_idx > actual_tasks) {
      actual_tasks = seed_idx;
    }
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
    int seedcov = 0;
    for (int chain_id=0; chain_id < chain->n; chain_id++) {
      const mem_seed_t *seed = &chain->seeds[chain_id];
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
  task->i_size = 0;
  task->o_size = 0;

  DLOG_IF(INFO, VLOG_IS_ON(4)) << "Total " << actual_tasks << " tasks in output";
  post_ts += getUs() - post_start_ts;
  DLOG_IF(INFO, VLOG_IS_ON(4)) << "Process output takes " 
    << getUs() - post_start_ts << " us";
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
  // empty filter & overflow filter
#if 1 
  int total_seed_num = 0;
  for (int chain_idx = 0; chain_idx < chains->n; chain_idx++)
    total_seed_num += chains->a[chain_idx].n;
  if (chains->n == 0 || total_seed_num == 0 || chains->n >= 2000) {
    kv_init(*alnregs);
    for (int j = 0; j < chains->n; j++) {
      mem_chain2aln(aux->opt,
                    aux->idx->bns,
                    aux->idx->pac,
                    seq->l_seq,
          (uint8_t*)seq->seq,
                    &chains->a[j],
                    alnregs);
    }
    return;
  }
#endif

  // preprocess
  int chain_idx = 0;
  int seed_idx = chains->a[0].n -1;
  kv_init(*alnregs);
  for ( ; chain_idx < chains->n; 
        seed_idx > 0 ? seed_idx-- : seed_idx = chains->a[++chain_idx].n-1) {
    mem_seed_t* seed_array = &chains->a[chain_idx].seeds[seed_idx];
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

  mem_alnreg_t** region_batch = task->region_batch;
  mem_chain_t** chain_batch = task->chain_batch;
  int chunk_size = FLAGS_chunk_size;
  mem_chainref_t* ref;
  seed_idx = chains->a[0].n - 1;
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
          task->o_data = (short*)realloc(task->o_data, task->max_o_size_ * sizeof(short));
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
    task->i_data = (char*)realloc(task->i_data, task->max_i_size_);
  }

  free(ref);
}

inline void finishUpOnCPU(RegionsRecord record, int start_seq) {
  int batch_num = record.batch_num;

  bseq1_t*      seqs   = record.seqs;
  mem_chain_v*  chains = record.chains;
  mem_alnreg_v* alnreg = record.alnreg;
  for (int i = start_seq; i < batch_num; i++) {
    kv_init(alnreg[i]);
    for (int j = 0; j < chains[i].n; j++) {
      mem_chain2aln(
          aux->opt, 
          aux->idx->bns, 
          aux->idx->pac,
          seqs[i].l_seq,
          (uint8_t*)seqs[i].seq,
          &chains[i].a[j],
          alnreg+i);
    }
  }
}


void ChainsToRegionsFPGA::compute(int wid) {
  DLOG(INFO) << "start FPGA worker #" << wid;

  int chunk_size = FLAGS_chunk_size;
  uint64_t g_prep_ts=0, g_write_ts=0, g_enq_ts=0, g_deq_ts=0, g_read_ts=0, g_post_ts=0;

  // Create SWTasks
  std::deque<SWTask*> task_queue;
  for (int i = 0; i < 2; i++) {
    SWTask* task = new SWTask(opencl_env, chunk_size);
    if (NULL == task) {
      std::string err_string = "Memory allocation failed";
      if (errno==12)
        err_string += " due to out-of-memory";
      else
        err_string += " due to internal failure"; 
      throw std::runtime_error(err_string);
    }
    task_queue.push_back(task);
  }

  bool flag_need_reads = false;
  bool flag_more_reads = true;

  ChainsRecord record; 
  while (flag_more_reads) { 
    // get one Chains record
    bool ready = this->getInput(record);
    uint64_t wait_input_start_ts = getUs();
    while (!this->isFinal() && !ready) {
      boost::this_thread::sleep_for(boost::chrono::microseconds(10));
      ready = this->getInput(record);
    }
    DLOG_IF(INFO, VLOG_IS_ON(2)) << "Wait for input for FPGA takes " << getUs() - wait_input_start_ts << " us";
    if(!ready) {
      n_active_--;
      if (n_active_ == 0 && cpu_stage_ != NULL) cpu_stage_->setUseAccx(false);
      flag_more_reads = false;
      break;
    }

    uint64_t prep_ts=0, enq_ts=0, deq_ts=0, post_ts=0;
    DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started ChainsToRegions() on FPGA ";
    DLOG_IF(INFO, VLOG_IS_ON(3)) << "SW start idx: " << record.start_idx;
    DLOG_IF(INFO, VLOG_IS_ON(3)) << "Record tag: " << record.tag;

    uint64_t start_ts = getUs();
    uint64_t pre_start_ts = getUs();

    int        start_idx = record.start_idx;
    int        batch_num = record.batch_num;
    bseq1_t*        seqs = record.seqs;
    mem_chain_v*  chains = record.chains;
    mem_alnreg_v* alnreg = new mem_alnreg_v[batch_num];
    if (NULL == alnreg) {
      std::string err_string = "Memory allocation failed";
      if (errno==12)
        err_string += " due to out-of-memory";
      else
        err_string += " due to internal failure"; 
      throw std::runtime_error(err_string);
    }

    RegionsRecord outputRecord;
    outputRecord.start_idx = start_idx;
    outputRecord.batch_num = batch_num;
    outputRecord.seqs = seqs;
    outputRecord.chains = chains;
    outputRecord.alnreg = alnreg;

    try {
      int chunk_id = 0;
      int task_num = 0;
      int kernel_buffer_idx = 0;
      task_queue.front()->start_seq = 0;

      for (int i = 0; i < batch_num; i++) {
        SWTask* task = task_queue.front();
        packReadData(aux, seqs+i, chains+i, alnreg+i,
                     task->i_data,
                     kernel_buffer_idx, task_num,
                     task);

        if (task_num >= chunk_size/2) {
          task->i_size = kernel_buffer_idx/sizeof(int);
          task->o_size = FPGA_RET_PARAM_NUM*task_num;
          if (task->o_size == 0) task->i_size = 0;
          task->end_seq = i;
          prep_ts += getUs() - pre_start_ts;
          DLOG_IF(INFO, VLOG_IS_ON(4)) << "Prepare data for " << task_num << " tasks takes "
                                       << getUs() - pre_start_ts << " us";

          // start sending task to FPGA
          DLOG_IF(INFO, VLOG_IS_ON(3)) << "Start chunk " << chunk_id;
          DLOG_IF(INFO, VLOG_IS_ON(4)) << "Seqs in chunk: [" << task->start_seq << ", " << task->end_seq << ").";
          uint64_t enq_start_ts = getUs();
          task->start(task_queue[1]);
          enq_ts += getUs() - enq_start_ts;

          // circulate the task
          task_queue.push_back(task);
          task_queue.pop_front();
          task = task_queue.front();

          // if task is available
          if (task->o_size > 0) { 
            DLOG_IF(INFO, VLOG_IS_ON(3)) << "Wait chunk " << chunk_id-1;
            uint64_t deq_start_ts = getUs();
            task->finish();
            deq_ts += getUs() - deq_start_ts;
            processOutput(task, post_ts);
          }
          chunk_id++;

          // reset 
          task_num = 0;
          kernel_buffer_idx = 0;
          task->start_seq = i;
          task->end_seq = 0;

          pre_start_ts = getUs();
        }
      }


      DLOG_IF(INFO, VLOG_IS_ON(4)) << "Starting the remaining tasks";

      // finish the remain reads even with small task number 
      if (task_num != 0) {
        SWTask* task = task_queue.front();
        task->i_size = kernel_buffer_idx/sizeof(int);
        task->o_size = FPGA_RET_PARAM_NUM*task_num;
        if (task->o_size == 0) task->i_size = 0;
        task->end_seq = batch_num;
        prep_ts += getUs() - pre_start_ts;

        DLOG_IF(INFO, VLOG_IS_ON(3)) << "Start chunk " << chunk_id;
        DLOG_IF(INFO, VLOG_IS_ON(4)) << "Seqs in chunk: [" << task->start_seq << ", " << task->end_seq << ").";
        uint64_t enq_start_ts = getUs();
        task->start(task_queue[1]);
        enq_ts += getUs() - enq_start_ts;
      }

      for (int iter = 0; iter < 2; iter++) {
        SWTask* task = task_queue.front();
        task_queue.push_back(task);
        task_queue.pop_front();

        task = task_queue.front();
        if (task->o_size > 0) { 
          DLOG_IF(INFO, VLOG_IS_ON(3)) << "Wait chunk " << chunk_id-1;
          uint64_t deq_start_ts = getUs();
          task->finish();
          deq_ts += getUs() - deq_start_ts;
          processOutput(task, post_ts);
        }
        chunk_id++;
      }
    }
    catch (fpgaHangError &e) {
      LOG_IF(WARNING, VLOG_IS_ON(1)) << "FPGA hangs, restart execution on CPU";
      DVLOG(1) << "batch start idx: " << start_idx << "\n"
               << "batch num: " << batch_num << "\n"
               << "tag: " << record.tag << "\n"
               << "task start seq: " << task_queue.front()->start_seq << "\n"
               << "task end seq: " << task_queue.front()->end_seq << "\n"
               << task_queue.front()->i_size << "\n"
               << task_queue.front()->o_size << "\n";

      n_active_--;
      if (n_active_ == 0 && cpu_stage_ != NULL) cpu_stage_->setUseAccx(false);
      SWTask *err_task = task_queue.front();
      finishUpOnCPU(outputRecord, err_task->start_seq);
      pushOutput(outputRecord);
      break;
    }
    catch (fpgaResultsError &e) {
      DLOG(WARNING) << "FPGA generated wrong results in smithwaterman kernel.";
      n_active_--;
      if (n_active_ == 0 && cpu_stage_ != NULL) cpu_stage_->setUseAccx(false);
      SWTask *err_task = task_queue.front();
      finishUpOnCPU(outputRecord, err_task->start_seq);
      pushOutput(outputRecord);
      break;
    }
   
    DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished ChainsToRegions() on FPGA for "
                                 << getUs() - start_ts << " us";
    //freeChains(chains, batch_num);

    DLOG_IF(INFO, VLOG_IS_ON(2)) << "[sw] Pre-processing takes " << prep_ts << " us";
    DLOG_IF(INFO, VLOG_IS_ON(2)) << "[sw] Enqueuing task takes " << enq_ts << " us";
    DLOG_IF(INFO, VLOG_IS_ON(2)) << "[sw] Dequeuing task takes " << deq_ts << " us";
    DLOG_IF(INFO, VLOG_IS_ON(2)) << "[sw] Post-processing takes " << post_ts << " us";
    g_prep_ts  += prep_ts;
    g_enq_ts   += enq_ts;
    g_deq_ts   += deq_ts;
    g_post_ts  += post_ts;

    pushOutput(outputRecord);
  }

  // delete everything
  while (!task_queue.empty()) {
    SWTask* task = task_queue.front();
    delete task;
    task_queue.pop_front();
  }
  DLOG_IF(INFO, VLOG_IS_ON(2)) << "[sw_global] Pre-processing takes " << g_prep_ts << " us";
  DLOG_IF(INFO, VLOG_IS_ON(2)) << "[sw_global] Enqueuing task takes " << g_enq_ts << " us";
  DLOG_IF(INFO, VLOG_IS_ON(2)) << "[sw_global] Dequeuing task takes " << g_deq_ts << " us";
  DLOG_IF(INFO, VLOG_IS_ON(2)) << "[sw_global] Post-processing takes " << g_post_ts << " us";
}


void SeqsToChainsFPGA::compute(int wid) {
  DLOG(INFO) << "start FPGA worker #" << wid;
  uint64_t g_prep_ts=0, g_write_ts=0, g_enq_ts=0, g_deq_ts=0, g_read_ts=0, g_post_ts=0;

  const int max_task_size = SMemTask::max_i_seq_num_;
  const int max_seq_len   = SMemTask::max_i_seq_len_;
  const int max_intv_alloc = SMemTask::max_intv_alloc_;

  // create SMemTasks
  std::deque<SMemTask*> task_queue;
  for (int i = 0; i < 2; i++) {
    task_queue.push_back( new SMemTask(opencl_env) );
  }

  bool flag_more_reads = true;

  SeqsRecord seqs_record;
  while (flag_more_reads) {
    // get one Seqs record
    bool ready = this->getInput(seqs_record);
    while (!this->isFinal() && ! ready) {
      boost::this_thread::sleep_for(boost::chrono::microseconds(10));
      ready = this->getInput(seqs_record);
    }
    if (!ready) {
      n_active_--;
      if (n_active_ == 0 && cpu_stage_ != NULL) cpu_stage_->setUseAccx(false);
      flag_more_reads = false;
      break;
    }

    uint64_t prep_ts=0, write_ts=0, enq_ts=0, deq_ts=0, read_ts=0, post_ts=0;
    uint64_t start_ts = getUs();
    DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started SeqsToChains() on FPGA";

    bseq1_t* srseqs = seqs_record.seqs;
    uint64_t start_idx = seqs_record.start_idx;
    int batch_num = seqs_record.batch_num;

    bwtintv_t** bwtintvs = (bwtintv_t**)malloc(batch_num*sizeof(bwtintv_t*)); 
    size_t* bwtintv_nums = (size_t*)malloc(batch_num*sizeof(size_t));

#if 0
    for (int i1 = 0; i1 < (batch_num+max_task_size-1)/max_task_size; i1++) {
      int task_size = std::min( max_task_size, batch_num-i1*max_task_size );
      for (int i2 = 0; i2 < task_size; i2++) {
        int i = i1*max_task_size+i2;
#if 0
        chains[i] = seq2chain(aux, &srseqs[i]);
#else
        bseq1_t *seqs = &srseqs[i];
        int j;
        mem_chain_v chn;
        for (j = 0; j < seqs->l_seq; ++j) // convert to 2-bit encoding if we have not done so
          seqs->seq[j] = seqs->seq[j] < 4? seqs->seq[j] : nst_nt4_table[(int)seqs->seq[j]];
#if 0
        // chn = mem_chain(aux->opt, aux->idx->bwt, aux->idx->bns, seqs->l_seq, (uint8_t*)seqs->seq, 0);         // the 0 should be reconsidered
#else
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
#endif
    SMemTask *task;
    for (int i1 = 0; i1 < (batch_num+max_task_size-1)/max_task_size; i1++) {
      task = task_queue.front();

      // prepare the task
      uint64_t inner_ts = getUs();
      int task_size = std::min( max_task_size, batch_num-i1*max_task_size );
      DLOG_IF(INFO, VLOG_IS_ON(4)) << "Task " << i1 << ": " << task_size;

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
          memset( &(task->i_seq_data[i2*max_seq_len]), UCHAR_MAX, max_seq_len*sizeof(uint8_t) );
          task->i_seq_len_data[i2] = 0;
        } else {
          memcpy( &(task->i_seq_data[i2*max_seq_len]), seqs->seq, seqs->l_seq*sizeof(uint8_t) );
          task->i_seq_len_data[i2] = (uint8_t)seqs->l_seq; 
        }
      }
      task->i_seq_num = task_size;
      task->i_seq_base_idx = i1*max_task_size;
      
      prep_ts += getUs() - inner_ts;
      DLOG_IF(INFO, VLOG_IS_ON(4)) << "Preparation takes " << getUs() - inner_ts << " us";

      // launch the task
      task->start( task_queue[1], write_ts, enq_ts );

      // circulation
      task_queue.push_back( task );
      task_queue.pop_front();
      task = task_queue.front();
       
      // finish the previous task
      if ( task->i_seq_num > 0 ) {
        DLOG_IF(INFO, VLOG_IS_ON(4)) << "Switch from Task " << i1 << " to Task " << i1-1;
        task->finish( deq_ts, read_ts );
          
        inner_ts = getUs();
        for (int i2 = 0; i2 < task->i_seq_num; i2++) {
          // unpack the data
          int i = task->i_seq_base_idx+i2;
          size_t n = task->o_num_data[i2];
          // postprocess
          bseq1_t *seqs = &srseqs[i];
          if (n > max_intv_alloc) {
            DLOG_IF(WARNING, VLOG_IS_ON(3)) << "Seq " << i << " in SeqsRecord[start_idx:" << seqs_record.start_idx << "] overflowed. Redo on CPU";
            smem_aux_t *smem_aux = smem_aux_init();
            mem_collect_intv_new(aux->opt, aux->idx->bwt, seqs->l_seq, (uint8_t*)seqs->seq, smem_aux);
            n = smem_aux->mem.n;
            bwtintvs[i] = (bwtintv_t*)malloc(n*sizeof(bwtintv_t));
            memcpy(bwtintvs[i], smem_aux->mem.a, n*sizeof(bwtintv_t));
            smem_aux_destroy(smem_aux);
          }
          else {
            bwtintvs[i] = (bwtintv_t*)malloc(n*sizeof(bwtintv_t));
            memcpy(bwtintvs[i], &(task->o_mem_data[i2*task->max_intv_alloc_]), n*sizeof(bwtintv_t));
          }
          bwtintv_nums[i] = n;
        }

        task->i_seq_num = 0;
        post_ts += getUs() - inner_ts;
        DLOG_IF(INFO, VLOG_IS_ON(4)) << "Post-processing takes " << getUs() - inner_ts << " us";
      }
    } // loop for batches

    // finish the remaining task
    // circulation
    task = task_queue.front();
    task_queue.push_back( task );
    task_queue.pop_front();
    task = task_queue.front();
     
    if ( task->i_seq_num > 0 ) {
      DLOG_IF(INFO, VLOG_IS_ON(4)) << "Switch to Task " << (int)( task->i_seq_base_idx/max_task_size );
      task->finish( deq_ts, read_ts );
        
      uint64_t inner_ts = getUs();
      for (int i2 = 0; i2 < task->i_seq_num; i2++) {
        // unpack the data
        int i = task->i_seq_base_idx+i2;
        size_t n = task->o_num_data[i2];
        bseq1_t *seqs = &srseqs[i];
        if (n > max_intv_alloc) {
          DLOG_IF(WARNING, VLOG_IS_ON(3)) << "Seq " << i << " in SeqsRecord[start_idx:" << seqs_record.start_idx << "] overflowed. Redo on CPU";
          smem_aux_t *smem_aux = smem_aux_init();
          mem_collect_intv_new(aux->opt, aux->idx->bwt, seqs->l_seq, (uint8_t*)seqs->seq, smem_aux);
          n = smem_aux->mem.n;
          bwtintvs[i] = (bwtintv_t*)malloc(n*sizeof(bwtintv_t));
          memcpy(bwtintvs[i], smem_aux->mem.a, n*sizeof(bwtintv_t));
          smem_aux_destroy(smem_aux);
        }
        else {
          bwtintvs[i] = (bwtintv_t*)malloc(n*sizeof(bwtintv_t));
          memcpy(bwtintvs[i], &(task->o_mem_data[i2*task->max_intv_alloc_]), n*sizeof(bwtintv_t));
        }
        bwtintv_nums[i] = n;
      }

      task->i_seq_num = 0;
      post_ts += getUs() - inner_ts;
    }

    ChainsRecord output;
    output.start_idx    = seqs_record.start_idx;
    output.batch_num    = batch_num;
    output.seqs         = srseqs;
    output.bwtintvs     = bwtintvs;
    output.bwtintv_nums = bwtintv_nums;
    output.tag          = 0;
 
    DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished SeqsToChains() on FPGA in " 
                                 << getUs() - start_ts << " us";
    DLOG_IF(INFO, VLOG_IS_ON(2)) << "[smem] Pre-processing takes " << prep_ts << " us";
    DLOG_IF(INFO, VLOG_IS_ON(2)) << "[smem] Writing buffer takes " << write_ts << " us";
    DLOG_IF(INFO, VLOG_IS_ON(2)) << "[smem] Enqueuing task takes " << enq_ts << " us";
    DLOG_IF(INFO, VLOG_IS_ON(2)) << "[smem] Dequeuing task takes " << deq_ts << " us";
    DLOG_IF(INFO, VLOG_IS_ON(2)) << "[smem] Reading buffer takes " << read_ts << " us";
    DLOG_IF(INFO, VLOG_IS_ON(2)) << "[smem] Post-processing takes " << post_ts << " us";
    g_prep_ts  += prep_ts;
    g_write_ts += write_ts;
    g_enq_ts   += enq_ts;
    g_deq_ts   += deq_ts;
    g_read_ts  += read_ts;
    g_post_ts  += post_ts;
    
    pushOutput(output);
  }

  // delete tasks
  while (!task_queue.empty()) {
    delete task_queue.front();
    task_queue.pop_front();
  }
  DLOG_IF(INFO, VLOG_IS_ON(2)) << "[smem_global] Pre-processing takes " << g_prep_ts << " us";
  DLOG_IF(INFO, VLOG_IS_ON(2)) << "[smem_global] Writing buffer takes " << g_write_ts << " us";
  DLOG_IF(INFO, VLOG_IS_ON(2)) << "[smem_global] Enqueuing task takes " << g_enq_ts << " us";
  DLOG_IF(INFO, VLOG_IS_ON(2)) << "[smem_global] Dequeuing task takes " << g_deq_ts << " us";
  DLOG_IF(INFO, VLOG_IS_ON(2)) << "[smem_global] Reading buffer takes " << g_read_ts << " us";
  DLOG_IF(INFO, VLOG_IS_ON(2)) << "[smem_global] Post-processing takes " << g_post_ts << " us";
}
