#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <vector>
#include <list>
#include "bwa/bwamem.h"
#include "bwa/bntseq.h"
#include "bwa/ksw.h"
#include "bwa/kvec.h"
#include "bwa/ksort.h"
#include "bwa/utils.h"
#include "bwa_wrapper.h"
#include "blaze/AccAgent.h"

#include "SWTask.h"
#include "SWRead.h"

#define FPGA_RET_PARAM_NUM 5
#define chunk_size 2000

#define USE_FPGA

// hw data structures
//extern "C"{
//void sw_top (int *a, int *output_a, int __inc);
//}

inline int Int2CharArray(char* arr, int idx, int num)
{
  arr[idx] = (char)(num&0xff);
  arr[idx+1] = (char)((num>>8)&0xff);
  arr[idx+2] = (char)((num>>16)&0xff);
  arr[idx+3] = (char)((num>>24)&0xff);
  return idx+4 ;
}

inline int Short2CharArray(char *arr, int idx, short num)
{
  arr[idx] = (char)(num & 0xff);
  arr[idx+1] = (char)((num>>8) & 0xff);
  return idx+2 ;
}

void SwFPGA(
    std::vector<ExtParam*> &tasks, 
    int batch_num, 
    mem_alnreg_v *av)
{
  uint64_t start_ts = blaze::getUs();

  int Buf1Len = 32 + 32*batch_num;
  char* buf1 = new char[Buf1Len];
  //------------------ store the public options at the beginning-----------------
  buf1[0] = (char)tasks[0]->oDel;
  buf1[1] = (char)tasks[0]->eDel;
  buf1[2] = (char)tasks[0]->oIns;
  buf1[3] = (char)tasks[0]->eIns;
  buf1[4] = (char)tasks[0]->penClip5;
  buf1[5] = (char)tasks[0]->penClip3;
  buf1[6] = (char)tasks[0]->w;
  Int2CharArray(buf1,8,batch_num);

  //-------------------pack the batch of parameters of each SW--------------------
  int i = 0;
  int LeftMaxIns = 0;
  int LeftMaxDel = 0;
  int RightMaxIns = 0;
  int RightMaxDel = 0;
  int TaskPos = 0 ;
  TaskPos = Buf1Len >> 2;
  int buf1idx = 32;
  while(i < batch_num)
  {
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(tasks[i]->leftQlen));
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(tasks[i]->leftRlen));
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(tasks[i]->rightQlen));
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(tasks[i]->rightRlen));
    buf1idx = Int2CharArray(buf1,buf1idx,TaskPos);
    TaskPos += ((((tasks[i]->leftQlen + tasks[i]->leftRlen + tasks[i]->rightQlen + tasks[i]->rightRlen)+1)/2)+3)/4;
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(tasks[i]->regScore));
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(tasks[i]->qBeg));
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(tasks[i]->h0));
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(tasks[i]->idx));
    LeftMaxIns = (short)((double)(tasks[i]->leftQlen*1+ tasks[i]->penClip5 -tasks[i]->oIns)/tasks[i]->eIns+1);
    LeftMaxDel = (short)((double)(tasks[i]->leftQlen*1 + tasks[i]->penClip5 -tasks[i]->oDel)/tasks[i]->eDel+1);
    RightMaxIns = (short)((double)(tasks[i]->rightQlen*1 + tasks[i]->penClip3 -tasks[i]->oIns)/tasks[i]->eIns+1);
    RightMaxDel = (short)((double)(tasks[i]->rightQlen*1 + tasks[i]->penClip3 -tasks[i]->oDel)/tasks[i]->eDel+1);         // 1 stands for tasks[i]->mat.max
    buf1idx = Short2CharArray(buf1,buf1idx,LeftMaxIns);
    buf1idx = Short2CharArray(buf1,buf1idx,LeftMaxDel);
    buf1idx = Short2CharArray(buf1,buf1idx,RightMaxIns);
    buf1idx = Short2CharArray(buf1,buf1idx,RightMaxDel);    
    buf1idx = Int2CharArray(buf1, buf1idx, i);
    i = i+1;
  }

  fprintf(stderr, "FPGA preparation used %dus until buf1\n", blaze::getUs()-start_ts); 
  char *buf2 = new char[(TaskPos<<2)-Buf1Len];
  int input_length = TaskPos<<2;
  int buf2idx = 0;
  i = 0;
  int j = 0;
  int TmpIntVar = 0;
  char Tmpchar = 0;  
  int Counter8 = 0;
 /* 
  for (int i = 0; i < batch_num ; i++ )
  {
    if(tasks[i]->leftQlen > 0){
      for (int j = 0;j < tasks[i]->leftQlen;j++ ){
        Counter8 = Counter8 + 1;
        buf2[buf2idx] = buf2[buf2idx]<<4|(tasks[i]->leftQs[j] & 0x0f);
        if(Counter8 %2 ==0){
          buf2idx ++; 
        }
      }
    }
    if(tasks[i]->rightQlen > 0){
      for (int j = 0;j < tasks[i]->rightQlen;j++ ){
        Counter8 = Counter8 + 1;
        buf2[buf2idx] = buf2[buf2idx]<<4 |(tasks[i]->rightQs[j] & 0x0f);
        if(Counter8 %2 ==0){
          buf2idx ++; 
        }
      }
    }
    if(tasks[i]->leftRlen > 0){
      for (int j = 0;j < tasks[i]->leftRlen;j++ ){
        Counter8 = Counter8 + 1;
        buf2[buf2idx] = buf2[buf2idx]<<4 |(tasks[i]->leftRs[j] & 0x0f);
        if(Counter8 %2 ==0){
          buf2idx ++; 
        }
      }
    }
    if(tasks[i]->rightRlen > 0){
      for (int j = 0;j < tasks[i]->rightRlen;j++ ){
        Counter8 = Counter8 + 1;
        buf2[buf2idx] = buf2[buf2idx]<<4 |(tasks[i]->rightRs[j] & 0x0f);
        if(Counter8 %2 ==0){
          buf2idx ++; 
        }
      }
    }
  }
 */ 
  while(i < batch_num)
  {
    if(tasks[i]->leftQlen > 0)
    {
      j = 0;
      while(j < tasks[i]->leftQlen)
      {
        Counter8 = Counter8 + 1;
        TmpIntVar = TmpIntVar <<4 | ((int)tasks[i]->leftQs[j] & 0x0f);
        if(Counter8 % 8 ==0)
          buf2idx = Int2CharArray(buf2,buf2idx,TmpIntVar);
        j = j + 1;
      }
    }
    if(tasks[i]->rightQlen > 0)
    {
      j = 0;
      while(j < tasks[i]->rightQlen)
      {
        Counter8 = Counter8 + 1;
        TmpIntVar = TmpIntVar <<4 | ((int)tasks[i]->rightQs[j] & 0x0f);
        if(Counter8 % 8 ==0)
          buf2idx = Int2CharArray(buf2,buf2idx,TmpIntVar);
        j = j + 1;
      }
    }
    if(tasks[i]->leftRlen > 0)
    {
      j = 0;
      while(j < tasks[i]->leftRlen)
      {
        Counter8 = Counter8 + 1;
        TmpIntVar = TmpIntVar <<4 | ((int)tasks[i]->leftRs[j] & 0x0f);
        if(Counter8 % 8 ==0)
          buf2idx = Int2CharArray(buf2,buf2idx,TmpIntVar);
        j = j + 1;
      }
    }
    if(tasks[i]->rightRlen > 0)
    {
      j = 0;
      while(j < tasks[i]->rightRlen)
      {
        Counter8 = Counter8 + 1;
        TmpIntVar = TmpIntVar <<4 | ((int)tasks[i]->rightRs[j] & 0x0f);
        if(Counter8 % 8 ==0)
          buf2idx = Int2CharArray(buf2,buf2idx,TmpIntVar);
        j = j + 1;
      }
    }
    if(Counter8 %8 != 0)
    {
      while(Counter8 %8 != 0 )
      {
        TmpIntVar = TmpIntVar << 4;
        Counter8 = Counter8 + 1;
      }
      buf2idx = Int2CharArray(buf2,buf2idx,TmpIntVar);
    }
    i = i + 1;
  }
  
  fprintf(stderr, "FPGA preparation used %dus until buf2\n", blaze::getUs()-start_ts); 
  int *data_ptr = new int[input_length/4];
  short* output_ptr = new short[FPGA_RET_PARAM_NUM*batch_num*2];

  memcpy(data_ptr,buf1,Buf1Len*sizeof(char) );
  memcpy(&data_ptr[Buf1Len/4],buf2,input_length-Buf1Len);                          
  delete buf1;
  delete buf2;             
  fprintf(stderr, "FPGA preparation used %dus\n", blaze::getUs()-start_ts); 

  // sw_top (data_ptr, (int *)output_ptr,batch_num);
  
  start_ts = blaze::getUs();
  blaze::Task_ptr fpga_task = agent->createTask(acc_id);
  if (!fpga_task) {
    throw blaze::internalError("Task is not created");
  }
  agent->writeInput(fpga_task, acc_id, data_ptr, 1, input_length/4, sizeof(int));
  agent->writeInput(fpga_task, acc_id, &batch_num, 1, 1, sizeof(int));
  agent->readOutput(fpga_task, output_ptr, FPGA_RET_PARAM_NUM*batch_num*4);
  fprintf(stderr, "FPGA kernel used %dus\n", blaze::getUs()-start_ts); 

  start_ts = blaze::getUs();
  for (int i = 0; i < batch_num; i++) {  
    
    mem_alnreg_t *newreg =tasks[i]->newreg;
    int task_idx = ((int)(output_ptr[1+FPGA_RET_PARAM_NUM*2*i])<<16) |
                    output_ptr[0+FPGA_RET_PARAM_NUM*2*i];
    int regs_idx = tasks[task_idx]->idx;

    newreg->qb = output_ptr[2+FPGA_RET_PARAM_NUM*2*i]; 
    newreg->rb = output_ptr[4+FPGA_RET_PARAM_NUM*2*i] + tasks[task_idx]->rBeg;
    newreg->qe = output_ptr[3+FPGA_RET_PARAM_NUM*2*i] + tasks[task_idx]->qBeg + tasks[task_idx]->seedLength;
    newreg->re = output_ptr[5+FPGA_RET_PARAM_NUM*2*i] + tasks[task_idx]->rBeg + tasks[task_idx]->seedLength;
    newreg->score = output_ptr[6+FPGA_RET_PARAM_NUM*2*i]; 
    newreg->truesc = output_ptr[7+FPGA_RET_PARAM_NUM*2*i]; 
    newreg->w = output_ptr[8+FPGA_RET_PARAM_NUM*2*i];
    // compute the seed cov
    newreg->seedcov=0;  // TODO:add the seedcov compute function
    for (int j = 0; j < tasks[i]->chain->n; ++j){
      const mem_seed_t *t = &tasks[i]->chain->seeds[j];
      if (t->qbeg >= newreg->qb && 
          t->qbeg + t->len <= newreg->qe && 
          t->rbeg >= newreg->rb && 
          t->rbeg + t->len <= newreg->re){
        newreg->seedcov += t->len; 
      }
    }

    kv_push(mem_alnreg_t,av[regs_idx],*newreg);
    tasks[i]->read_obj->finish();
  }
  
  fprintf(stderr, "FPGA output used %dus\n", blaze::getUs()-start_ts); 
  delete [] data_ptr;
  delete [] output_ptr;
}

#define MAX_BAND_TRY  2

void extendOnCPU(
    std::vector<ExtParam*> &tasks,
    int numoftask,
    mem_opt_t *opt,
    mem_alnreg_v *av) 
{
  for (int i = 0; i < numoftask; i++) {
    int aw[2];int max_off[2];
    int tmpidx = tasks[i]->idx;
    aw[0]=opt->w;
    aw[1]=opt->w;

    mem_alnreg_t* newreg = tasks[i]->newreg;

    if(tasks[i]->qBeg){
      int qle,tle,gtle,gscore;
      for (int j=0;j<MAX_BAND_TRY;j++){
        int prev = newreg->score;
        aw[0] = opt->w<<j;
        newreg->score = ksw_extend2(
            tasks[i]->leftQlen,
            tasks[i]->leftQs,
            tasks[i]->leftRlen,
            tasks[i]->leftRs,
            5, opt->mat,
            tasks[i]->oDel,
            tasks[i]->eDel,
            tasks[i]->oIns,
            tasks[i]->eIns,
            aw[0],
            tasks[i]->penClip5,
            tasks[i]->zdrop,
            tasks[i]->h0,
            &qle, &tle, &gtle,
            &gscore, &max_off[0]);
        if(newreg->score == prev||max_off[0]<(aw[0]>>1)+(aw[0]>>2)) break;
      }
      // local extension
      if (gscore <= 0 || gscore <= newreg->score - opt->pen_clip5) { 
        newreg->qb = tasks[i]->qBeg - qle;
        newreg->rb = tasks[i]->rBeg - tle;
        newreg->truesc = newreg->score;
      }
      else { // to-end extension
        newreg->qb =0 ;
        newreg->rb =tasks[i]->rBeg -gtle;
        newreg->truesc = gscore;
      }
    }
    else {
      newreg->score = newreg->truesc = tasks[i]->h0;
      newreg->qb = 0;
      newreg->rb = tasks[i]->rBeg;
    }
    if (tasks[i]->rightQlen) {
      int qle,tle,gtle,gscore,sc0 = newreg->score;

      for (int j = 0; j < MAX_BAND_TRY; j++) {
        int prev = newreg->score;
        aw[1] = opt->w<<j;
        newreg->score = ksw_extend2(
            tasks[i]->rightQlen,
            tasks[i]->rightQs,
            tasks[i]->rightRlen,
            tasks[i]->rightRs,
            5, opt->mat,
            tasks[i]->oDel,
            tasks[i]->eDel,
            tasks[i]->oIns,
            tasks[i]->eIns,
            aw[0],
            tasks[i]->penClip5,
            tasks[i]->zdrop,
            sc0,
            &qle, &tle, &gtle,
            &gscore, &max_off[1]);
        if (newreg->score == prev ||
            max_off[1] < (aw[1]>>1) + (aw[1]>>2)) 
          break;
      }
      if (gscore <= 0 || gscore <= newreg->score - opt->pen_clip5) {
        newreg->qe = tasks[i]->l_query - tasks[i]->rightQlen + qle;
        newreg->re = tasks[i]->rBeg + tasks[i]->seedLength + tle;
        newreg->truesc += newreg->score-sc0;
      }
      else{
        newreg->qe = tasks[i]->l_query;
        newreg->re = tasks[i]->rBeg + tasks[i]->seedLength + gtle;
        newreg->truesc +=gscore-sc0;
      }
    }
    else{
      newreg->qe = tasks[i]->l_query;
      newreg->re = tasks[i]->rBeg + tasks[i]->seedLength;
    }
    newreg->w = aw[0] > aw[1]? aw[0] : aw[1];
    // compute the seed cov
    newreg->seedcov=0;  // TODO:add the seedcov compute function

    for (int j = 0; j < tasks[i]->chain->n; ++j){
      const mem_seed_t *t = &tasks[i]->chain->seeds[j];
      if (t->qbeg >= newreg->qb && 
          t->qbeg + t->len <= newreg->qe && 
          t->rbeg >= newreg->rb && 
          t->rbeg + t->len <= newreg->re){
        newreg->seedcov += t->len; 
      }
    }
    // Add current aligned region to corresponding reads
    kv_push(mem_alnreg_t, av[tmpidx], *newreg);
    
    // Notify the read that the current task is finished
    tasks[i]->read_obj->finish();
  }
}

static inline int cal_max_gap(
    const mem_opt_t *opt, 
    int qlen
) {
	int l_del = (int)((double)(qlen * opt->a - opt->o_del) / opt->e_del + 1.);
	int l_ins = (int)((double)(qlen * opt->a - opt->o_ins) / opt->e_ins + 1.);
	int l = l_del > l_ins? l_del : l_ins;
	l = l > 1? l : 1;
	return l < opt->w<<1? l : opt->w<<1;
}

void prepareChainRef(
    const ktp_aux_t *aux,
    const bseq1_t *seqs,
    const mem_chain_v *chains,
    mem_chainref_t** &chain_refs,
    int batch_num)
{
  int64_t l_pac = aux->idx->bns->l_pac;

  chain_refs = (mem_chainref_t**)malloc(batch_num*sizeof(mem_chainref_t*));

  for (int i = 0; i < batch_num; i++) {

    // input for SmithWaterman for each chain
    chain_refs[i] = (mem_chainref_t*)malloc(chains[i].n*sizeof(mem_chainref_t));

    for (int j = 0; j < chains[i].n; j++) {
      // Prepare the maxspan and rseq for each seed
      mem_chain_t *c = &chains[i].a[j]; 

      int64_t rmax[2], tmp, max = 0;
      uint8_t *rseq = 0;
      uint64_t *srt;

      if (c->n == 0) {
        //chain_refs[i][j].rmax[0] = rmax[0];
        //chain_refs[i][j].rmax[1] = rmax[1];
        //chain_refs[i][j].rseq = rseq;
        //chain_refs[i][j].srt = srt;
        
        continue;
      }
      // get the max possible span
      rmax[0] = l_pac<<1; rmax[1] = 0;
      for (int k = 0; k < c->n; ++k) {
        int64_t b, e;
        const mem_seed_t *t = &c->seeds[k];
        b = t->rbeg - (t->qbeg + cal_max_gap(aux->opt, t->qbeg));
        e = t->rbeg + t->len + ((seqs[i].l_seq - t->qbeg - t->len)
            + cal_max_gap(aux->opt, seqs[i].l_seq - t->qbeg - t->len));
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

      chain_refs[i][j].rmax[0] = rmax[0];
      chain_refs[i][j].rmax[1] = rmax[1];
      chain_refs[i][j].rseq = rseq;
      chain_refs[i][j].srt = srt;
    } 
  } 
}

void freeTasks(std::vector<ExtParam*> &sw_task_v, int task_num) {

  for (int i = 0; i < task_num; i++) {
    delete [] sw_task_v[i]->leftQs;
    delete [] sw_task_v[i]->leftRs;
    delete [] sw_task_v[i]->rightQs;
    delete [] sw_task_v[i]->rightRs;
    delete sw_task_v[i];
  }
}

void freeChainRef(
    const mem_chain_v *chains, 
    mem_chainref_t **chain_refs, 
    int batch_num) 
{
   for (int i = 0 ;i < batch_num; i++) {
     for (int j = 0;j < chains[i].n; j++) {
        free(chain_refs[i][j].srt);
        free(chain_refs[i][j].rseq);
        //free(&preresult[i][j]) ;   
     }
     free(chain_refs[i]);
   }  
}

void mem_chain2aln_hw(
    ktp_aux_t *aux,
    const bseq1_t *seqs,
    const mem_chain_v* chains,
    mem_alnreg_v *av,
    int batch_num)
{
  // For statistics
  uint64_t swFPGA_time = 0;
  uint64_t extCPU_time = 0;
  int      swFPGA_num  = 0;
  int      extCPU_num  = 0;

  uint64_t start_ts = blaze::getUs();

  // Initialize per-chain references for SmithWaterman
  mem_chainref_t **chain_refs;
  prepareChainRef(aux, seqs, chains, chain_refs, batch_num);

  mem_alnreg_t* newregs = new mem_alnreg_t[batch_num];

  // Initialize output aligned regions
  for (int i=0;i<batch_num; i++){
    kv_init(av[i]);
  }

  // Initialize batch of SW tasks
  std::vector<ExtParam*> sw_task_v(chunk_size);

  // Initialize batch of SWRead objects
  std::list<SWRead*> read_batch;
  for (int i = 0; i < batch_num; i++) {
    SWRead *SWRead_ptr = new SWRead(i, aux, 
        seqs+i, chain_refs[i], chains+i, av+i);

    read_batch.push_back(SWRead_ptr); 
  }

  fprintf(stderr, "Preparation takes %dus\n", blaze::getUs()-start_ts);

  int task_num = 0;
  while (!read_batch.empty()) {
    
    std::list<SWRead*>::iterator iter = read_batch.begin();
    while (iter != read_batch.end()) {

      uint64_t start_ts;
      ExtParam* param_task;

      switch ((*iter)->nextTask(param_task, newregs+(*iter)->index())) {
        case SWRead::TaskStatus::Successful:
          sw_task_v[task_num] = param_task;
          task_num++;
          if (task_num >= chunk_size) {
            start_ts = blaze::getUs();
            #ifdef USE_FPGA
            SwFPGA(sw_task_v,task_num,av);
            #else
            extendOnCPU(sw_task_v, task_num, aux->opt, av);
            #endif
            swFPGA_time += blaze::getUs() - start_ts;
            swFPGA_num ++;

            freeTasks(sw_task_v, task_num);
            task_num = 0;
          }
          iter++;
          break;

        case SWRead::TaskStatus::Pending:
          // No more tasks, must do extend before proceeding
          start_ts = blaze::getUs();
          extendOnCPU(sw_task_v, task_num, aux->opt, av);
          extCPU_time += blaze::getUs() - start_ts;
          extCPU_num ++;

          freeTasks(sw_task_v, task_num);
          task_num = 0;
          break;

        case SWRead::TaskStatus::Finished:
          // Read is finished, remove from batch
          iter = read_batch.erase(iter);
          break;

        default: ;
      }
    }
  }
  fprintf(stderr, "%d tasks is batched, %d is not\n", swFPGA_num, extCPU_num);
  fprintf(stderr, "Batched tasks takes %dus\n", swFPGA_time);
  fprintf(stderr, "Normal tasks takes %dus\n", extCPU_time);

  delete [] newregs;
  freeChainRef(chains, chain_refs, batch_num);  
}
