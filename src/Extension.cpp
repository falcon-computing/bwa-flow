#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <vector>
#include <queue>
#include <list>
#include <boost/smart_ptr.hpp>
#include <boost/thread/thread.hpp>

#include "bwa_wrapper.h"
#include "Extension.h"
#include "FPGAAgent.h"
#include "SWTask.h"
#include "SWRead.h"
#include "util.h"

#define FPGA_RET_PARAM_NUM 5
#define USE_FPGA
//#define SMITHWATERMAN_SIM

#ifdef SMITHWATERMAN_SIM
// hw data structures
extern "C"{
void sw_top (int *a, int *output_a, int __inc);
}
#endif

extern FPGAAgent* agent;

void packData(
    int stage_cnt,
    ExtParam** &tasks,
    int batch_num,
    mem_opt_t *opt)
{
  int Buf1Len = 32 + 32*batch_num;
  int data_size = Buf1Len >> 2;
  for (int i = 0; i <batch_num ; i++){
    data_size +=((((tasks[i]->leftQlen + 
                    tasks[i]->leftRlen + 
                    tasks[i]->rightQlen + 
                    tasks[i]->rightRlen) + 1) / 
                2)+3)/4; 
  } 
  char* buf1= new char[data_size << 2];
  //------------------ store the public options at the beginning-----------------
  buf1[0] = (char)opt->o_del;
  buf1[1] = (char)opt->e_del;
  buf1[2] = (char)opt->o_ins;
  buf1[3] = (char)opt->e_ins;
  buf1[4] = (char)opt->pen_clip5;
  buf1[5] = (char)opt->pen_clip3;
  buf1[6] = (char)opt->w;
  *(int*)(&buf1[8])= batch_num;
  int i = 0;
  int LeftMaxIns = 0;
  int LeftMaxDel = 0;
  int RightMaxIns = 0;
  int RightMaxDel = 0;
  int TaskPos = 0 ;
  TaskPos = Buf1Len >> 2;
  int buf1idx = 32;

  for (int i = 0; i <batch_num ; i++)
  {
    *((short*)(&buf1[buf1idx]))= (short)tasks[i]->leftQlen; buf1idx +=2;
    *((short*)(&buf1[buf1idx]))= (short)tasks[i]->leftRlen; buf1idx +=2;
    *((short*)(&buf1[buf1idx]))= (short)tasks[i]->rightQlen; buf1idx +=2;
    *((short*)(&buf1[buf1idx]))= (short)tasks[i]->rightRlen; buf1idx +=2;
    *((int*)(&buf1[buf1idx]))= TaskPos; buf1idx +=4;
    TaskPos += ((((tasks[i]->leftQlen + tasks[i]->leftRlen + tasks[i]->rightQlen + tasks[i]->rightRlen)+1)/2)+3)/4;
    *((short*)(&buf1[buf1idx]))= (short)tasks[i]->regScore; buf1idx +=2;
    *((short*)(&buf1[buf1idx]))= (short)tasks[i]->qBeg; buf1idx +=2;
    *((short*)(&buf1[buf1idx]))= (short)tasks[i]->h0; buf1idx +=2;
    *((short*)(&buf1[buf1idx]))= (short)tasks[i]->idx; buf1idx +=2;
    LeftMaxIns = (short)((double)(tasks[i]->leftQlen* 1 + opt->pen_clip5 -opt->o_ins)/opt->e_ins+1);
    LeftMaxDel = (short)((double)(tasks[i]->leftQlen* 1 + opt->pen_clip5 -opt->o_del)/opt->e_del+1);
    RightMaxIns = (short)((double)(tasks[i]->rightQlen *1 + opt->pen_clip3 -opt->o_ins)/opt->e_ins+1);
    RightMaxDel = (short)((double)(tasks[i]->rightQlen* 1 + opt->pen_clip3 -opt->o_del)/opt->e_del+1);         // 1 stands for tasks[i]->mat.max
    *((short*)(&buf1[buf1idx]))= LeftMaxIns; buf1idx +=2;
    *((short*)(&buf1[buf1idx]))= LeftMaxDel; buf1idx +=2;
    *((short*)(&buf1[buf1idx]))= RightMaxIns; buf1idx +=2;
    *((short*)(&buf1[buf1idx]))= RightMaxDel; buf1idx +=2;
    *((int*)(&buf1[buf1idx]))= i; buf1idx +=4;
  }
  i = 0;
  int j = 0;
  int TmpIntVar = 0;
  int counter8 = 0;

  while(i < batch_num) {
    if(tasks[i]->leftQlen > 0) {
      j = 0;
      while(j < tasks[i]->leftQlen) {
        counter8 = counter8 + 1;
        TmpIntVar = TmpIntVar <<4 | tasks[i]->leftQs[j] ;
        if(counter8 % 8 ==0){
          *(uint32_t*) (&buf1[buf1idx])= TmpIntVar;
          buf1idx += 4;
        }
        j = j + 1;
      }
    }
    if(tasks[i]->rightQlen > 0) {
      j = 0;
      while(j < tasks[i]->rightQlen) {
        counter8 = counter8 + 1;
        TmpIntVar = TmpIntVar <<4 | tasks[i]->rightQs[j] ;
        if(counter8 % 8 ==0){
          *(uint32_t*) (&buf1[buf1idx])= TmpIntVar;
          buf1idx += 4;
        }
        j = j + 1;
      }
    }
    if(tasks[i]->leftRlen > 0) {
      j = 0;
      while(j < tasks[i]->leftRlen) {
        counter8 = counter8 + 1;
        TmpIntVar = TmpIntVar <<4 | tasks[i]->leftRs[j] ;
        if(counter8 % 8 ==0){
          *(uint32_t*) (&buf1[buf1idx])= TmpIntVar;
          buf1idx += 4;
        }
        j = j + 1;
      }
    }
    if(tasks[i]->rightRlen > 0) {
      j = 0;
      while(j < tasks[i]->rightRlen) {
        counter8 = counter8 + 1;
        TmpIntVar = TmpIntVar <<4 | tasks[i]->rightRs[j] ;
        if(counter8 % 8 ==0){
          *(uint32_t*) (&buf1[buf1idx])= TmpIntVar;
          buf1idx += 4;
        }
        j = j + 1;
      }
    }
    if(counter8 %8 != 0) {
      while(counter8 %8 != 0 ) {
        TmpIntVar = TmpIntVar << 4;
        counter8 = counter8 + 1;
      }
      *(uint32_t*) (&buf1[buf1idx])= TmpIntVar;
      buf1idx += 4;
    }
    delete [] tasks[i]->leftQs;
    delete [] tasks[i]->leftRs;
    i = i + 1;
  }
  agent->writeInput(buf1, data_size*sizeof(int), stage_cnt);

  delete [] buf1;
}

void SwFPGA(
    int stage_cnt,
    ExtParam** &tasks,
    int batch_num,
    mem_opt_t *opt)
{
  short* output_ptr = new short[FPGA_RET_PARAM_NUM*batch_num*2];
  uint64_t start_ts = getUs();

  agent->start(batch_num, stage_cnt);
  agent->readOutput(output_ptr, FPGA_RET_PARAM_NUM*batch_num*4, stage_cnt);

  //fprintf(stderr, "FPGA kernel used %dus\n", getUs()-start_ts); 

  start_ts = getUs();
  for (int i = 0; i < batch_num; i++) {  
    
    int task_idx = ((int)(output_ptr[1+FPGA_RET_PARAM_NUM*2*i])<<16) |
                    output_ptr[0+FPGA_RET_PARAM_NUM*2*i];
    int regs_idx = tasks[task_idx]->idx;

    mem_alnreg_t *newreg = tasks[task_idx]->newreg;

    newreg->qb = output_ptr[2+FPGA_RET_PARAM_NUM*2*i]; 
    newreg->rb = output_ptr[4+FPGA_RET_PARAM_NUM*2*i] + tasks[task_idx]->rBeg;
    newreg->qe = output_ptr[3+FPGA_RET_PARAM_NUM*2*i] + tasks[task_idx]->qBeg + tasks[task_idx]->seedLength;
    newreg->re = output_ptr[5+FPGA_RET_PARAM_NUM*2*i] + tasks[task_idx]->rBeg + tasks[task_idx]->seedLength;
    newreg->score = output_ptr[6+FPGA_RET_PARAM_NUM*2*i]; 
    newreg->truesc = output_ptr[7+FPGA_RET_PARAM_NUM*2*i]; 
    newreg->w = output_ptr[8+FPGA_RET_PARAM_NUM*2*i];

    // compute the seed cov
    newreg->seedcov=0;  // TODO:add the seedcov compute function
    for (int j = 0; j < tasks[task_idx]->chain->n; ++j){
      const mem_seed_t *t = &tasks[task_idx]->chain->seeds[j];
      if (t->qbeg >= newreg->qb && 
          t->qbeg + t->len <= newreg->qe && 
          t->rbeg >= newreg->rb && 
          t->rbeg + t->len <= newreg->re){
        newreg->seedcov += t->len; 
      }
    }
    tasks[task_idx]->read_obj->finish();
    delete tasks[task_idx];
  }
 
  //fprintf(stderr, "FPGA output used %dus\n", getUs()-start_ts); 
  delete [] output_ptr;
}

void extendOnCPU(
    ExtParam** tasks,
    int numoftask,
    mem_opt_t *opt)
{
  int aw[2];
  int max_off[2];
  aw[0] = opt->w;
  aw[1] = opt->w;

  for (int i = 0; i < numoftask; i++) {
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
            opt->o_del,
            opt->e_del,
            opt->o_ins,
            opt->e_ins,
            aw[0],
            opt->pen_clip5,
            opt->zdrop,
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
      int qle,tle,gtle,gscore;
      int sc0 = newreg->score;

      for (int j = 0; j < MAX_BAND_TRY; j++) {
        int prev = newreg->score;
        aw[1] = opt->w<<j;
        newreg->score = ksw_extend2(
            tasks[i]->rightQlen,
            tasks[i]->rightQs,
            tasks[i]->rightRlen,
            tasks[i]->rightRs,
            5, opt->mat,
            opt->o_del,
            opt->e_del,
            opt->o_ins,
            opt->e_ins,
            aw[0],
            opt->pen_clip3,
            opt->zdrop,
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
    // Notify the read that the current task is finished
    tasks[i]->read_obj->finish();

    // free the tasks
    delete [] tasks[i]->leftQs;
    delete [] tasks[i]->leftRs;
    delete tasks[i];
  }
}

#if 0
void mem_chain2aln_hw(
    ktp_aux_t *aux,
    const bseq1_t *seqs,
    const mem_chain_v* chains,
    mem_alnreg_v *av,
    int batch_num)
{
  // For statistics
  uint64_t getTask_time = 0;
  uint64_t swFPGA_time  = 0;
  uint64_t extCPU_time  = 0;
  int      swFPGA_num   = 0;
  int      extCPU_num   = 0;
  int      join_num     = 0;


  // Initialize output aligned regions
  for (int i=0;i<batch_num; i++){
    kv_init(av[i]);
  }

  // Initialize batch of SWRead objects
  std::list<SWRead*> read_batch;
  for (int i = 0; i < batch_num; i++) {
    SWRead *read_ptr = new SWRead(i, aux, 
        seqs+i, chains+i, av+i);

    read_batch.push_back(read_ptr); 
  }

  // Create thread for ping pong
  const int stage_num = 2;
  int stage_cnt = 0;
  //boost::thread_group stage_workers;
  std::queue<boost::shared_ptr<boost::thread> > stage_workers;

  // Initialize batch of SW tasks
  int task_num[stage_num] = {0};
  ExtParam **sw_task_v[stage_num];
  for (int i = 0; i < stage_num; i++) {
    task_num[i] = 0;
    sw_task_v[i] = new ExtParam*[chunk_size];
  }

  uint64_t start_ts = getUs();
  uint64_t last_ts;
  while (!read_batch.empty()) {
    
    std::list<SWRead*>::iterator iter = read_batch.begin();
    while (iter != read_batch.end()) {

      uint64_t start_ts;
      ExtParam* param_task;
      uint64_t nt_time = getUs();
      SWRead::TaskStatus status = (*iter)->nextTask(param_task);
      getTask_time += getUs()-nt_time; 

      switch (status) {

        case SWRead::TaskStatus::Successful:
          sw_task_v[stage_cnt][task_num[stage_cnt]] = param_task;
          task_num[stage_cnt]++;
          if (task_num[stage_cnt] >= chunk_size) {
            start_ts = getUs();
#ifdef USE_FPGA
            uint64_t pd_ts = getUs();
            packData(stage_cnt,
                sw_task_v[stage_cnt],
                task_num[stage_cnt],
                aux->opt);

            if (stage_workers.size() >= stage_num - 1) {
              stage_workers.front()->join();
              stage_workers.pop();
            }
            fprintf(stderr, "the time since last kernal call is %dus\n", getUs()-last_ts);
            fprintf(stderr, "the time used in nextTask is %dus\n", getTask_time);
            fprintf(stderr, "the time used in packData is %dus\n", getUs()-pd_ts);        

            last_ts = getUs();
            getTask_time = 0 ;
            boost::shared_ptr<boost::thread> worker(new 
                boost::thread(&SwFPGA,
                  stage_cnt,
                  sw_task_v[stage_cnt],
                  task_num[stage_cnt],
                  aux->opt));
#else
            if (stage_workers.size() >= stage_num - 1) {
              stage_workers.front()->join();
              stage_workers.pop();
            }
            boost::shared_ptr<boost::thread> worker(new 
                boost::thread(&extendOnCPU,
                  sw_task_v[stage_cnt],
                  task_num[stage_cnt],
                  aux->opt));
#endif
            stage_workers.push(worker);

            task_num[stage_cnt] = 0;

            stage_cnt = (stage_cnt + 1) % stage_num;
            swFPGA_time += getUs() - start_ts;
            swFPGA_num ++;
          }
          iter++;
          break;

        case SWRead::TaskStatus::Pending:
          start_ts = getUs();
          if (!stage_workers.empty()) {
            stage_workers.front()->join();
            stage_workers.pop();
            join_num ++;
          }
          else {
            // No more tasks, must do extend before proceeding
            extendOnCPU(sw_task_v[stage_cnt],
                task_num[stage_cnt],
                aux->opt);

            task_num[stage_cnt] = 0;
            iter++;

            extCPU_num ++;
          }
          extCPU_time += getUs() - start_ts;
          break;

        case SWRead::TaskStatus::Finished:
          // Read is finished, remove from batch
          delete *iter;
          iter = read_batch.erase(iter);
          break;

        default: ;
      }
    }
  }
  for (int i = 0; i < stage_num; i++) {
    delete [] sw_task_v[i];
  }
  fprintf(stderr, "%d force join tasks\n", join_num);
  fprintf(stderr, "%d Batched tasks takes %dus\n", swFPGA_num, swFPGA_time);
  fprintf(stderr, "%d Normal tasks takes %dus\n", extCPU_num, extCPU_time);
  fprintf(stderr, "All tasks takes %dus\n", getUs()-start_ts);
}
#endif