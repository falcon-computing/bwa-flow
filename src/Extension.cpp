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
#include "Extension.h"
#include "FPGAAgent.h"
#include "SWTask.h"
#include "SWRead.h"
#include "util.h"

#define MAX_BAND_TRY  2

#ifdef SMITHWATERMAN_SIM
// hw data structures
extern "C"{
void sw_top (int *a, int *output_a, int __inc);
}
#endif

void extendOnFPGAPackInput(
    FPGAAgent* agent,
    int stage_cnt,
    ExtParam** &tasks,
    int batch_num,
    mem_opt_t *opt)
{
  uint64_t start_ts = getUs();

  int Buf1Len = 32 + 32*batch_num;
  int data_size = Buf1Len >> 2;
  for (int i = 0; i <batch_num ; i++){
    data_size +=((((tasks[i]->leftQlen + 
                    tasks[i]->leftRlen + 
                    tasks[i]->rightQlen + 
                    tasks[i]->rightRlen) + 1) / 
                2)+3)/4; 
  } 
  char* buf1 = (char*)malloc(data_size << 2);

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
  VLOG(3) << "packData takes " 
    << getUs() - start_ts << " us";

  start_ts = getUs();
  agent->writeInput(buf1, data_size*sizeof(int), stage_cnt);
  agent->start(batch_num, stage_cnt);

  delete [] buf1;
  VLOG(3) << "FPGA input takes " 
    << getUs() - start_ts << " us";
}

void extendOnFPGAProcessOutput(
    FPGAAgent* agent,
    int stage_cnt,
    ExtParam** &tasks,
    int batch_num,
    mem_opt_t *opt)
{
  short* output_ptr = new short[FPGA_RET_PARAM_NUM*batch_num*2];

  uint64_t start_ts = getUs();
  agent->wait(stage_cnt);
  VLOG(3) << "Wait for FPGA takes " 
    << getUs() - start_ts << " us";

  start_ts = getUs();
  agent->readOutput(output_ptr, FPGA_RET_PARAM_NUM*batch_num*4, stage_cnt);

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

    tasks[task_idx]->read_obj->finish();
    delete tasks[task_idx];
  }
  delete [] output_ptr;
  VLOG(3) << "FPGA output " << stage_cnt << " used " 
    << getUs() - start_ts << " us";
}

boost::mutex driver_mutex;
boost::condition_variable driver_cond;

void fpga_driver(FPGAAgent* agent)
{
   VLOG(1) << "Starting a thread working on fpga_driver";
    int pack_stage_cnt = 0;
    int stage_cnt_save = 0;
    int old_stage_cnt = 0;
    int pingpong_flag = 0;
    const int stage_num = 5;
    int chunk_size =FLAGS_chunk_size; 
    boost::unique_lock<boost::mutex> lock(driver_mutex);
    lock.unlock();
  try{
    while(true){
      while (pend_depth ==0){
            lock.lock();
            driver_cond.wait(lock);
            lock.unlock();
      }
      if(pend_flag[pack_stage_cnt]==true) {
        VLOG(3) << "Pack input on stage " << pack_stage_cnt <<" ";
        extendOnFPGAPackInput(
           agent,
           pingpong_flag,
           task_batch[pack_stage_cnt],
           chunk_size,
           aux->opt);
        stage_cnt_save = pack_stage_cnt;
      }
      pingpong_flag = 1 - pingpong_flag;
      if (agent->pending(pingpong_flag)){
            extendOnFPGAProcessOutput(
                    agent,
                    pingpong_flag,
                    task_batch[old_stage_cnt],
                    chunk_size,
                    aux->opt);
            {
              boost::lock_guard<boost::mutex> guard(driver_mutex);
              pend_flag[old_stage_cnt] =false;
              pend_depth --;
            }
            driver_cond.notify_one();
            old_stage_cnt = stage_cnt_save;
      }
      pack_stage_cnt =( pack_stage_cnt +1)%stage_num ;
    }
  }
  catch (boost::thread_interrupted &e) {
     VLOG(3) << "fpga driver stopped";
  }
}

void extendOnCPU(
    ExtParam** tasks,
    int numoftask,
    mem_opt_t *opt)
{
  uint64_t start_ts = getUs();
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

    // Notify the read that the current task is finished
    tasks[i]->read_obj->finish();

    // free the tasks
    delete [] tasks[i]->leftQs;
    delete [] tasks[i]->leftRs;
    delete tasks[i];
  }
  VLOG(3) << "SmithWaterman kernel on CPU used " 
    << getUs() - start_ts << " us";
}

