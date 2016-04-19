#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <vector>

#include "bwa/bwamem.h"
#include "bwa/bntseq.h"
#include "bwa/ksw.h"
#include "bwa/kvec.h"
#include "bwa/ksort.h"
#include "bwa/utils.h"
#include "bwa_wrapper.h"
#include "blaze/AccAgent.h"

#define FPGA_RET_PARAM_NUM 5
#define chunk_size 2000

//#define USE_FPGA

// hw data structures
//extern "C"{
//void sw_top (int *a, int *output_a, int __inc);
//}

class ExtParam
{
  public:
    uint8_t *leftQs;
    int     leftQlen ;
    uint8_t *leftRs;
    int     leftRlen ;
    uint8_t *rightQs;
    int     rightQlen ;
    uint8_t *rightRs;
    int     rightRlen ;
    int     w ;
    int8_t *mat;
    int     oDel ;
    int     eDel ;
    int     oIns ;
    int     eIns ;
    int     penClip5 ;
    int     penClip3 ;
    int     zdrop   ;
    int     h0 ;
    int     regScore ;
    int     qBeg    ;
    int     idx    ;
    int64_t     rBeg   ;    // just for testing on cpu
    int     seedLength ;
    int     l_query ;
};

class ExtRet
{
  public:
    int qBeg;
    long rBeg;
    int qEnd;
    long rEnd;
    int score;
    int trueScore;
    int width;
    int idx;
};


class increRes{
public:
  bool isfinished;
  int start; 
  int end;
  increRes(bool _isfinished,int _start,int _end): 
    isfinished(_isfinished),start(_start),end(_end)
    {};
};


class preResultofSw{
public:
   int64_t *rmax;
   uint8_t *rseq; uint64_t *srt;
   preResultofSw()
   {
     rmax = 0;
     rseq = 0;
     srt = 0;
   }
   preResultofSw(int64_t *_rmax,uint8_t *_rseq,uint64_t *_srt)
    :rseq(_rseq),srt(_srt)
   {
     rmax = new int64_t[2];
     rmax[0]= _rmax[0];
     rmax[1]=_rmax[1];
   }
   preResultofSw(const preResultofSw& preresult)
   {
     rmax = preresult.rmax;
     rseq = preresult.rseq;
     srt = preresult.srt;
   }
  /* ~preResultofSw()
   {
     free (rmax);
     free (rseq);
     free (srt);
   }  */ 
};

// The GetTask function get the extension parameters for each seed
static inline void GetTask(
    const mem_seed_t *seed, 
    mem_opt_t *opt,
    std::vector<ExtParam*> &sw_task_v,
    const uint8_t *query,
    int l_query,
    int64_t rmax_0, 
    int64_t rmax_1,
    uint8_t *rseq,
    int idx, 
    int &taskidx
) {
  if (seed->qbeg > 0 ||
      seed->qbeg + seed->len != l_query) {

    ExtParam *SwTask = sw_task_v[taskidx];

    int i = 0;
    SwTask->leftQlen = seed->qbeg;
    if(SwTask->leftQlen > 0)
    {
      SwTask->leftQs = new uint8_t[SwTask->leftQlen];
      for(i = 0;i < SwTask->leftQlen; i++)
        SwTask->leftQs[i] = query[SwTask->leftQlen-1-i];
      SwTask->leftRlen =(int)( seed->rbeg - rmax_0) ;
      SwTask->leftRs = new uint8_t[SwTask->leftRlen];
      for(i = 0; i<SwTask->leftRlen; i++)
        SwTask->leftRs[i] = rseq[SwTask->leftRlen-1-i];
    }
    else
    {
      SwTask->leftQs = NULL;
      SwTask->leftRlen = 0;
      SwTask->leftRs  = NULL;
    }

    int qe = seed->qbeg + seed->len;
    SwTask->rightQlen = l_query - qe;
    if(SwTask->rightQlen > 0)
    {
      SwTask->rightQs   = new uint8_t[SwTask->rightQlen];
      for(int i = 0;i<SwTask->rightQlen;i++)
        SwTask->rightQs[i] = query[i+qe];
      int64_t re = seed->rbeg + seed->len - rmax_0;
      SwTask->rightRlen =(int) (rmax_1 - rmax_0 -re) ;
      SwTask->rightRs = new uint8_t[SwTask->rightRlen];
      for(int i = 0;i<SwTask->rightRlen; i++)
        SwTask->rightRs[i] = rseq[i+(int)re];
    }
    else
    {
      SwTask->rightQs = NULL;
      SwTask->rightRlen = 0;
      SwTask->rightRs = NULL;
    }
    SwTask->w = opt->w ;
    SwTask->mat = opt->mat ;
    SwTask->oDel = opt->o_del ;
    SwTask->oIns = opt->o_ins ;
    SwTask->eDel = opt->e_del ;
    SwTask->eIns = opt->e_ins ;
    SwTask->penClip5 = opt->pen_clip5 ;
    SwTask->penClip3 = opt->pen_clip3 ;
    SwTask->zdrop = opt->zdrop ;
    SwTask->h0 = seed->len*opt->a ;
    SwTask->regScore = seed->len*opt->a ;
    SwTask->qBeg = seed->qbeg ;
    SwTask->rBeg = seed->rbeg ;     // for testing
    SwTask->seedLength = seed->len ;
    SwTask->idx = idx ;
    SwTask->l_query = l_query ;

    //sw_task_v[taskidx] = SwTask;

    taskidx += 1;
  }
}

// update the coordinates
static void updateCoordinates(
    int coordinates[][2],
    const MemChainVector *chains,
    int i
) {
  if(coordinates[i][1]>0){
    coordinates[i][1]-=1;
  }   
  else if(coordinates[i][1]==0){
    if(coordinates[i][0]==chains[i].n-1){
      coordinates[i][1]= -1;
    }
    else{
      coordinates[i][0]+=1;
      coordinates[i][1]=chains[i].a[coordinates[i][0]].n-1; 
    }
  }
}

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
    mem_alnreg_t *newregs,
    const MemChainVector *chains,
    int coordinates[][2],
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
//  int buf1idx = 8;
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
    LeftMaxIns = (int)((double)(tasks[i]->leftQlen*1+ tasks[i]->penClip5 -tasks[i]->oIns)/tasks[i]->eIns+1);
    LeftMaxDel = (int)((double)(tasks[i]->leftQlen*1 + tasks[i]->penClip5 -tasks[i]->oDel)/tasks[i]->eDel+1);
    RightMaxIns = (int)((double)(tasks[i]->rightQlen*1 + tasks[i]->penClip3 -tasks[i]->oIns)/tasks[i]->eIns+1);
    RightMaxDel = (int)((double)(tasks[i]->rightQlen*1 + tasks[i]->penClip3 -tasks[i]->oDel)/tasks[i]->eDel+1);         // 1 stands for tasks[i]->mat.max
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(LeftMaxIns));
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(LeftMaxDel));
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(RightMaxIns));
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(RightMaxDel));    
    //buf1idx = Int2CharArray(buf1,buf1idx,tasks[i]->idx);
    buf1idx = Int2CharArray(buf1, buf1idx, i);
    i = i+1;
  }

  char *buf2 = new char[(TaskPos<<2)-Buf1Len];
  int input_length = TaskPos<<2;
  int buf2idx = 0;
  i = 0;
  int j = 0;
  int TmpIntVar = 0;
  int Counter8 = 0;
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

  int *data_ptr = new int[input_length/4];
  short* output_ptr = new short[FPGA_RET_PARAM_NUM*batch_num*2];

  memcpy(data_ptr,buf1,Buf1Len*sizeof(char) );
  memcpy(&data_ptr[Buf1Len/4],buf2,input_length-Buf1Len);                          
  delete buf1;
  delete buf2;             
  //fprintf(stderr, "FPGA preparation used %dus\n", blaze::getUs()-start_ts); 

  // sw_top (data_ptr, (int *)output_ptr,batch_num);
  
  start_ts = blaze::getUs();
  blaze::Task_ptr fpga_task = agent->createTask(acc_id);
  if (!fpga_task) {
    throw blaze::internalError("Task is not created");
  }
  agent->writeInput(fpga_task, acc_id, data_ptr, 1, input_length/4, sizeof(int));
  agent->writeInput(fpga_task, acc_id, &batch_num, 1, 1, sizeof(int));
  agent->readOutput(fpga_task, output_ptr, FPGA_RET_PARAM_NUM*batch_num*4);
  //fprintf(stderr, "FPGA kernel used %dus\n", blaze::getUs()-start_ts); 

  start_ts = blaze::getUs();
  for (int i = 0; i < batch_num; i++) {  
    
    int task_idx = ((int)(output_ptr[1+FPGA_RET_PARAM_NUM*2*i])<<16) |
                    output_ptr[0+FPGA_RET_PARAM_NUM*2*i];
    int regs_idx = tasks[task_idx]->idx;

    newregs[regs_idx].qb = output_ptr[2+FPGA_RET_PARAM_NUM*2*i]; 
    newregs[regs_idx].rb = output_ptr[4+FPGA_RET_PARAM_NUM*2*i] + tasks[task_idx]->rBeg;
    newregs[regs_idx].qe = output_ptr[3+FPGA_RET_PARAM_NUM*2*i] + tasks[task_idx]->qBeg + tasks[task_idx]->seedLength;
    newregs[regs_idx].re = output_ptr[5+FPGA_RET_PARAM_NUM*2*i] + tasks[task_idx]->rBeg + tasks[task_idx]->seedLength;
    newregs[regs_idx].score = output_ptr[6+FPGA_RET_PARAM_NUM*2*i]; 
    newregs[regs_idx].truesc = output_ptr[7+FPGA_RET_PARAM_NUM*2*i]; 
    newregs[regs_idx].w = output_ptr[8+FPGA_RET_PARAM_NUM*2*i];
    // compute the seed cov
    newregs[regs_idx].seedcov=0;  // TODO:add the seedcov compute function
    for (int j = 0; j < chains[regs_idx].a[coordinates[regs_idx][0]].n; ++j){
          const mem_seed_t *t =&chains[regs_idx].a[coordinates[regs_idx][0]].seeds[j];
          if(t->qbeg >= newregs[regs_idx].qb && 
             t->qbeg + t->len <= newregs[regs_idx].qe && 
             t->rbeg >= newregs[regs_idx].rb && 
             t->rbeg + t->len <= newregs[regs_idx].re){
             newregs[regs_idx].seedcov += t->len; 
          }
        }
    kv_push(mem_alnreg_t,av[regs_idx],newregs[regs_idx]);
    // increment the coordinates record
    updateCoordinates(coordinates,chains,regs_idx);
  }
  
  //fprintf(stderr, "FPGA output used %dus\n", blaze::getUs()-start_ts); 
  delete [] data_ptr;
  delete [] output_ptr;
}

#define MAX_BAND_TRY  2

void extendOnCPU(
    std::vector<ExtParam*> &tasks,
    int numoftask,
    mem_opt_t *opt,
    mem_alnreg_t *newregs,
    const MemChainVector *chains,
    int coordinates[][2],
    mem_alnreg_v *av) 
{
  for (int i = 0; i < numoftask; i++) {
    int aw[2];int max_off[2];
    int tmpidx = tasks[i]->idx;
    aw[0]=opt->w;
    aw[1]=opt->w;

    if(tasks[i]->qBeg){
      int qle,tle,gtle,gscore;
      for (int j=0;j<MAX_BAND_TRY;j++){
        int prev = newregs[tmpidx].score;
        aw[0] = opt->w<<j;
        newregs[tmpidx].score = ksw_extend2(
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
        if(newregs[tmpidx].score == prev||max_off[0]<(aw[0]>>1)+(aw[0]>>2)) break;
      }
      // local extension
      if (gscore <= 0 || gscore <= newregs[tmpidx].score - opt->pen_clip5) { 
        newregs[tmpidx].qb = tasks[i]->qBeg - qle;
        newregs[tmpidx].rb = tasks[i]->rBeg - tle;
        newregs[tmpidx].truesc = newregs[tmpidx].score;
      }
      else { // to-end extension
        newregs[tmpidx].qb =0 ;
        newregs[tmpidx].rb =tasks[i]->rBeg -gtle;
        newregs[tmpidx].truesc = gscore;
      }
    }
    else {
      newregs[tmpidx].score = newregs[tmpidx].truesc = tasks[i]->h0;
      newregs[tmpidx].qb = 0;
      newregs[tmpidx].rb = tasks[i]->rBeg;
    }
    if (tasks[i]->rightQlen) {
      int qle,tle,gtle,gscore,sc0 = newregs[tmpidx].score;

      for (int j = 0; j < MAX_BAND_TRY; j++) {
        int prev = newregs[tmpidx].score;
        aw[1] = opt->w<<j;
        newregs[tmpidx].score = ksw_extend2(
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
        if (newregs[tmpidx].score == prev ||
            max_off[1] < (aw[1]>>1) + (aw[1]>>2)) 
          break;
      }
      if (gscore <= 0 || gscore <= newregs[tmpidx].score - opt->pen_clip5) {
        newregs[tmpidx].qe = tasks[i]->l_query - tasks[i]->rightQlen + qle;
        newregs[tmpidx].re = tasks[i]->rBeg + tasks[i]->seedLength + tle;
        newregs[tmpidx].truesc += newregs[tmpidx].score-sc0;
      }
      else{
        newregs[tmpidx].qe = tasks[i]->l_query;
        newregs[tmpidx].re = tasks[i]->rBeg + tasks[i]->seedLength + gtle;
        newregs[tmpidx].truesc +=gscore-sc0;
      }
    }
    else{
      newregs[tmpidx].qe = tasks[i]->l_query;
      newregs[tmpidx].re = tasks[i]->rBeg + tasks[i]->seedLength;
    }
    newregs[tmpidx].w = aw[0] > aw[1]? aw[0] : aw[1];
    // compute the seed cov
    newregs[tmpidx].seedcov=0;  // TODO:add the seedcov compute function
    for (int j=0; j<chains[tmpidx].a[coordinates[tmpidx][0]].n; ++j){
      const mem_seed_t *t =&chains[tmpidx].a[coordinates[tmpidx][0]].seeds[j];
      if(t->qbeg >= newregs[tmpidx].qb && 
          t->qbeg + t->len <= newregs[tmpidx].qe && 
          t->rbeg >= newregs[tmpidx].rb && 
          t->rbeg + t->len <= newregs[tmpidx].re){
        newregs[tmpidx].seedcov += t->len; 
      }
    }
    kv_push(mem_alnreg_t,av[tmpidx],newregs[tmpidx]);
    // increment the coordinates record
    updateCoordinates(coordinates,chains,tmpidx);
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

static inline void initCoordinates(
    int coordinates[][2],
    const MemChainVector *chains,
    int numofreads
) { 
  for (int i =0;i<numofreads;i++){
    coordinates[i][0]=0;
    if(&chains[i]!=NULL&&chains[i].n>0){
      coordinates[i][1]=chains[i].a[0].n-1;
    }
    else coordinates[i][1]= -1;
  } 
}

increRes incrementCoordinates(
    int coordinates[][2],
    const MemChainVector *chains,
    int numofreads,
    int start,
    int end
) {
  bool isfinished = true;
  int curstart = start;
  int curend = end;
  while(curstart < curend && isfinished == true){
    if(&chains[curstart]==NULL) curstart +=1;
    else if(coordinates[curstart][1]<0) curstart +=1;
    else if(coordinates[curstart][1]==0 &&coordinates[curstart][0]==chains[curstart].n-1)
      curstart +=1;
    else
      isfinished = false;
  }
  bool endflag = true;
  while(curstart<curend-1&&endflag==true){
    if(&chains[curend-1]==NULL) curend -=1;
    else if(coordinates[curend-1][1]<0) curend-=1;
    else if(coordinates[curend-1][1]==0&&coordinates[curend-1][0]==chains[curend-1].n-1)
      curend-=1;
    else endflag = false;
  } 
  int i = curstart;
  while (i <curend){
    if(&chains[i]!=NULL){
      if(coordinates[i][1]>0)
         coordinates[i][1]-=1;   
      else if(coordinates[i][1]==0){
        if(coordinates[i][0]==chains[i].n-1)
          coordinates[i][1]= -1;
        else{
          coordinates[i][0]+=1;
          coordinates[i][1]=chains[i].a[coordinates[i][0]].n-1; 
        }
      }
    } 
    i = i+1;
  } 
  increRes ret(isfinished,curstart,curend);
  return ret;
}

static inline int testExtension(
    mem_opt_t *opt,
    mem_seed_t& seed,
    mem_alnreg_v& alnregv, 
    int l_query
) {
  long rdist = -1;
  int qdist = -1;  
  int maxgap = -1;
  int mindist = -1;
  int w = -1;
  int breakidx = 0;
  bool isbreak = false;

  int i;
  for (i = 0; i < alnregv.n; i++) {
    if (seed.rbeg < alnregv.a[i].rb ||
        seed.rbeg + seed.len> alnregv.a[i].re||
        seed.qbeg < alnregv.a[i].qb ||
        seed.qbeg + seed.len > alnregv.a[i].qe)
    {
      continue; 
    }
    if (seed.len - alnregv.a[i].seedlen0 > .1 * l_query){
      continue;
    }
    qdist   = seed.qbeg - alnregv.a[i].qb;
    rdist   = seed.rbeg - alnregv.a[i].rb;
    mindist = (qdist < rdist) ? qdist : (int)rdist;
    maxgap  = cal_max_gap(opt,mindist);
    w = (maxgap < alnregv.a[i].w) ? maxgap : alnregv.a[i].w;

    if (qdist - rdist < w &&
        rdist - qdist < w) 
    {
      break;
    }

    qdist   = alnregv.a[i].qe -(seed.qbeg + seed.len);
    rdist   = alnregv.a[i].re - (seed.rbeg + seed.len);
    mindist = (qdist < rdist) ? qdist : (int)rdist;
    maxgap  = cal_max_gap(opt, mindist);
    w = (maxgap < alnregv.a[i].w) ? maxgap : alnregv.a[i].w;
    if (qdist - rdist < w &&
        rdist - qdist < w) 
    {
      break;
    }
  }
  return i;
}

static inline int checkoverlap(
    int startidx,
    mem_seed_t& seed,
    mem_chain_t& chain,
    uint64_t *srt
) {
  int i;
  for (i = startidx; i < chain.n; i++) {
    const mem_seed_t* targetseed;
    if(srt[i]==0) {
      continue;
    }
    targetseed = &chain.seeds[(uint32_t)srt[i]];
    if(targetseed->len < seed.len* 0.95) {
      continue;
    }
    if(seed.qbeg <= targetseed->qbeg && 
       seed.qbeg + seed.len - targetseed->qbeg >= seed.len>>2 &&
       targetseed->qbeg - seed.qbeg != targetseed->rbeg-seed.rbeg) {
      break;
    }
    if(targetseed->qbeg <= seed.qbeg &&
       targetseed->qbeg + targetseed->len - seed.qbeg >= seed.len>>2 &&
       seed.qbeg - targetseed->qbeg != seed.rbeg - targetseed->rbeg) 
    {
      break;
    }
  }
  return i;
}

void prepareForSw(
    int batch_num,
    const MemChainVector *chains,
    const bseq1_t *seqs,
    preResultofSw **preResultofSw_m,
    ktp_aux_t *aux)
{
  int64_t l_pac = aux->idx->bns->l_pac;

  for (int i=0; i<batch_num; i++)   //loop for each seq
  {
    preResultofSw* preResultofSw_v = new preResultofSw[chains[i].n];  // preresult vector for each chain
    for (int j=0;j<chains[i].n;j++){
      mem_chain_t *c = &chains[i].a[j]; // ------------prepare the maxspan and rseq for each seed
      int64_t rmax[2], tmp, max = 0;
      uint8_t *rseq = 0;
      uint64_t *srt;
      if (c->n == 0) {
        preResultofSw preresult(rmax,rseq,srt);
        preResultofSw_v[j]= preresult; 
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
      preResultofSw preresult(rmax,rseq,srt);
      preResultofSw_v[j]=preresult;
    } 
    preResultofSw_m[i]= preResultofSw_v;
  } 
}

void freeTask(std::vector<ExtParam*> &sw_task_v, int task_num) {

  for (int i = 0; i < task_num; i++) {
    free(sw_task_v[i]->leftQs);
    free(sw_task_v[i]->leftRs);
    free(sw_task_v[i]->rightQs);
    free(sw_task_v[i]->rightRs);
  }
}

void freePreresult(
    preResultofSw **preresult, 
    const MemChainVector *chains, 
    int batch_num
) {
   for (int i =0 ;i <batch_num ;i++){
     for (int j = 0;j < chains[i].n;j++){
        free(preresult[i][j].rmax);
        free(preresult[i][j].srt);
        free(preresult[i][j].rseq);
        //free(&preresult[i][j]) ;   
     }
     free(preresult[i]);
   }  
}

void mem_chain2aln_hw(
    ktp_aux_t *aux,
    const bseq1_t *seqs,
    const MemChainVector* chains,
    mem_alnreg_v *av,
    int batch_num)
{
  // For statistics
  uint64_t swFPGA_time = 0;
  uint64_t extCPU_time = 0;
  int      swFPGA_num  = 0;
  int      extCPU_num  = 0;

  int numoftask = 0;
  int testCount = 0;
  preResultofSw **preResultofSw_m = new preResultofSw*[batch_num];  // preresult matrix for each read
  int max_off[2], aw[2]; // aw: actual bandwidth used in extension
  
  uint64_t start_ts = blaze::getUs();
  prepareForSw(batch_num,chains,seqs,preResultofSw_m,aux);    // prepare process

  int coordinates[batch_num][2];    // the chain and seed record for each read
  initCoordinates(coordinates,chains,batch_num); 
  int extensionflags;
  int overlapflags;

  int start = 0;
  int end = batch_num;
  int taskidx = 0;

  mem_alnreg_t *newregs=new mem_alnreg_t[batch_num];

  for (int i=0;i<batch_num; i++){
    kv_init(av[i]);
  }
  std::vector<ExtParam*> sw_task_v(chunk_size);
  for (int i = 0; i < chunk_size; i++) {
    sw_task_v[i] = new ExtParam;
  }
  fprintf(stderr, "Preparation takes %dus\n", blaze::getUs()-start_ts);

  int i = start;
  int old_iter_end = 0;
  int64_t total_task = 0;
  while (true) {
    // while loop to get the tasks
    while (coordinates[i][1] >= 0) {
      int chain_idx = coordinates[i][0];
      int seed_idx = coordinates[i][1];
      uint32_t sorted_idx = (uint32_t)(preResultofSw_m[i][chain_idx].srt[seed_idx]);

      // get next available seed in the current read
      mem_seed_t* seedarray = &chains[i].a[chain_idx].seeds[sorted_idx];

      // test if extension has been made before
      extensionflags = testExtension(aux->opt, *seedarray, av[i], seqs[i].l_seq);

      // Check the overlapping in the former processed seeds        
      if (extensionflags < av[i].n) {                                               
        overlapflags = checkoverlap(
            coordinates[i][1]+1,
            *seedarray,
            chains[i].a[coordinates[i][0]],
            preResultofSw_m[i][chain_idx].srt);   
      }

      // There is no need to compute this seed
      if (extensionflags < av[i].n && 
          overlapflags == chains[i].a[coordinates[i][0]].n) 
      {
        // Mark this seed as un computed
        preResultofSw_m[i][chain_idx].srt[seed_idx] = 0;         

        updateCoordinates(coordinates, chains, i);
      }
      else {
        // initialize the newregs
        newregs[i].score = seedarray->len*aux->opt->a;
        newregs[i].truesc = seedarray->len*aux->opt->a;
        newregs[i].qb = 0;
        newregs[i].rb = seedarray->rbeg;
        newregs[i].qe = seqs[i].l_seq;
        newregs[i].re = seedarray->rbeg + seedarray->len; 
        newregs[i].rid = chains[i].a[coordinates[i][0]].rid;
        newregs[i].seedlen0 = seedarray->len; 
        newregs[i].frac_rep = chains[i].a[coordinates[i][0]].frac_rep;
        newregs[i].w = aux->opt->w;

        if (seedarray->qbeg > 0 ||
            seedarray->qbeg + seedarray->len != seqs[i].l_seq)
        { 
          // need to do extension
          GetTask(
              seedarray,
              aux->opt,
              sw_task_v,
              (const uint8_t*)seqs[i].seq,
              seqs[i].l_seq,
              preResultofSw_m[i][chain_idx].rmax[0],
              preResultofSw_m[i][chain_idx].rmax[1],
              preResultofSw_m[i][chain_idx].rseq,
              i, taskidx);
        }
        else {
          // no need to extend ,just push to alnreg_v 
          newregs[i].seedcov=0;

          for (int j = 0; j < chains[i].a[coordinates[i][0]].n; j++) {
            const mem_seed_t *t =&chains[i].a[coordinates[i][0]].seeds[j];
            if (t->qbeg >= newregs[i].qb && 
                t->qbeg + t->len <= newregs[i].qe && 
                t->rbeg >= newregs[i].rb && 
                t->rbeg + t->len <= newregs[i].re) {
              newregs[i].seedcov += t->len; 
            }
          }
          kv_push(mem_alnreg_t,av[i],newregs[i]);       
          updateCoordinates(coordinates,chains,i);
        }
        // go on to next read
        break;
      }
    }
    
    if (taskidx >= chunk_size) {
      uint64_t start_ts = blaze::getUs();
#ifdef USE_FPGA
      SwFPGA(sw_task_v,taskidx,newregs,chains,coordinates,av);        
      //         fprintf(stderr,"FPGA computed %d tasks for %dus\n",taskidx,blaze::getUs()-start_ts);
      //printf("FPGA computed %d tasks for %dus\n",taskidx,blaze::getUs()-start_ts);
#else
      extendOnCPU(sw_task_v,taskidx,aux->opt,newregs,chains,coordinates,av);
      //         fprintf(stderr,"CPU computed %d tasks for %dus\n",taskidx,blaze::getUs()-start_ts);
      //printf("CPU computed %d tasks for %dus\n",taskidx,blaze::getUs()-start_ts);
#endif
      swFPGA_time += blaze::getUs() - start_ts;
      swFPGA_num ++;

      freeTask(sw_task_v, taskidx);

      total_task += taskidx;
      taskidx = 0;
      old_iter_end = i ;
      if(i == end-1){
        i = start;
        continue;
      }
      else{
        i= i+1 ;
        continue ;
      }
    }
    else if (old_iter_end == start && i == end -1) {
      if (taskidx > 0) {
        uint64_t start_ts = blaze::getUs();
        extendOnCPU(sw_task_v,taskidx,aux->opt,newregs,chains,coordinates,av);
        extCPU_time += blaze::getUs() - start_ts;
        extCPU_num ++;
        // fprintf(stderr,"CPU computed %d tasks\n",taskidx);
        freeTask(sw_task_v, taskidx);
        total_task += taskidx;
        taskidx = 0;
        old_iter_end = i ;
        i = start;
        continue;
      }
      else break;
    }
    else if (old_iter_end > start && i == old_iter_end -1) {
      if (taskidx > 0) {
        uint64_t start_ts = blaze::getUs();
        extendOnCPU(sw_task_v,taskidx,aux->opt,newregs,chains,coordinates,av);
        extCPU_time += blaze::getUs() - start_ts;
        extCPU_num ++;
        //   fprintf(stderr,"CPU computed %d tasks\n",taskidx);
        freeTask(sw_task_v, taskidx);
        total_task += taskidx;
        taskidx = 0;
        old_iter_end = i ;
        i = i+1;
      }
      else break;
    }
    else if (i == end-1) {
      i = start;
      continue;
    }
    else {
      i = i+1 ;
    }
  }
  fprintf(stderr, "%d tasks is batched, %d is not\n", swFPGA_num, extCPU_num);
  fprintf(stderr, "Batched tasks takes %dus\n", swFPGA_time);
  fprintf(stderr, "Normal tasks takes %dus\n", extCPU_time);

  for (int i = 0; i < sw_task_v.size(); i++) {
    delete sw_task_v[i];
  }
  freePreresult(preResultofSw_m,chains,batch_num);  
}
