#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif
#include "bwa/bwamem.h"
#include "bwa/bntseq.h"
#include "bwa/ksw.h"
#include "bwa/kvec.h"
#include "bwa/ksort.h"
#include "bwa/utils.h"
#include "bwa_wrapper.h"
#include "SWClient.h"
#include "blaze/AccAgent.h"
#define FPGA_RET_PARAM_NUM 5
#include <vector>
#define CHECK_FPGA

blaze::AccAgent* agent;

// hw data structures
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
};

// The GetTask function get the extension parameters for each seed

void GetTask(const mem_seed_t *seed, mem_opt_t *opt,std::vector<ExtParam> *sw_task_v,
             const uint8_t *query,int l_query,int64_t rmax_0, int64_t rmax_1,uint8_t *rseq ,
             int idx, int *taskidx) 
{
  if(seed->qbeg>0||(seed->qbeg+seed->len!=l_query))
  {
    ExtParam *SwTask = new ExtParam;
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
        SwTask->rightRs[i] = rseq[i+re];
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
   (*sw_task_v).push_back(*SwTask);
    *taskidx = *taskidx +1;
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

void SwFPGA(std::vector<ExtParam> SwTask, int BatchNum, mem_alnreg_t *newregs)
{
  uint64_t start_ts = blaze::getUs();

  int Buf1Len = 32 + 32*BatchNum;
  char* buf1 = new char[Buf1Len];
  //------------------ store the public options at the beginning-----------------
  buf1[0] = (char)SwTask[0].oDel;
  buf1[1] = (char)SwTask[0].eDel;
  buf1[2] = (char)SwTask[0].oIns;
  buf1[3] = (char)SwTask[0].eIns;
  buf1[4] = (char)SwTask[0].penClip5;
  buf1[5] = (char)SwTask[0].penClip3;
  buf1[6] = (char)SwTask[0].w;
  Int2CharArray(buf1,8,BatchNum);

  //-------------------pack the batch of parameters of each SW--------------------
  int i = 0;
  int LeftMaxIns = 0;
  int LeftMaxDel = 0;
  int RightMaxIns = 0;
  int RightMaxDel = 0;
  int TaskPos = 0 ;
  TaskPos = Buf1Len >> 2;
  int buf1idx = 32;
  while(i < BatchNum)
  {
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(SwTask[i].leftQlen));
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(SwTask[i].leftRlen));
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(SwTask[i].rightQlen));
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(SwTask[i].rightRlen));
    buf1idx = Int2CharArray(buf1,buf1idx,TaskPos);
    TaskPos += ((((SwTask[i].leftQlen + SwTask[i].leftRlen + SwTask[i].rightQlen + SwTask[i].rightRlen)+1)/2)+3)/4;
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(SwTask[i].regScore));
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(SwTask[i].qBeg));
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(SwTask[i].h0));
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(SwTask[i].idx));
    LeftMaxIns = (int)((double)(SwTask[i].leftQlen*4+ SwTask[i].penClip5 -SwTask[i].oIns)/SwTask[i].eIns+1);
    LeftMaxDel = (int)((double)(SwTask[i].leftQlen*4 + SwTask[i].penClip5 -SwTask[i].oDel)/SwTask[i].eDel+1);
    RightMaxIns = (int)((double)(SwTask[i].rightQlen*4 + SwTask[i].penClip3 -SwTask[i].oIns)/SwTask[i].eIns+1);
    RightMaxDel = (int)((double)(SwTask[i].rightQlen*4 + SwTask[i].penClip3 -SwTask[i].oDel)/SwTask[i].eDel+1);         // 1 stands for SwTask[i].mat.max
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(LeftMaxIns));
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(LeftMaxDel));
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(RightMaxIns));
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(RightMaxDel));    
    //buf1idx = Int2CharArray(buf1,buf1idx,SwTask[i].idx);
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
  while(i < BatchNum)
  {
    if(SwTask[i].leftQlen > 0)
    {
      j = 0;
      while(j < SwTask[i].leftQlen)
      {
        Counter8 = Counter8 + 1;
        TmpIntVar = TmpIntVar <<4 | ((int)SwTask[i].leftQs[j] & 0x0f);
        if(Counter8 % 8 ==0)
          buf2idx = Int2CharArray(buf2,buf2idx,TmpIntVar);
        j = j + 1;
      }
    }
    if(SwTask[i].rightQlen > 0)
    {
      j = 0;
      while(j < SwTask[i].rightQlen)
      {
        Counter8 = Counter8 + 1;
        TmpIntVar = TmpIntVar <<4 | ((int)SwTask[i].rightQs[j] & 0x0f);
        if(Counter8 % 8 ==0)
          buf2idx = Int2CharArray(buf2,buf2idx,TmpIntVar);
        j = j + 1;
      }
    }
    if(SwTask[i].leftRlen > 0)
    {
      j = 0;
      while(j < SwTask[i].leftRlen)
      {
        Counter8 = Counter8 + 1;
        TmpIntVar = TmpIntVar <<4 | ((int)SwTask[i].leftRs[j] & 0x0f);
        if(Counter8 % 8 ==0)
          buf2idx = Int2CharArray(buf2,buf2idx,TmpIntVar);
        j = j + 1;
      }
    }
    if(SwTask[i].rightRlen > 0)
    {
      j = 0;
      while(j < SwTask[i].rightRlen)
      {
        Counter8 = Counter8 + 1;
        TmpIntVar = TmpIntVar <<4 | ((int)SwTask[i].rightRs[j] & 0x0f);
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
  short* output_ptr = new short[FPGA_RET_PARAM_NUM*BatchNum*2];

  memcpy(data_ptr,buf1,Buf1Len*sizeof(char) );
  memcpy(&data_ptr[Buf1Len/4],buf2,input_length-Buf1Len);                          
  delete buf1;
  delete buf2;             
  fprintf(stderr, "FPGA preparation used %dus\n", blaze::getUs()-start_ts); 

  start_ts = blaze::getUs();
  blaze::Task_ptr task = agent->createTask(acc_id);
  if (!task) {
    throw blaze::internalError("Task is not created");
  }
  agent->writeInput(task, acc_id, data_ptr, 1, input_length/4, sizeof(int));
  agent->writeInput(task, acc_id, &BatchNum, 1, 1, sizeof(int));

  agent->readOutput(task, output_ptr, FPGA_RET_PARAM_NUM*BatchNum*4);
  fprintf(stderr, "FPGA kernel used %dus\n", blaze::getUs()-start_ts); 

  start_ts = blaze::getUs();
  for (int i = 0; i < BatchNum; i++) {  
    
    int task_idx = ((int)(output_ptr[1+FPGA_RET_PARAM_NUM*2*i])<<16) |
                    output_ptr[0+FPGA_RET_PARAM_NUM*2*i];
    int regs_idx = SwTask[task_idx].idx;

    newregs[regs_idx].qb = output_ptr[2+FPGA_RET_PARAM_NUM*2*i]; 
    newregs[regs_idx].rb = output_ptr[4+FPGA_RET_PARAM_NUM*2*i] + SwTask[task_idx].rBeg;
    newregs[regs_idx].qe = output_ptr[3+FPGA_RET_PARAM_NUM*2*i] + SwTask[task_idx].qBeg + SwTask[task_idx].seedLength;
    newregs[regs_idx].re = output_ptr[5+FPGA_RET_PARAM_NUM*2*i] + SwTask[task_idx].rBeg + SwTask[task_idx].seedLength;
    newregs[regs_idx].score = output_ptr[6+FPGA_RET_PARAM_NUM*2*i]; 
    newregs[regs_idx].truesc = output_ptr[7+FPGA_RET_PARAM_NUM*2*i]; 
    newregs[regs_idx].w = output_ptr[8+FPGA_RET_PARAM_NUM*2*i];
  }
  fprintf(stderr, "FPGA output used %dus\n", blaze::getUs()-start_ts); 

  delete [] data_ptr;
  delete [] output_ptr;
}

#define MAX_BAND_TRY  2

void extendOnCPU(std::vector<ExtParam> tasks,int numoftask,mem_opt_t *opt,mem_alnreg_t *newregs) {
  for (int i=0;i<numoftask;i++){
    int aw[2];int max_off[2];
    int tmpidx = tasks[i].idx;
    aw[0]=opt->w;
    aw[1]=opt->w;
    if(tasks[i].qBeg){
      int qle,tle,gtle,gscore;
      for (int j=0;j<MAX_BAND_TRY;j++){
        int prev = newregs[tmpidx].score;
        aw[0] = opt->w<<j;
        newregs[tmpidx].score = ksw_extend2(tasks[i].leftQlen,tasks[i].leftQs,tasks[i].leftRlen,
        tasks[i].leftRs,5,opt->mat,tasks[i].oDel,tasks[i].eDel,tasks[i].oIns,tasks[i].eIns,
        aw[0],tasks[i].penClip5,tasks[i].zdrop,tasks[i].h0,&qle,&tle,&gtle,&gscore,&max_off[0]);
        if(newregs[tmpidx].score == prev||max_off[0]<(aw[0]>>1)+(aw[0]>>2)) break;
      }
      if (gscore <= 0 || gscore <= newregs[tmpidx].score - opt->pen_clip5) { // local extension
        newregs[tmpidx].qb = tasks[i].qBeg- qle;
        newregs[tmpidx].rb = tasks[i].rBeg- tle;
        newregs[tmpidx].truesc = newregs[tmpidx].score;
      }
      else { // to-end extension
        newregs[tmpidx].qb =0 ;
        newregs[tmpidx].rb =tasks[i].rBeg -gtle;
        newregs[tmpidx].truesc = gscore;
      }
    }
    else{
      newregs[tmpidx].score=newregs[tmpidx].truesc=tasks[i].h0;
      newregs[tmpidx].qb=0;
      newregs[tmpidx].rb=tasks[i].rBeg;
    }
    if(tasks[i].rightQlen){
      int qle,tle,gtle,gscore,sc0 = newregs[tmpidx].score;
      for (int j =0;j< MAX_BAND_TRY;++j){
        int prev = newregs[tmpidx].score;
        aw[1] = opt->w<<j;
        newregs[tmpidx].score= ksw_extend2(tasks[i].rightQlen,tasks[i].rightQs,tasks[i].rightRlen,
        tasks[i].rightRs,5,opt->mat,tasks[i].oDel,tasks[i].eDel,tasks[i].oIns,tasks[i].eIns,
        aw[0],tasks[i].penClip5,tasks[i].zdrop,sc0,&qle,&tle,&gtle,&gscore,&max_off[1]);
        if(newregs[tmpidx].score == prev||max_off[1]<(aw[1]>>1)+(aw[1]>>2)) break;
      }
      if (gscore <= 0 || gscore <= newregs[tmpidx].score - opt->pen_clip5) {
        newregs[tmpidx].qe = tasks[i].l_query-tasks[i].rightQlen + qle;
        newregs[tmpidx].re = tasks[i].rBeg + tasks[i].seedLength + tle;
        newregs[tmpidx].truesc +=newregs[tmpidx].score-sc0;
      }
      else{
        newregs[tmpidx].qe = tasks[i].l_query;
        newregs[tmpidx].re = tasks[i].rBeg + tasks[i].seedLength + gtle;
        newregs[tmpidx].truesc +=gscore-sc0;
      }
    }
    else{
      newregs[tmpidx].qe = tasks[i].l_query;
      newregs[tmpidx].re = tasks[i].rBeg + tasks[i].seedLength;
    }
    newregs[tmpidx].w = aw[0] > aw[1]? aw[0] : aw[1];
 }

}

static inline int cal_max_gap(const mem_opt_t *opt, int qlen)
{
	int l_del = (int)((double)(qlen * opt->a - opt->o_del) / opt->e_del + 1.);
	int l_ins = (int)((double)(qlen * opt->a - opt->o_ins) / opt->e_ins + 1.);
	int l = l_del > l_ins? l_del : l_ins;
	l = l > 1? l : 1;
	return l < opt->w<<1? l : opt->w<<1;
}

bool initCoordinates(int coordinates[][2],const MemChainVector *chains,int numofreads)
{ 
  bool isfinished = true;
  for (int i =0;i<numofreads;i++){
    coordinates[i][0]=0;
    if(&chains[i]!=NULL&&chains[i].n>0){
      coordinates[i][1]=chains[i].a[0].n-1;
       isfinished = false;
    }
    else coordinates[i][1]= -1;
  } 
  return isfinished;
}



increRes incrementCoordinates(int coordinates[][2],const MemChainVector *chains,int numofreads,int start,int end)
{
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

int testExtension(mem_opt_t *opt,mem_seed_t seed,mem_alnreg_v alnregv, int l_query)
{
  long rdist = -1;
  int qdist = -1;  
  int maxgap = -1;
  int mindist = -1;
  int w = -1;
  int breakidx = 0;
  bool isbreak = false;
  int i =0 ;
  for (i =0; i<alnregv.n;i++){
   if(seed.rbeg < alnregv.a[i].rb ||
      seed.rbeg + seed.len> alnregv.a[i].re||
      seed.qbeg < alnregv.a[i].qb ||
      seed.qbeg + seed.len >alnregv.a[i].qe )
      {
     continue; 
   }
   if(seed.len - alnregv.a[i].seedlen0 > .1 * l_query){
     continue;
  }
   qdist = seed.qbeg - alnregv.a[i].qb;
   rdist = seed.rbeg - alnregv.a[i].rb;
   mindist =(qdist<rdist)?qdist:(int)rdist;
   maxgap = cal_max_gap(opt,mindist);
   w = (maxgap < alnregv.a[i].w)? maxgap:alnregv.a[i].w ;
   if((qdist-rdist)<w &&(rdist-qdist)<w) {
     break;
   }
   qdist = alnregv.a[i].qe -(seed.qbeg + seed.len);
   rdist = alnregv.a[i].re - (seed.rbeg + seed.len);
   mindist =(qdist<rdist)?qdist:(int)rdist;
   maxgap = cal_max_gap(opt,mindist);
   w = (maxgap<alnregv.a[i].w)? maxgap:alnregv.a[i].w ;
   if((qdist-rdist)<w &&(rdist-qdist)<w) {
     break;
   }
  }
  return i ;
}

int checkoverlap(int startidx,mem_seed_t seed,mem_chain_t chain,uint64_t *srt)
{
  int i=startidx;
  for(i= startidx;i<chain.n;++i){
    const mem_seed_t *targetseed;
    if(srt[i]==0) continue;
    targetseed =&chain.seeds[(uint32_t)srt[i]];
    if(targetseed->len < seed.len* 0.95) continue;
    if(seed.qbeg <= targetseed->qbeg && 
       seed.qbeg + seed.len - targetseed->qbeg >= seed.len>>2 &&
       targetseed->qbeg-seed.qbeg != targetseed->rbeg-seed.rbeg) break;
    if(targetseed->qbeg <= seed.qbeg &&
       targetseed->qbeg + targetseed->len - seed.qbeg >= seed.len>>2 &&
       seed.qbeg-targetseed->qbeg != seed.rbeg-targetseed->rbeg) break;
  }
   return i ;
}

void mem_chain2aln_hw(
    ktp_aux_t *aux,
    const bseq1_t *seqs,
    const MemChainVector* chains,
    mem_alnreg_v *av,
    int batch_num)
{
  int numoftask = 0;
  int testCount = 0;
 // std::vector<std::vector<preResultofSw>> preResultofSw_m; 
  preResultofSw **preResultofSw_m = new preResultofSw*[batch_num];  // preresult matrix for each read
  int64_t l_pac = aux->idx->bns->l_pac;
  int rid, max_off[2], aw[2]; // aw: actual bandwidth used in extension
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
      rseq = bns_fetch_seq(aux->idx->bns, aux->idx->pac, &rmax[0], c->seeds[0].rbeg, &rmax[1], &rid);
      assert(c->rid == rid);

      srt = (uint64_t *)malloc(c->n * 8);
      for (int l = 0; l < c->n; ++l)
        srt[l] = (uint64_t)c->seeds[l].score<<32 | l;
      ks_introsort_64(c->n, srt);
      preResultofSw preresult(rmax,rseq,srt);
      preResultofSw_v[j]=preresult;
     // preResultofSw_v.push_back(preresult);    
    } 
  //  preResultofSw_m.push_back(preResultofSw_v);
    preResultofSw_m[i]= preResultofSw_v;
  } 
  int coordinates[batch_num][2];
  bool isfinished= initCoordinates(coordinates,chains,batch_num); 
  int extensionflags[batch_num];
  int overlapflags[batch_num];
  mem_seed_t **seedarray = new mem_seed_t*[batch_num];
  bool regflags[batch_num];

  int start = 0;
  int end = batch_num ;
  int taskidx = 0;
  mem_alnreg_t *newregs=new mem_alnreg_t[batch_num];
  newregs = (mem_alnreg_t*)malloc(batch_num*sizeof(mem_alnreg_t));
  for (int i=0;i<batch_num; i++){
    kv_init(av[i]);
  }
  while(!isfinished){
    taskidx = 0;
    int64_t pre_st= blaze::getUs();
    std::vector<ExtParam> sw_task_v;
    for (int i = start; i < end; i++) {
      regflags[i]= false;
      if(coordinates[i][1]>=0){
        seedarray[i]=& chains[i].
        a[coordinates[i][0]].
        seeds[(uint32_t)(preResultofSw_m[i][coordinates[i][0]].srt[coordinates[i][1]])];
        extensionflags[i]=testExtension(aux->opt,*seedarray[i],av[i],seqs[i].l_seq);
       // overlapflags[i]= -1;
        if(extensionflags[i]<av[i].n){
          overlapflags[i]= checkoverlap(coordinates[i][1]+1,*seedarray[i],chains[i].a[coordinates[i][0]],
              preResultofSw_m[i][coordinates[i][0]].srt);   
             if(overlapflags[i] == chains[i].a[coordinates[i][0]].n){
              preResultofSw_m[i][coordinates[i][0]].srt[coordinates[i][1]]= 0;
              continue;
          }
        }
          regflags[i]=true;
          newregs[i].score = seedarray[i]->len*aux->opt->a;
          newregs[i].truesc = seedarray[i]->len*aux->opt->a;
          newregs[i].qb = 0;
          newregs[i].rb = seedarray[i]->rbeg;
          newregs[i].qe = seqs[i].l_seq;
          newregs[i].re = seedarray[i]->rbeg + seedarray[i]->len; 
          newregs[i].rid = chains[i].a[coordinates[i][0]].rid;
          newregs[i].seedlen0 = seedarray[i]->len; 
          newregs[i].frac_rep = chains[i].a[coordinates[i][0]].frac_rep;
          newregs[i].w = aux->opt->w;
          GetTask(seedarray[i],aux->opt,&sw_task_v,(const uint8_t*)seqs[i].seq,
                  seqs[i].l_seq,preResultofSw_m[i][coordinates[i][0]].rmax[0] ,
                  preResultofSw_m[i][coordinates[i][0]].rmax[1],
                  preResultofSw_m[i][coordinates[i][0]].rseq,
                  i,&taskidx);
        
      }
    }
    int64_t cost_pre = blaze::getUs()-pre_st;
    //printf("prepare tasks used %dus for %d tasks\n",cost_pre,taskidx);
    ExtRet* SwResults = new ExtRet[taskidx]; 
    //ExtRet* SwResultsCPU = new ExtRet[taskidx]; 
    
    if (taskidx >= 500) {
      int64_t start_ts_hw = blaze::getUs();
      SwFPGA(sw_task_v,taskidx,newregs);
      int64_t cost_hw = blaze::getUs()-start_ts_hw;

#ifdef CHECK_FPGA
      int* fpga_qb     = new int[taskidx];
      int64_t* fpga_rb = new int64_t[taskidx];
      int* fpga_qe     = new int[taskidx];
      int64_t* fpga_re = new int64_t[taskidx];
      int* fpga_score  = new int[taskidx];
      int* fpga_truesc = new int[taskidx];
      int* fpga_w      = new int[taskidx];

      // dump FPGA results
      for (int i=0; i<taskidx; i++) {
        int regs_idx = sw_task_v[i].idx;
        fpga_qb[i]     = newregs[regs_idx].qb;
        fpga_rb[i]     = newregs[regs_idx].rb;
        fpga_qe[i]     = newregs[regs_idx].qe;
        fpga_re[i]     = newregs[regs_idx].re;
        fpga_score[i]  = newregs[regs_idx].score;
        fpga_truesc[i] = newregs[regs_idx].truesc;
        fpga_w[i]      = newregs[regs_idx].w;
      }
#endif

      int64_t start_ts_sw = blaze::getUs();
      extendOnCPU(sw_task_v,taskidx,aux->opt,newregs);
      int64_t cost_sw = blaze::getUs()-start_ts_sw;

#ifdef CHECK_FPGA
      // check with CPU results
      for (int i=0; i<taskidx; i++) {
        int regs_idx = sw_task_v[i].idx;
        if (fpga_qb[i] != newregs[regs_idx].qb) printf("#task=%d, #%d qb mismatch: %d!=%d\n", taskidx, regs_idx, fpga_qb[i], newregs[regs_idx].qb);
        if (fpga_rb[i] != newregs[regs_idx].rb) printf("#task=%d, #%d rb mismatch: %d!=%d\n", taskidx, regs_idx, fpga_rb[i], newregs[regs_idx].rb);
        if (fpga_qe[i] != newregs[regs_idx].qe) printf("#task=%d, #%d qe mismatch: %d!=%d\n", taskidx, regs_idx, fpga_qe[i], newregs[regs_idx].qe);
        if (fpga_re[i] != newregs[regs_idx].re) printf("#task=%d, #%d re mismatch: %d!=%d\n", taskidx, regs_idx, fpga_re[i], newregs[regs_idx].re);
        if (fpga_score[i] != newregs[regs_idx].score) printf("#task=%d, #%d score mismatch: %d!=%d\n", taskidx, regs_idx, fpga_score[i], newregs[regs_idx].score);
        if (fpga_truesc[i] != newregs[regs_idx].truesc) printf("#task=%d, #%d truesc mismatch: %d!=%d\n", taskidx, regs_idx, fpga_truesc[i], newregs[regs_idx].truesc);
        if (fpga_w[i] != newregs[regs_idx].w) printf("#task=%d, #%d w mismatch: %d!=%d\n", taskidx, regs_idx, fpga_w[i], newregs[regs_idx].w);
      }
      delete [] fpga_qb;
      delete [] fpga_rb;
      delete [] fpga_qe;
      delete [] fpga_re;
      delete [] fpga_score;
      delete [] fpga_truesc;
      delete [] fpga_w;
#endif

      fprintf(stderr, "hw used %dus and sw used %dus in %d tasks\n",
          cost_hw, cost_sw, taskidx); 
    }
    else{
      extendOnCPU(sw_task_v,taskidx,aux->opt,newregs);
    }
    for (int i = start; i < end; i++) {
      if(regflags[i]==true){
        newregs[i].seedcov=0;  // TODO:add the seedcov compute function
        for (int j=0; j<chains[i].a[coordinates[i][0]].n; ++j){
          const mem_seed_t *t =&chains[i].a[coordinates[i][0]].seeds[j];
          if(t->qbeg >= newregs[i].qb && 
             t->qbeg + t->len <= newregs[i].qe && 
             t->rbeg >= newregs[i].rb && 
             t->rbeg + t->len <= newregs[i].re){
             newregs[i].seedcov += t->len; 
          }
        }
        kv_push(mem_alnreg_t,av[i],newregs[i]);
      }
    }
    increRes nextiter = incrementCoordinates(coordinates,chains,batch_num,start,end);
    isfinished = nextiter.isfinished;
    start = nextiter.start;
    end = nextiter.end; 
  }
}
