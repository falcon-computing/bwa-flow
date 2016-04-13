#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

//#include "kstring.h"
#include "bwa/bwamem.h"
#include "bwa/bntseq.h"
#include "bwa/ksw.h"
#include "bwa/kvec.h"
#include "bwa/ksort.h"
#include "bwa/utils.h"
#include "bwa_wrapper.h"
#include "blaze/AccAgent.h"
#include "SWClient.h"

#define FPGA_RET_PARAM_NUM 5


#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif
#include <vector>
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
      SwTask->leftRlen = seed->rbeg - rmax_0 ;
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
      int re = seed->rbeg + seed->len - rmax_0;
      SwTask->rightRlen = rmax_1 - rmax_0 -re ;
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

void SwFPGA(std::vector<ExtParam> SwTask,int BatchNum,ExtRet *SwResult)
{
  SWClient client;
  
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
    LeftMaxIns = (int)((double)(SwTask[i].leftQlen*1+ SwTask[i].penClip5 -SwTask[i].oIns)/SwTask[i].eIns+1);
    LeftMaxDel = (int)((double)(SwTask[i].leftQlen*1 + SwTask[i].penClip5 -SwTask[i].oDel)/SwTask[i].eDel+1);
    RightMaxIns = (int)((double)(SwTask[i].rightQlen*1 + SwTask[i].penClip3 -SwTask[i].oIns)/SwTask[i].eIns+1);
    RightMaxIns = (int)((double)(SwTask[i].rightQlen*1 + SwTask[i].penClip3 -SwTask[i].oDel)/SwTask[i].eDel+1);         // 1 stands for SwTask[i].mat.max
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(LeftMaxIns));
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(LeftMaxDel));
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(RightMaxIns));
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(RightMaxDel));    // dont know why dont change to short types in the upper lines, but the scala codes did that
    buf1idx = Int2CharArray(buf1,buf1idx,SwTask[i].idx);
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

   int *data_ptr = (int*)client.createInput(0,
                         1,input_length/4,sizeof(int),BLAZE_INPUT);
  int* taskNum_ptr = (int*)client.createInput(1,
                         1,1,sizeof(int),BLAZE_INPUT);
  *taskNum_ptr = BatchNum;
  memcpy(data_ptr,buf1,Buf1Len*sizeof(char) );
  memcpy(&data_ptr[Buf1Len/4],buf2,input_length-Buf1Len);                          
 // FILE *fout = fopen("dump_yh.dat","wb");
 // fwrite(&BatchNum,1,sizeof(int),fout);
 // fwrite(&input_length,1,sizeof(int),fout);
 // fwrite(data_ptr,input_length/4,sizeof(int),fout);
 // fclose(fout);
  delete buf1;
  delete buf2;             
  client.start();          
  short* output_ptr = (short*)client.getOutputPtr(0);
 // ExtRet* results =(ExtRet*) malloc(BatchNum*sizeof(ExtRet));
  
  for (int i=0;i<BatchNum;i++){
     SwResult[i].idx=(int)(output_ptr[1+FPGA_RET_PARAM_NUM*2*i])<<16|
     (int)output_ptr[0+FPGA_RET_PARAM_NUM*2*i];
     SwResult[i].qBeg = output_ptr[2+FPGA_RET_PARAM_NUM*2*i];
     SwResult[i].qEnd = output_ptr[3+FPGA_RET_PARAM_NUM*2*i];
     SwResult[i].rBeg = output_ptr[4+FPGA_RET_PARAM_NUM*2*i];
     SwResult[i].rEnd = output_ptr[5+FPGA_RET_PARAM_NUM*2*i];
     SwResult[i].score= output_ptr[6+FPGA_RET_PARAM_NUM*2*i];
     SwResult[i].trueScore= output_ptr[7+FPGA_RET_PARAM_NUM*2*i];
     SwResult[i].width= output_ptr[8+FPGA_RET_PARAM_NUM*2*i];
   }
}
void SwFPGA_old(ExtParam *SwTask,int BatchNum,ExtRet *SwResult)
{
  SWClient client;
  int Buf1Len = 32 + 32*BatchNum;
  char* buf1 = new char[Buf1Len];
  //------------------ store the public options at the beginning-----------------
  buf1[0] = (char)SwTask->oDel;
  buf1[1] = (char)SwTask->eDel;
  buf1[2] = (char)SwTask->oIns;
  buf1[3] = (char)SwTask->eIns;
  buf1[4] = (char)SwTask->penClip5;
  buf1[5] = (char)SwTask->penClip3;
  buf1[6] = (char)SwTask->w;
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
    LeftMaxIns = (int)((double)(SwTask[i].leftQlen*1+ SwTask[i].penClip5 -SwTask[i].oIns)/SwTask[i].eIns+1);
    LeftMaxDel = (int)((double)(SwTask[i].leftQlen*1 + SwTask[i].penClip5 -SwTask[i].oDel)/SwTask[i].eDel+1);
    RightMaxIns = (int)((double)(SwTask[i].rightQlen*1 + SwTask[i].penClip3 -SwTask[i].oIns)/SwTask[i].eIns+1);
    RightMaxIns = (int)((double)(SwTask[i].rightQlen*1 + SwTask[i].penClip3 -SwTask[i].oDel)/SwTask[i].eDel+1);         // 1 stands for SwTask[i].mat.max
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(LeftMaxIns));
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(LeftMaxDel));
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(RightMaxIns));
    buf1idx = Short2CharArray(buf1,buf1idx,(short)(RightMaxDel));    // dont know why dont change to short types in the upper lines, but the scala codes did that
    buf1idx = Int2CharArray(buf1,buf1idx,SwTask[i].idx);
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

   int *data_ptr = (int*)client.createInput(0,
                         1,input_length/4,sizeof(int),BLAZE_INPUT);
  int* taskNum_ptr = (int*)client.createInput(1,
                         1,1,sizeof(int),BLAZE_INPUT);
  *taskNum_ptr = BatchNum;
  memcpy(data_ptr,buf1,Buf1Len*sizeof(char) );
  memcpy(&data_ptr[Buf1Len/4],buf2,input_length-Buf1Len);                          // TODO:maybe there is another way to concat the two array, need to be changed
  FILE *fout = fopen("dump_yh.dat","wb");
  fwrite(&BatchNum,1,sizeof(int),fout);
  fwrite(&input_length,1,sizeof(int),fout);
  fwrite(data_ptr,input_length/4,sizeof(int),fout);
  fclose(fout);
  delete buf1;
  delete buf2;             
 // delete buf2fpga;
  client.start();          
  short* output_ptr = (short*)client.getOutputPtr(0);
  ExtRet* results =(ExtRet*) malloc(BatchNum*sizeof(ExtRet));
  
  for (int i=0;i<BatchNum;i++){
     results[i].idx=(int)(output_ptr[1+FPGA_RET_PARAM_NUM*2*i])<<16|
     (int)output_ptr[0+FPGA_RET_PARAM_NUM*2*i];
     results[i].qBeg = output_ptr[2+FPGA_RET_PARAM_NUM*2*i];
     results[i].qEnd = output_ptr[3+FPGA_RET_PARAM_NUM*2*i];
     results[i].rBeg = output_ptr[4+FPGA_RET_PARAM_NUM*2*i];
     results[i].rEnd = output_ptr[5+FPGA_RET_PARAM_NUM*2*i];
     results[i].score= output_ptr[6+FPGA_RET_PARAM_NUM*2*i];
     results[i].trueScore= output_ptr[7+FPGA_RET_PARAM_NUM*2*i];
     results[i].width= output_ptr[8+FPGA_RET_PARAM_NUM*2*i];
  
   }
  for (int i = 0;i<BatchNum;i++){
     SwResult[results[i].idx].idx=results[i].idx;
     SwResult[results[i].idx].qBeg=results[i].qBeg ;
     SwResult[results[i].idx].qEnd=results[i].qEnd+ SwTask[i].qBeg + SwTask[i].seedLength;
     SwResult[results[i].idx].rBeg=results[i].rBeg+ SwTask[i].rBeg ;
     SwResult[results[i].idx].rEnd=results[i].rEnd+ SwTask[i].rBeg + SwTask[i].seedLength;
     SwResult[results[i].idx].score=results[i].score;
     SwResult[results[i].idx].trueScore=results[i].trueScore;
     SwResult[results[i].idx].width=results[i].width;
  }                  
  free(results);
}



void seq2intv(ktp_aux_t *aux,bseq1_t *seqs,smem_aux_t *SMEM)
{
  int i ;
  for (i = 0; i < seqs->l_seq; ++i) // convert to 2-bit encoding if we have not done so
    seqs->seq[i] = seqs->seq[i] < 4? seqs->seq[i] : nst_nt4_table[(int)seqs->seq[i]];
  mem_collect_intv(aux->opt,aux->idx->bwt,seqs->l_seq,(uint8_t *)seqs->seq,SMEM);
  SMEM->id_read = seqs->id;
}


MemChainVector seq2chain(ktp_aux_t *aux,bseq1_t *seqs)
{
  int i;
  MemChainVector chain_wid;
  mem_chain_v chn;
  for (i = 0; i < seqs->l_seq; ++i) // convert to 2-bit encoding if we have not done so
    seqs->seq[i] = seqs->seq[i] < 4? seqs->seq[i] : nst_nt4_table[(int)seqs->seq[i]];

  chn = mem_chain(aux->opt, aux->idx->bwt, aux->idx->bns, seqs->l_seq, (uint8_t*)seqs->seq, 0);         // the 0 should be reconsidered
  chn.n = mem_chain_flt(aux->opt, chn.n, chn.a);
  mem_flt_chained_seeds(aux->opt, aux->idx->bns, aux->idx->pac, seqs->l_seq, (uint8_t*)seqs->seq, chn.n, chn.a);
  chain_wid.id_read = seqs->id;
  chain_wid.n = chn.n;
  chain_wid.m = chn.m;
  chain_wid.a = chn.a;
  return chain_wid;
}


void chain2reg(ktp_aux_t *aux,bseq1_t *seqs,MemChainVector chn,mem_alnreg_v *alnreg)
{
  int id = chn.id_read;
  int i;
  kv_init(*alnreg);
  for (i = 0; i < chn.n; ++i) {
    mem_chain_t *p = &chn.a[i];
    if (bwa_verbose >= 4) err_printf("* ---> Processing chain(%d) <---\n", i);
    mem_chain2aln(aux->opt, aux->idx->bns, aux->idx->pac, seqs[id].l_seq, (uint8_t*)seqs[id].seq, p, alnreg);
    free(chn.a[i].seeds);
  }
  free(chn.a);
  alnreg->n = mem_sort_dedup_patch(aux->opt, aux->idx->bns, aux->idx->pac, (uint8_t*)seqs[id].seq, alnreg->n, alnreg->a);
  if (bwa_verbose >= 4) {
    err_printf("* %ld chains remain after removing duplicated chains\n", alnreg->n);
    for (i = 0; i < alnreg->n; ++i) {
      mem_alnreg_t *p = &alnreg->a[i];
      printf("** %d, [%d,%d) <=> [%ld,%ld)\n", p->score, p->qb, p->qe, (long)p->rb, (long)p->re);
    }
  }
  for (i = 0; i < alnreg->n; ++i) {
    mem_alnreg_t *p = &alnreg->a[i];
    if (p->rid >= 0 && aux->idx->bns->anns[p->rid].is_alt)
      p->is_alt = 1;
  }
}


void reg2sam(ktp_aux_t *aux,bseq1_t *seqs,int batch_num,int64_t n_processed,mem_alnreg_v *alnreg)
{
  int i ;
  mem_pestat_t pes[4];
  mem_pestat(aux->opt, aux->idx->bns->l_pac, batch_num, alnreg, pes);               // the operation between worker1 and worker2
  for (i=0; i<batch_num/2; i++)
  {
    mem_sam_pe(aux->opt,aux->idx->bns,aux->idx->pac,pes,n_processed+i,&seqs[i<<1],&alnreg[i<<1]);
    free(alnreg[i<<1|0].a); free(alnreg[i<<1|1].a);
  }
  free(alnreg);
  // print
  for (i = 0; i < batch_num; ++i)
  {
    if (seqs[i].sam) err_fputs(seqs[i].sam, stdout);
    free(seqs[i].name); free(seqs[i].comment);
    free(seqs[i].seq); free(seqs[i].qual); free(seqs[i].sam);
  }
  free(seqs);
}

void reg_dump(mem_alnreg_v *alnreg,mem_alnreg_v *alnreg_hw,int batch_num)
{
  // print the software result
  FILE *fp_old = xopen("reg_sw.txt","wb");
  FILE *fp_new = xopen("reg_hw.txt","wb");
  int i = 0;
  int j = 0;
  for(j=0;j<batch_num;j++)
  {
    for(i=0;i<alnreg->n;i++)
    {
      err_fwrite(&alnreg->a[i].rb,sizeof(int64_t),1,fp_old);
      err_fwrite(&alnreg->a[i].re,sizeof(int64_t),1,fp_old);
      err_fwrite(&alnreg->a[i].qb,sizeof(int),1,fp_old);
      err_fwrite(&alnreg->a[i].qe,sizeof(int),1,fp_old);
      err_fwrite(&alnreg->a[i].rid,sizeof(int),1,fp_old);
      err_fwrite(&alnreg->a[i].score,sizeof(int),1,fp_old);
      err_fwrite(&alnreg->a[i].sub,sizeof(int),1,fp_old);
      err_fwrite(&alnreg->a[i].alt_sc,sizeof(int),1,fp_old);
      err_fwrite(&alnreg->a[i].csub,sizeof(int),1,fp_old);
      err_fwrite(&alnreg->a[i].sub_n,sizeof(int),1,fp_old);
      err_fwrite(&alnreg->a[i].w,sizeof(int),1,fp_old);
      err_fwrite(&alnreg->a[i].seedcov,sizeof(int),1,fp_old);
      err_fwrite(&alnreg->a[i].secondary,sizeof(int),1,fp_old);
      err_fwrite(&alnreg->a[i].secondary_all,sizeof(int),1,fp_old);
      err_fwrite(&alnreg->a[i].seedlen0,sizeof(int),1,fp_old);
      //         		err_fwrite(&alnreg->a[i].n_comp,sizeof(int),1,fp_old);
      //         		err_fwrite(&alnreg->a[i].is_alt,sizeof(int),1,fp_old);
      err_fwrite(&alnreg->a[i].frac_rep,sizeof(int),1,fp_old);
      err_fwrite(&alnreg->a[i].hash,sizeof(uint64_t),1,fp_old);
      err_fwrite("\r\n",1,2,fp_old);
    }
    // print the hardware result
    for(i=0;i<alnreg_hw->n;i++)
    {
      err_fwrite(&alnreg_hw->a[i].rb,sizeof(int64_t),1,fp_new);
      err_fwrite(&alnreg_hw->a[i].re,sizeof(int64_t),1,fp_new);
      err_fwrite(&alnreg_hw->a[i].qb,sizeof(int),1,fp_new);
      err_fwrite(&alnreg_hw->a[i].qe,sizeof(int),1,fp_new);
      err_fwrite(&alnreg_hw->a[i].rid,sizeof(int),1,fp_new);
      err_fwrite(&alnreg_hw->a[i].score,sizeof(int),1,fp_new);
      err_fwrite(&alnreg_hw->a[i].sub,sizeof(int),1,fp_new);
      err_fwrite(&alnreg_hw->a[i].alt_sc,sizeof(int),1,fp_new);
      err_fwrite(&alnreg_hw->a[i].csub,sizeof(int),1,fp_new);
      err_fwrite(&alnreg_hw->a[i].sub_n,sizeof(int),1,fp_new);
      err_fwrite(&alnreg_hw->a[i].w,sizeof(int),1,fp_new);
      err_fwrite(&alnreg_hw->a[i].seedcov,sizeof(int),1,fp_new);
      err_fwrite(&alnreg_hw->a[i].secondary,sizeof(int),1,fp_new);
      err_fwrite(&alnreg_hw->a[i].secondary_all,sizeof(int),1,fp_new);
      err_fwrite(&alnreg_hw->a[i].seedlen0,sizeof(int),1,fp_new);
      //              err_fwrite(&alnreg_hw->a[i].n_comp,sizeof(int),1,fp_new);
      //        		err_fwrite(&alnreg_hw->a[i].is_alt,sizeof(int),1,fp_new);
      err_fwrite(&alnreg_hw->a[i].frac_rep,sizeof(int),1,fp_new);
      err_fwrite(&alnreg_hw->a[i].hash,sizeof(uint64_t),1,fp_new);
      err_fwrite("\r\n",1,2,fp_new);
    }
  }
  err_fflush(fp_old);
  err_fclose(fp_old);
  err_fflush(fp_new);
  err_fclose(fp_new);
}


#define MAX_BAND_TRY  2

void extendOnCPU(std::vector<ExtParam> tasks,ExtRet *results,int numoftask,mem_opt_t *opt)
{
  for (int i=0;i<numoftask;i++){
    int aw[2];int max_off[2];
    if(tasks[i].qBeg){
      int qle,tle,gtle,gscore;
      for (int j=0;j<MAX_BAND_TRY;j++){
        int prev = results[i].score;
        aw[0] = opt->w<<j;
        results[i].score = ksw_extend2(tasks[i].leftQlen,tasks[i].leftQs,tasks[i].leftRlen,
        tasks[i].leftRs,5,opt->mat,tasks[i].oDel,tasks[i].eDel,tasks[i].oIns,tasks[i].eIns,
        aw[0],tasks[i].penClip5,tasks[i].zdrop,tasks[i].h0,&qle,&tle,&gtle,&gscore,&max_off[0]);
        if(results[i].score == prev||max_off[0]<(aw[0]>>1)+(aw[0]>>2)) break;
      }
      if (gscore <= 0 || gscore <= results[i].score - opt->pen_clip5) { // local extension
        results[i].qBeg = tasks[i].qBeg- qle;
        results[i].rBeg = tasks[i].rBeg- tle;
        results[i].trueScore = results[i].score;
      }
      else { // to-end extension
        results[i].qBeg =0 ;
        results[i].rBeg =tasks[i].rBeg -gtle;
        results[i].trueScore = gscore;
      }
    }
    else{
      results[i].score=results[i].trueScore=tasks[i].h0;
      results[i].qBeg=0;
      results[i].rBeg=tasks[i].rBeg;
    }
    if(tasks[i].rightQlen){
      int qle,tle,gtle,gscore,sc0 = results[i].score;
      for (int j =0;j< MAX_BAND_TRY;++j){
        int prev = results[i].score;
        aw[1] = opt->w<<j;
        results[i].score= ksw_extend2(tasks[i].rightQlen,tasks[i].rightQs,tasks[i].rightRlen,
        tasks[i].rightRs,5,opt->mat,tasks[i].oDel,tasks[i].eDel,tasks[i].oIns,tasks[i].eIns,
        aw[0],tasks[i].penClip5,tasks[i].zdrop,sc0,&qle,&tle,&gtle,&gscore,&max_off[1]);
        if(results[i].score == prev||max_off[1]<(aw[1]>>1)+(aw[1]>>2)) break;
      }
      if (gscore <= 0 || gscore <= results[i].score - opt->pen_clip5) {
        results[i].qEnd = 150-tasks[i].rightQlen + qle;
        results[i].rEnd = tasks[i].rBeg + tasks[i].seedLength + tle;
        results[i].trueScore +=results[i].score-sc0;
      }
      else{
        results[i].qEnd = 150;
        results[i].rEnd = tasks[i].rBeg + tasks[i].seedLength + gtle;
        results[i].trueScore +=gscore-sc0;
      }
    }
    else{
      results[i].qEnd = 150;
      results[i].rEnd = tasks[i].rBeg + tasks[i].seedLength;
    }
    results[i].width = aw[0] > aw[1]? aw[0] : aw[1];
    results[i].idx= tasks[i].idx;
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


class preResultofSw{
public:
   int64_t *rmax;
   uint8_t *rseq;
   uint64_t *srt;
   preResultofSw(int64_t *_rmax,uint8_t *_rseq,uint64_t *_srt)
    :rseq(_rseq),srt(_srt)
   {
     rmax = new int64_t[2];
     rmax[0]= _rmax[0];
     rmax[1]=_rmax[1];
   };
};


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

class increRes{
public:
  bool isfinished;
  int start; 
  int end;
  increRes(bool _isfinished,int _start,int _end): 
    isfinished(_isfinished),start(_start),end(_end)
    {};
};

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

int testExtension(mem_opt_t *opt,mem_seed_t seed,mem_alnreg_v alnregv)
{
  long rdist = -1;
  int qdist = -1;  
  int maxgap = -1;
  int mindist = -1;
  int w = -1;
  int breakidx = 0;
  bool isbreak = false;
  int i = 0;
  while(i<alnregv.n&&!isbreak){
    if(seed.rbeg>=alnregv.a[i].rb&&(seed.rbeg+seed.len)<=alnregv.a[i].re&&
      seed.qbeg>=alnregv.a[i].qb&&(seed.qbeg+seed.len)<=alnregv.a[i].qe){
      qdist = seed.qbeg - alnregv.a[i].qb;
      rdist = seed.rbeg - alnregv.a[i].rb;
      mindist =(qdist<rdist)?qdist:(int)rdist;
      maxgap = cal_max_gap(opt,mindist);
      w = (maxgap<opt->w)? maxgap:opt->w ;
      if((qdist-rdist)<w &&(rdist-qdist)<w) {
        breakidx = i;
        isbreak = true;
      }
      if(!isbreak){
        qdist = alnregv.a[i].qe -(seed.qbeg + seed.len);
        rdist = alnregv.a[i].re - (seed.rbeg + seed.len);
        mindist =(qdist<rdist)?qdist:(int)rdist;
        maxgap = cal_max_gap(opt,mindist);
        w = (maxgap<opt->w)? maxgap:opt->w ;
        if((qdist-rdist)<w &&(rdist-qdist)<w) {
          breakidx = i;
          isbreak = true;
        } 
         
      }
    }
    i+=1;
  }
  if(isbreak) return breakidx;
  else return i ; 
}

int checkoverlap(int startidx,mem_seed_t seed,mem_chain_t chain,uint64_t *srt)
{
  int breakidx = chain.n;
  int i = startidx;
  bool isbreak = false;
  while(i<chain.n &&!isbreak){
    if(srt[i]!=0){
      mem_seed_t targetseed= chain.seeds[(uint32_t)srt[i]];
      if(targetseed.len >=seed.len*0.95){
        if(seed.qbeg<=targetseed.qbeg&&(seed.qbeg+seed.len-targetseed.qbeg)>=(seed.len>>2)&&
           (targetseed.qbeg-seed.qbeg)!=(targetseed.qbeg-seed.rbeg)){
          breakidx = i;
          isbreak = true;
        }
        if(!isbreak&&targetseed.qbeg<=seed.qbeg&&( targetseed.qbeg+targetseed.len-seed.qbeg)>=
            (seed.len>>2)&&(seed.qbeg-targetseed.qbeg)!=(seed.rbeg-targetseed.rbeg)){
          breakidx = i;
          isbreak = true;
        }
      }   
    }
    i +=1;
  }
  if(isbreak) return breakidx;
  else return i;
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
  std::vector<ExtParam> sw_task_v;
  std::vector<std::vector<preResultofSw>> preResultofSw_m; 
  int64_t l_pac = aux->idx->bns->l_pac;
  int rid, max_off[2], aw[2]; // aw: actual bandwidth used in extension
  for (int i=0; i<batch_num; i++)   //loop for each seq
  {
    std::vector<preResultofSw> preResultofSw_v;
    for (int j=0;j<chains[i].n;j++){
      mem_chain_t *c = &chains[i].a[j]; // ------------prepare the maxspan and rseq for each seed
      int64_t rmax[2], tmp, max = 0;
      uint8_t *rseq = 0;
      uint64_t *srt;
      if (c->n == 0) continue;
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
      preResultofSw_v.push_back(preresult);    
    } 
    preResultofSw_m.push_back(preResultofSw_v);
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
    int i =start;
    while (i <end){
      regflags[i]= false;
      if(coordinates[i][1]>=0){
        seedarray[i]=& chains[i].
        a[coordinates[i][0]].
        seeds[(uint32_t)(preResultofSw_m[i][coordinates[i][0]].srt[coordinates[i][1]])];
        extensionflags[i]=testExtension(aux->opt,*seedarray[i],av[i]);
        overlapflags[i]= -1;
        if(extensionflags[i]<av[i].n)
          overlapflags[i]= checkoverlap(coordinates[i][1]+1,*seedarray[i],chains[i].a[coordinates[i][0]],
              preResultofSw_m[i][coordinates[i][0]].srt);   
        if(extensionflags[i]<av[i].n&&overlapflags[i]==chains[i].a[coordinates[i][0]].n)
              preResultofSw_m[i][coordinates[i][0]].srt[coordinates[i][1]]= 0;
        else{
          regflags[i]=true;
          newregs[i].score = seedarray[i]->len*aux->opt->a;
          newregs[i].truesc = seedarray[i]->len*aux->opt->a;
          newregs[i].qb = 0;
          newregs[i].rb = seedarray[i]->rbeg;
          newregs[i].qe = seqs[i].l_seq;
          newregs[i].re = seedarray[i]->rbeg + seedarray[i]->len; 
          newregs[i].rid = chains[i].a[coordinates[i][0]].rid;
          newregs[i].seedlen0 = seedarray[i]->len; 
         /* ExtParam sw_task_temp;
          GetTask(seedarray[i],aux->opt,&sw_task_temp,(const uint8_t*)seqs[i].seq,
                  seqs[i].l_seq,preResultofSw_m[i][coordinates[i][0]].rmax[0] ,
                  preResultofSw_m[i][coordinates[i][0]].rmax[1],
                  preResultofSw_m[i][coordinates[i][0]].rseq,
                  i);
          sw_task_v.push_back(sw_task_temp);
          taskidx += 1;*/                                  //the old GetTask function
          GetTask(seedarray[i],aux->opt,&sw_task_v,(const uint8_t*)seqs[i].seq,
                  seqs[i].l_seq,preResultofSw_m[i][coordinates[i][0]].rmax[0] ,
                  preResultofSw_m[i][coordinates[i][0]].rmax[1],
                  preResultofSw_m[i][coordinates[i][0]].rseq,
                  i,&taskidx);
        }
      }
      i = i+1;
    }
    int64_t cost_pre = blaze::getUs()-pre_st;
    printf("prepare tasks used %dus\n",cost_pre);
    ExtRet* SwResults = new ExtRet[taskidx]; 
    ExtRet* SwResultsCPU = new ExtRet[taskidx]; 
    if(taskidx>=20){                                       // start the sw compute both on fpga and cpu   
      int64_t start_ts_hw = blaze::getUs();
      SwFPGA(sw_task_v,taskidx,SwResults);
      int64_t cost_hw = blaze::getUs()-start_ts_hw;
      int64_t start_ts_sw = blaze::getUs();
      extendOnCPU(sw_task_v,SwResultsCPU,taskidx,aux->opt);
      int64_t cost_sw = blaze::getUs()-start_ts_sw;
      printf("hw used %dus and sw used %dus in %d tasks\n",cost_hw,cost_sw,taskidx); 
    }
    else{
      extendOnCPU(sw_task_v,SwResults,taskidx,aux->opt);
    }
    i = 0;  
    while(i<taskidx){
      int tmpidx = SwResults[i].idx;
      newregs[tmpidx].qb = SwResults[i].qBeg;
      newregs[tmpidx].rb = SwResults[i].rBeg + seedarray[tmpidx]->rbeg;
      newregs[tmpidx].qe = SwResults[i].qEnd + seedarray[tmpidx]->qbeg +seedarray[tmpidx]->len;
      newregs[tmpidx].re = SwResults[i].rEnd + seedarray[tmpidx]->rbeg +seedarray[tmpidx]->len; 
      newregs[tmpidx].score = SwResults[i].score; 
      newregs[tmpidx].truesc = SwResults[i].trueScore; 
      newregs[tmpidx].w = SwResults[i].width;
      i = i+1;
    }
    i = start;
    while(i<end){
      if(regflags[i]==true){
        newregs[i].seedcov=0;  // TODO:add the seedcov compute function
        kv_push(mem_alnreg_t,av[i],newregs[i]);
      }
      i = i+1;
    }
    increRes nextiter = incrementCoordinates(coordinates,chains,batch_num,start,end);
    isfinished = nextiter.isfinished;
    start = nextiter.start;
    end = nextiter.end; 
  }
}



