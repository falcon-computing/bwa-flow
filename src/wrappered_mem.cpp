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

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

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

void GetTask(const mem_seed_t *seed, mem_opt_t *opt,ExtParam *SwTask , const uint8_t *query,int l_query,int64_t rmax_0, int64_t rmax_1,uint8_t *rseq,int idx)   //      the rmax_0 is rmax[0] and the rseq is the retrieved reference sequence
{
  if(seed->qbeg>0||(seed->qbeg+seed->len!=l_query))
  {
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
      for(int i = 0;i<=SwTask->rightQlen;i++)
        SwTask->rightQs[i] = query[i+qe];

      int re = seed->rbeg + seed->len - rmax_0;
      SwTask->rightRlen = rmax_1 - rmax_0 -re ;
      SwTask->rightRs = new uint8_t[SwTask->rightRlen];
      for(int i = 0;i<=SwTask->rightRlen; i++)
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

void SwFPGA(ExtParam *SwTask,int BatchNum,ExtRet *SwResult)
{
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
    buf1idx = Int2CharArray(buf1,buf1idx,TaskPos);
    i = i+1;
  }

  char *buf2 = new char[(TaskPos<<2)-Buf1Len];
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

  char *buf2fpga = new char[TaskPos<<2];
  strcat(buf2fpga,buf1);
  strcat(buf2fpga,buf2);                          // TODO:maybe there is another way to concat the two array, need to be changed
  delete buf1;
  delete buf2;
  // SwExtendFPGA(buf2fpga,BatchNum,SwResult);        //  This function send the packed data to fpga and get the result
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

void extendOnCPU(ExtParam *tasks,ExtRet *results,int numoftask,mem_opt_t *opt)
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
     // delete(tasks[i].leftQs);
     // delete(tasks[i].leftRs);
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
    // delete(tasks[i].rightQs);
    // delete(tasks[i].rightRs);
    }
    else{
      results[i].qEnd = 150;
      results[i].rEnd = tasks[i].rBeg + tasks[i].seedLength;
    }
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


void mem_chain2aln_hw(
    ktp_aux_t *aux,
    const bseq1_t *seqs,
    const MemChainVector* chains,
    mem_alnreg_v *av,
    int batch_num)
{
  int numoftask = 0;
  int testCount = 0;
  ExtParam *SwTask = new ExtParam[batch_num*10];                // maybe allocate the memory first or use vector
  memset(SwTask,0,batch_num*10*sizeof(ExtParam));
  ExtRet *SwResult = new ExtRet[batch_num*10];
  memset(SwResult,0,batch_num*10*sizeof(ExtRet));
  for(int i=0; i<batch_num; i++)    // loop for each seq
  {
    kv_init(av[i]);
    int z=0;
    for (int j=0; j<chains[i].n; j++) {    // loop for each chain
      mem_chain_t *c = &chains[i].a[j];
      // ----------------------------------prepare the maxspan and rseq for each seed
      int rid, max_off[2], aw[2]; // aw: actual bandwidth used in extension
      int64_t l_pac = aux->idx->bns->l_pac, rmax[2], tmp, max = 0;
      const mem_seed_t *s;
      uint8_t *rseq = 0;
      uint64_t *srt;
      if (c->n == 0) continue;
      // get the max possible span
      rmax[0] = l_pac<<1; rmax[1] = 0;
      for (int k = 0; k < c->n; ++k) {
        int64_t b, e;
        const mem_seed_t *t = &c->seeds[k];
        b = t->rbeg - (t->qbeg + cal_max_gap(aux->opt, t->qbeg));
        e = t->rbeg + t->len + ((seqs[i].l_seq - t->qbeg - t->len) + cal_max_gap(aux->opt, seqs[i].l_seq - t->qbeg - t->len));
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

      for (int m = c->n - 1; m >= 0; --m) {   // loop for each seed
        mem_alnreg_t *a;
        testCount +=1;
        s = &c->seeds[(uint32_t)srt[m]];
        // test extension here
        for (z=0;z<av[i].n;++z){
          mem_alnreg_t *p = &av[i].a[z];
          int64_t rd;
          int qd,w,max_gap;
          if (s->rbeg < p->rb || s->rbeg + s->len > p->re ||
          s->qbeg < p->qb || s->qbeg + s->len > p->qe) continue; // not fully contained
          if (s->len - p->seedlen0 > .1 * seqs[i].l_seq) continue;
          qd = s->qbeg - p->qb; rd = s->rbeg - p->rb;
          max_gap = cal_max_gap(aux->opt, qd < rd? qd : rd); // the maximal gap allowed in regions ahead of the seed
          w = max_gap < p->w? max_gap : p->w; // bounded by the band width
          if (qd - rd < w && rd - qd < w) break; // the seed is "around" a previous hit
          // similar to the previous four lines, but this time we look at the region behind
          qd = p->qe - (s->qbeg + s->len); rd = p->re - (s->rbeg + s->len);
          max_gap = cal_max_gap(aux->opt, qd < rd? qd : rd);
          w = max_gap < p->w? max_gap : p->w;
          if (qd - rd < w && rd - qd < w) break;
        }
        if (z<av[i].n){
          for (z = m+1;z<c->n;++z){
            const mem_seed_t *t;
            if (srt[z] == 0) continue;
            t = &c->seeds[(uint32_t)srt[z]];
            if (t->len < s->len * .95) continue; // only check overlapping if t is long enough; TODO: more efficient by early stopping
            if (s->qbeg <= t->qbeg && s->qbeg + s->len - t->qbeg
            >= s->len>>2 && t->qbeg - s->qbeg != t->rbeg - s->rbeg) break;
            if (t->qbeg <= s->qbeg && t->qbeg + t->len - s->qbeg
            >= s->len>>2 && s->qbeg - t->qbeg != s->rbeg - t->rbeg) break;
          }
          if (z == c->n) { // no overlapping seeds; then skip extension
            srt[m] = 0; // mark that seed extension has not been performed
            continue;
          }
        }
        a = kv_pushp(mem_alnreg_t, av[i]);
        memset(a, 0, sizeof(mem_alnreg_t));
        a->w = aw[0] = aw[1] = aux->opt->w;
        a->score = a->truesc = -1;
        a->rid = c->rid;
        GetTask(s,aux->opt,&SwTask[numoftask],(const uint8_t*)seqs[i].seq,seqs[i].l_seq,rmax[0],rmax[1],rseq,numoftask);
        numoftask += 1;
        //TODO: add the seedcov compute later
        /*a->seedcov = 0;
        for (int n=0;n <c->n; ++n){
          const mem_seed_t *t = &c->seeds[n];
          if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe && t->rbeg >= a->rb && t->rbeg + t->len <= a->re)
            a->seedcov += t->len; // this is not very accurate, but for approx. mapQ, this is good enough
        }
        a->w = aw[0]>aw[1]?aw[0]:aw[1];*/
        a->seedlen0 = s->len;
        a->frac_rep = c->frac_rep;
      }
      free(srt); free(rseq);
    }
  }
  // numoftask = 0;
  // send to fpga and get the results, for now use cpu
  extendOnCPU(SwTask,SwResult,numoftask,aux->opt);
  //---------------------------------
  numoftask = 0;
  for (int i=0; i<batch_num; i++)
    {
      for (int j=0; j<chains[i].n; j++){
      mem_chain_t *c = &chains[i].a[j];
        for (int m=chains[i].a[j].n - 1; m>=0; --m){
            av[i].a[j].score = SwResult[numoftask].score;
            av[i].a[j].qb = SwResult[numoftask].qBeg;
            av[i].a[j].rb = SwResult[numoftask].rBeg;
            av[i].a[j].qe = SwResult[numoftask].qEnd;
            av[i].a[j].re = SwResult[numoftask].rEnd;
            av[i].a[j].truesc = SwResult[numoftask].trueScore;
            numoftask += 1;
            for (int n=0;n <c->n; ++n){
              const mem_seed_t *t = &c->seeds[n];
              if (t->qbeg >= av[i].a[j].qb && t->qbeg + t->len <= av[i].a[j].qe && t->rbeg >= av[i].a[j].rb && t->rbeg + t->len <= av[i].a[j].re)
                av[i].a[j].seedcov += t->len; // this is not very accurate, but for approx. mapQ, this is good enough
            }
        }
      }
    }
}



