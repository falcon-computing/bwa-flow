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
#include "bwamem.h"
#include "bntseq.h"
#include "ksw.h"
#include "kvec.h"
#include "ksort.h"
#include "utils.h"
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

void GetTask(mem_seed_t *seed, mem_opt_t *opt,ExtParam *SwTask , const uint8_t *query,int l_query,int64_t rmax_0, int64_t rmax_1,uint8_t *rseq,int idx)   //      the rmax_0 is rmax[0] and the rseq is the retrieved reference sequence
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
    SwTask->idx = idx ;
  }
}

int Int2CharArray(char* arr, int idx, int num)
{
  arr[idx] = (char)(num&0xff);
  arr[idx+1] = (char)((num>>8)&0xff);
  arr[idx+2] = (char)((num>>16)&0xff);
  arr[idx+3] = (char)((num>>24)&0xff);
  return idx+4 ;
}

int Short2CharArray(char *arr, int idx, short num)
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

void mem_chain2aln_hw(
    ktp_aux_t *aux,
    const bseq1_t *seqs, 
    const MemChainVector* chains, 
    mem_alnreg_v *av,
    int batch_num)
{
  for(int i=0; i<batch_num; i++)
  {
    kv_init(av[i]);
    for (int j=0; j<chains[i].n; j++) {
      mem_chain_t *p = &chains[i].a[j];

      // call mem_chain2aln to compute baseline
      mem_chain2aln(aux->opt, aux->idx->bns, aux->idx->pac, 
          seqs[i].l_seq, (uint8_t*)seqs[i].seq, p, av+i);
    }
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
