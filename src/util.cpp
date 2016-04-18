#include <stdio.h>
#include <stdlib.h>

#include "bwa/bwamem.h"
#include "bwa/bntseq.h"
#include "bwa/utils.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#include "util.h"

template <typename T>
void assertEQ(T const &a, T const &b) {
  if (a != b) {
    throw resultsError("");
  }
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
      err_fwrite(&alnreg->a[i].frac_rep,sizeof(int),1,fp_old);
      err_fwrite(&alnreg->a[i].hash,sizeof(uint64_t),1,fp_old);
      err_fwrite("\r\n",1,2,fp_old);
    }
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

void regionsCompare(
    mem_alnreg_v *alnreg_base,
    mem_alnreg_v *alnreg_test,
    int num_seqs)
{
  try {
    for (int i = 0; i < num_seqs; i++) {
      // Check if the number of aligned regions match
      assertEQ(alnreg_base[i].n, alnreg_test[i].n);
      for (int j = 0; j < alnreg_base[i].n; j++) {
        mem_alnreg_t* region_base = &alnreg_base[i].a[j];
        mem_alnreg_t* region_test = &alnreg_test[i].a[j];

        assertEQ(region_base->rb, region_test->rb);
        assertEQ(region_base->re, region_test->re);
        assertEQ(region_base->qb, region_test->qb);
        assertEQ(region_base->qe, region_test->qe);
        assertEQ(region_base->rid, region_test->rid);
        assertEQ(region_base->score, region_test->score);
        assertEQ(region_base->sub, region_test->sub);
        assertEQ(region_base->alt_sc, region_test->alt_sc);
        assertEQ(region_base->csub, region_test->csub);
        assertEQ(region_base->sub_n, region_test->sub_n);
        assertEQ(region_base->w, region_test->w);
        assertEQ(region_base->seedcov, region_test->seedcov);
        assertEQ(region_base->secondary, region_test->secondary);
        assertEQ(region_base->secondary_all, region_test->secondary_all);
        assertEQ(region_base->seedlen0, region_test->seedlen0);
        assertEQ(region_base->frac_rep, region_test->frac_rep);
        assertEQ(region_base->hash, region_base->hash);
      }
    }
    fprintf(stderr, "Correct alignment results.\n");
  }
  catch (resultsError &e) {
    fprintf(stderr, "Wrong alignment results.\n");
  }
}
