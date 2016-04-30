#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include "bwa/bwamem.h"
#include "bwa_wrapper.h"
#include "util.h"

template <typename T>
void assertEQ(int i, int j, const char* msg, T const &a, T const &b) {
  if (a != b) {
    std::stringstream ss;
    ss << "[" << i << "," << j << "] " << msg << ":"
       << a << "!=" << b << "\n";
    throw resultsError(ss.str().c_str());
  }
}

void regionsCompare(
    mem_alnreg_v *alnreg_base,
    mem_alnreg_v *alnreg_test,
    int num_seqs)
{
  try {
    for (int i = 0; i < num_seqs; i++) {
      // Check if the number of aligned regions match
      assertEQ(i, -1, "n", alnreg_base[i].n, alnreg_test[i].n);
      for (int j = 0; j < alnreg_base[i].n; j++) {
        mem_alnreg_t* region_base = &alnreg_base[i].a[j];
        mem_alnreg_t* region_test = &alnreg_test[i].a[j];

        assertEQ(i, j, "rb", region_base->rb, region_test->rb);
        assertEQ(i, j, "re", region_base->re, region_test->re);
        assertEQ(i, j, "qb", region_base->qb, region_test->qb);
        assertEQ(i, j, "qe", region_base->qe, region_test->qe);
        assertEQ(i, j, "rid", region_base->rid, region_test->rid);
        assertEQ(i, j, "score", region_base->score, region_test->score);
        assertEQ(i, j, "sub", region_base->sub, region_test->sub);
        assertEQ(i, j, "alt_sc", region_base->alt_sc, region_test->alt_sc);
        assertEQ(i, j, "csub", region_base->csub, region_test->csub);
        assertEQ(i, j, "sub_n", region_base->sub_n, region_test->sub_n);
        assertEQ(i, j, "w", region_base->w, region_test->w);
        assertEQ(i, j, "seedcov", region_base->seedcov, region_test->seedcov);
        assertEQ(i, j, "secondary", region_base->secondary, region_test->secondary);
        assertEQ(i, j, "secondary_all", region_base->secondary_all, region_test->secondary_all);
        assertEQ(i, j, "seedlen0", region_base->seedlen0, region_test->seedlen0);
        assertEQ(i, j, "frac_rep", region_base->frac_rep, region_test->frac_rep);
        assertEQ(i, j, "hash", region_base->hash, region_base->hash);
      }
    }
    fprintf(stderr, "Correct alignment results.\n");
  }
  catch (resultsError &e) {
    fprintf(stderr, "Wrong alignment results: %s\n", e.what());
  }
}

void smemCompare(
    smem_aux_t **smem_base,
    smem_aux_t **smem_test,
    int num 
)
{
  try {
     for (int i=0; i<num; i++){
     bwtintv_v *intv_base = &smem_base[i]->mem;
     bwtintv_v *intv_test = &smem_test[i]->mem;
     assertEQ (i,  -1,"n",intv_base->n,intv_test->n);
      for (int j =0;j < intv_base->n; j++){
        bwtintv_t intv_t_base = intv_base->a[j];
        bwtintv_t intv_t_test = intv_test->a[j];
        assertEQ(i,j, "x[0]",intv_t_base.x[0],intv_t_test.x[0]);
        assertEQ(i,j, "x[1]",intv_t_base.x[1],intv_t_test.x[1]);
        assertEQ(i,j, "x[2]",intv_t_base.x[2],intv_t_test.x[2]);
        assertEQ(i,j, "info",intv_t_base.info,intv_t_test.info);
      }
    } 
  fprintf(stderr,"Correct smem results\n");
 }
 catch (resultsError &e){
  fprintf(stderr,"wrong smem results\n");
 }
}

