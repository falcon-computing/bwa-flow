#include <stdio.h>
#include <stdlib.h>
#include "bwa/bwamem.h"

#include "util.h"

template <typename T>
void assertEQ(T const &a, T const &b) {
  if (a != b) {
    throw resultsError("");
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
