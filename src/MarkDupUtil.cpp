#include <cstring>

#include "util.h"

#define DUP_FLAG 1024
#define NP_FLAG 256

splitLine_t * bamToSplitLine(ktp_aux_t* aux, bam1_t* bam_record) {
  splitLine* sline = getSplitLine();
  kstring_t ks = { 0, 0, NULL };
  sam_format1(aux->h, bam_record, &ks);
  sline->bufLen = ks.l;
  strcpy(sline->buffer, ks.s);
  free(ks.s);
  splitSplitLine(sline, 12);
  return sline;
}

splitLine_t * readSeq(ktp_aux_t* aux, bseq1_t seq) {
  if (seq.bams->bams[0] != NULL) {
    if (seq.bams->bams[0]->core.flag & NP_FLAG) {
      LOG(INFO)<<"alignment not primary.";
    }
    splitLine* line = bamToSplitLine(aux, seq.bams->bams[0]);
    if (line->bufLen < 1) {
      return NULL;
    }
    return line;
  }
  return NULL;
}

bool checkSplitLineDup(splitLine_t * sline) {
  return (bool)(sline->flag & DUP_FLAG);
}

void markDupSeq(bseq1_t* seq){
  for(int i = 0; i < seq->bams->l; i++) {
    if (i >= 2) {
      DLOG(INFO)<<"more than 1 alignments marked: "<<i;
    }
    if(seq->bams->bams[i] != NULL){
      seq->bams->bams[i]->core.flag = (seq->bams->bams[i]->core.flag | DUP_FLAG);
    }
  }
  return;
}
