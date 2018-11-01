#include "mdPipeline.h"
#include <cstring>

#define DUP_FLAG 1024
#define NP_FLAG 256

splitLine_t * bamToSplitLine(ktp_aux_t* aux, bam1_t* bam_record) {
  splitLine* sline = getSplitLine();
  kstring_t ks = { 0, 0, NULL };
  sam_format1(aux->h, bam_record, &ks);
  sline->bufLen = ks.l;
  sline->buffer = ks.s;
  splitSplitLine(sline, 12);
  return sline;
}

splitLine_t * readSeq(ktp_aux_t* aux, bseq1_t seq) {
  //LOG(INFO)<<"warning: enter readSeq function.";
  if (seq.bams->bams[0] != NULL) {
    //LOG(INFO)<<"start reading alignemnt.";
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
  if(seq->bams->bams[0] != NULL){
    //DLOG(INFO)<<"mark duplicate alignment.";
    seq->bams->bams[0]->core.flag = (seq->bams->bams[0]->core.flag | DUP_FLAG);
    return;
  }
  return;
}

void Markdup::InitializeState(ktp_aux_t* aux) {
  state = makeState();
  state->seqLens = (UINT32*)calloc(1, sizeof(UINT32));
  state->seqOffs = (UINT64*)calloc(1, sizeof(UINT64));
  state->seqs[strdup("*")] = 0;
  state->seqLens[0] = padLength(0);
  state->seqOffs[0] = 0;
#ifdef USE_HTSLIB
  UINT64 totalLen = 0;
  for(int i = 0; i < aux->h->n_targets; i++){
    char * seqID = aux->h->target_name[i];;
    UINT32 seqLen = (UINT32)aux->h->target_len[i];
    UINT64 seqOff = totalLen;
    totalLen += (UINT64)(seqLen + 1);
    if(i % 32768 == 1){
      state->seqLens = (UINT32*)realloc(state->seqLens, (i + 32768)*sizeof(UINT32));
      state->seqOffs = (UINT64*)realloc(state->seqOffs, (i + 32768)*sizeof(UINT64));
    }
    state->seqs[strdup(seqID)] = i;
    state->seqLens[i] = seqLen;
    state->seqOffs[i] = seqOff;
  }
  int binCount = (totalLen >> BIN_SHIFT);
  if (binCount >= (1 << 15)){
    //Error Too many sequences in header of input sam file.
  }
  state->binCount = binCount;
  state->sigArraySize = (binCount * 2 + 1) * (binCount * 2 + 1) + 1;
  state->sigs = (sigSet_t *) malloc(state->sigArraySize * sizeof(sigSet_t));
  if (state->sigs == NULL) fatalError("samblaster: Unable to allocate signature set array.");
  for (UINT32 i=0; i<state->sigArraySize; i++) hashTableInit(&(state->sigs[i])); 
#else
#endif
}

SeqsRecord Markdup::compute(SeqsRecord const & input) {
  //mtx_.lock();
  DLOG(INFO) << "Started MarkDup()";
  //DLOG(INFO) << "start comptuting md.";
  int batch_num = input.batch_num;
  bseq1_t* seqs = input.seqs;
  splitLine_t* line = readSeq(aux, seqs[0]);
  //DLOG(INFO) << "read first line.";
  if (line != NULL) {
    //DLOG(INFO) << "fisrt read is not NULL.";
  }
  splitLine_t* head = line;
  int count = 1;
  splitLine_t* last = line;
  for (int i = 1; i < batch_num; i++) {
    splitLine_t* nextLine =  readSeq(aux, seqs[i]);
    if (nextLine == NULL) {
      break;
    }
    if (strcmp(line->fields[QNAME], nextLine->fields[QNAME]) != 0) {
      //DLOG(INFO) << "before md " << line->fields[QNAME] << " " << nextLine->fields[QNAME];
      markDupsDiscordants(line, state);
      splitLine_t* tmp = line;
      for (int j = 0; j < count; j++) {
        if (checkSplitLineDup(tmp)) {
          markDupSeq(seqs + i - count + j);
        }
        tmp = tmp->next;
      }
      last = line = nextLine;
      count = 1;
    }
    else{
      last->next = nextLine;
      last = nextLine;
      count += 1;
    }
  } 

  markDupsDiscordants(line, state);
  splitLine_t* tmp = line;
  for(int j = 0; j < count; j++){
    if (checkSplitLineDup(tmp)) {
      markDupSeq(seqs + batch_num - count + j);
    }
    tmp = tmp->next;
  }

  tmp = head;
  while(tmp->next != NULL){
    splitLine_t* release = tmp;
    tmp = tmp->next;
    deleteSplitLine(release);
  }
  deleteSplitLine(tmp);
  DLOG(INFO) << "Finished MarkDup()";
  //mtx_.unlock();
  return input;
}
