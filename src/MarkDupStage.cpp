#include "MarkDupStage.h"
#include <cstring>

#include "util.h"

#define DUP_FLAG 1024
#define NP_FLAG 256

static splitLine_t * bamToSplitLine(ktp_aux_t* aux, bam1_t* bam_record) {
  splitLine* sline = getSplitLine();
  kstring_t ks = { 0, 0, NULL };
  sam_format1(aux->h, bam_record, &ks);
  sline->bufLen = ks.l;
  strcpy(sline->buffer, ks.s);
  free(ks.s);
  splitSplitLine(sline, 12);
  return sline;
}

static splitLine_t * readSeq(ktp_aux_t* aux, bseq1_t seq) {
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

static bool checkSplitLineDup(splitLine_t * sline) {
  return (bool)(sline->flag & DUP_FLAG);
}

static void markDupSeq(bseq1_t* seq) {
  for(int i = 0; i < seq->bams->l; i++) {
    if (i >= 2) {
      DLOG(INFO)<<"more than 1 alignments marked: "<<i;
    }
    if(seq->bams->bams[i] != NULL) {
      seq->bams->bams[i]->core.flag = (seq->bams->bams[i]->core.flag | DUP_FLAG);
    }
  }
  return;
}

void MarkDupStage::InitializeState(ktp_aux_t* aux) {
  state = makeState();
  state->seqLens = (UINT32*)calloc(1, sizeof(UINT32));
  state->seqOffs = (UINT64*)calloc(1, sizeof(UINT64));
  state->seqs[strdup("*")] = 0;
  state->seqLens[0] = padLength(0);
  state->seqOffs[0] = 0;
#ifdef USE_HTSLIB
  UINT64 totalLen = 0;
  for(int i = 0; i < aux->h->n_targets; i++) {
    char * seqID = aux->h->target_name[i];;
    UINT32 seqLen = (UINT32)aux->h->target_len[i];
    UINT64 seqOff = totalLen;
    totalLen += (UINT64)(seqLen + 1);
    if(i % 32768 == 1) {
      state->seqLens = (UINT32*)realloc(state->seqLens, (i + 32768)*sizeof(UINT32));
      state->seqOffs = (UINT64*)realloc(state->seqOffs, (i + 32768)*sizeof(UINT64));
    }
    state->seqs[strdup(seqID)] = i;
    state->seqLens[i] = seqLen;
    state->seqOffs[i] = seqOff;
  }
  int binCount = (totalLen >> BIN_SHIFT);
  if (binCount >= (1 << 15)) {
    //Error Too many sequences in header of input sam file.
  }
  state->binCount = binCount;
  state->sigArraySize = (binCount * 2 + 1) * (binCount * 2 + 1) + 1;
  //state->sigs = (sigSet_t *) malloc(state->sigArraySize * sizeof(sigSet_t));
  state->sigs = (sigSet_t *) new sigSet_t[state->sigArraySize];
  if (state->sigs == NULL) fatalError("samblaster: Unable to allocate signature set array.");
  //for (UINT32 i=0; i<state->sigArraySize; i++) state->sigs[i].hashTableInit(); 
#else
#endif
}

SeqsRecord MarkDupStage::compute(SeqsRecord const & input) {
//uint64_t read_seq_time = 0;
//uint64_t mark_dup_time = 0;
#if 0
  while (true) {
    SeqsRecord input;
    bool ready = this->getInput(input);
    while (!this->isFinal() && !ready) {
      boost::this_thread::sleep_for(boost::chrono::microseconds(10));
      ready = this->getInput(input);
    }
    if (!ready) {
      break;
    }
#endif
    //uint64_t all_start = getUs();
    DLOG(INFO) << "Started MarkDup()";
    //uint64_t read_seq_s = getUs();
#if 0 //for bamsRecord input
    std::vector<SeqsRecord> * seqsRecord = input.records_list;
#endif
    //for SeqsReord input
    std::vector<SeqsRecord> tmp_vec;
    tmp_vec.push_back(input);
    std::vector<SeqsRecord> * seqsRecord = &tmp_vec;
    //uint64_t read_seq_e = getUs();
    //read_seq_time += (read_seq_e - read_seq_s);
    for (int k =0; k < (*seqsRecord).size(); k++) {
      int batch_num = (*seqsRecord)[k].batch_num;
      splitLine_t* line = readSeq(aux, (*seqsRecord)[k].seqs[0]);
      if (line != NULL) {
        //DLOG(INFO) << "fisrt read is not NULL.";
      }
      splitLine_t* head = line;
      int count = 1;
      splitLine_t* last = line;
      for (int i = 1; i < batch_num; i++) {
        //read_seq_s = getUs();
        splitLine_t* nextLine =  readSeq(aux, (*seqsRecord)[k].seqs[i]);
        //read_seq_e = getUs();
        //read_seq_time += (read_seq_e - read_seq_s);
        if (nextLine == NULL) {
          break;
        }
        if (strcmp(line->fields[QNAME], nextLine->fields[QNAME]) != 0) {
          //DLOG(INFO) << "before md " << line->fields[QNAME] << " " << nextLine->fields[QNAME];
          mtx_.lock();
          markDupsDiscordants(line, state);
          mtx_.unlock();
          splitLine_t* tmp = line;
          for (int j = 0; j < count; j++) {
            if (checkSplitLineDup(tmp)) {
              //uint64_t mark_dup_s = getUs();
              markDupSeq(&((*seqsRecord)[k].seqs[i- count +j]));
              //uint64_t mark_dup_e = getUs();
              //mark_dup_time += (mark_dup_e - mark_dup_s);
            }
            splitLine* release = tmp;
            tmp = tmp->next;
            deleteSplitLine(release);
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

      mtx_.lock();
      markDupsDiscordants(line, state);
      mtx_.unlock();
      splitLine_t* tmp = line;
      for(int j = 0; j < count; j++) {
        if (checkSplitLineDup(tmp)) {
          //uint64_t mark_dup_s = getUs();
          markDupSeq(&((*seqsRecord)[k].seqs[batch_num - count + j]));
          //uint64_t mark_dup_e = getUs();
          //mark_dup_time += (mark_dup_e - mark_dup_s);
        }
        splitLine* release = tmp;
        tmp = tmp->next;
        deleteSplitLine(release);
      }
    }
    //pushOutput(input); // for mapPartitionStage
    DLOG(INFO) << "Finished MarkDup()";
    //uint64_t all_end = getUs();
    //uint64_t all_diff = all_end - all_start;
    //DLOG(INFO) << "MdStage AllTime = " << all_diff;
    //DLOG(INFO) << "ReadSeq Time = " << read_seq_time;
    //DLOG(INFO) << "MarkDup Time = " << mark_dup_time;
    return input;
  // } // for MapPartitionStage
#if 0
  tmp = head;
  while(tmp->next != NULL) {
    splitLine_t* release = tmp;
    tmp = tmp->next;
    deleteSplitLine(release);
  }
  deleteSplitLine(tmp);
#endif
  //mtx_.unlock();
}
