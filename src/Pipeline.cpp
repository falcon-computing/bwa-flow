#include <algorithm>
#include <boost/asio.hpp>
#include <boost/function_types/result_type.hpp>
#include <boost/make_shared.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread.hpp>
#include <fstream>
#include <list>
#include <queue>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>

#include "bwa/utils.h"
#include "kflow/Queue.h"

#ifdef USE_HTSLIB
#include "htslib/ksort.h"
#endif

#include "bwa_wrapper.h"
#include "config.h"
#include "Pipeline.h"  
#include "util.h"
#include "allocation_wrapper.h"

// Comparator function for bam1_t records
#ifdef USE_HTSLIB
bool bam1_lt(const bam1_t* a, const bam1_t* b) {
  return ((uint64_t)a->core.tid<<32|(a->core.pos+1)<<1|bam_is_rev(a)) 
       < ((uint64_t)b->core.tid<<32|(b->core.pos+1)<<1|bam_is_rev(b));
}
#endif

static inline void trim_readno(kstring_t *s)
{
	if (s->l > 2 && s->s[s->l-2] == '/' && isdigit(s->s[s->l-1]))
		s->l -= 2, s->s[s->l] = 0;
}

void getKseqBatch(int chunk_size, 
    int *n_, 
    void *ks1_, void *ks2_,
    kseq_new_t *ks_buffer)
{
  kseq_t *ks = (kseq_t*)ks1_, *ks2 = (kseq_t*)ks2_;
  int size = 0, t = 0;

  while (true) {
    if (kseq_read_new(&ks_buffer[t], ks) <0) {
      break;
    }
    size += strlen(ks_buffer[t++].seq.s);
    if(ks2) {
      if(kseq_read_new(&ks_buffer[t], ks2) < 0) {
        fprintf(stderr,"[W::%s] the 2nd file has fewer sequences.\n", __func__);
        break;
      }
      size += strlen(ks_buffer[t++].seq.s);
    }
    if (size >= chunk_size && (t&1) ==0) break;
  }
  if (size == 0) {
    if (ks2 && kseq_read(ks2) >= 0)
      fprintf(stderr, "[W::%s] the 1st file has fewer sequences.\n", __func__);
  }
  *n_ = t;
}

const int kseq_buffer_size = 8;
kestrelFlow::Queue<kseq_new_t*, kseq_buffer_size+1> kseq_queue;

void KseqsRead::compute() {
  uint64_t num_seqs_produced = 0;
  // initialize kseq_queue, TODO calculate the size instead of the magic number
  for (int i = 0; i < kseq_buffer_size; i++) {
    kseq_new_t *ks_new = (kseq_new_t*)calloc(140000, sizeof(kseq_new_t));
    kseq_queue.push(ks_new);
  }
  while (true) {
    uint64_t start_ts = getUs();

    int batch_num = 0;
    kseq_new_t *ks_buffer ;
    kseq_queue.pop(ks_buffer);

    // Get the kseq batch
    getKseqBatch(10000000, &batch_num, aux->ks, aux->ks2, ks_buffer);
    if (batch_num == 0) break;

    DLOG_IF(INFO, VLOG_IS_ON(1)) << "Read " << batch_num << " seqs in "
            << getUs() - start_ts << " us";

    KseqsRecord record;
    record.ks_buffer = ks_buffer;
    record.batch_num = batch_num;
    record.start_idx = num_seqs_produced;
    pushOutput(record);
    num_seqs_produced += batch_num;
  }
}

SeqsRecord KseqsToBseqs::compute(KseqsRecord const & input) {

  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started KseqsToBseqs() for one input";
  uint64_t start_ts = getUs();
  
  int batch_num = input.batch_num;
  kseq_new_t *ks_buffer = input.ks_buffer;
  bseq1_t* seqs = (bseq1_t*)calloc(batch_num, sizeof(bseq1_t));

  for (int i=0; i<batch_num; i++) { 
    trim_readno(&ks_buffer[i].name);
    seqs[i].name = strdup(ks_buffer[i].name.s);
    seqs[i].comment = ks_buffer[i].comment.l ? strdup(ks_buffer[i].comment.s) : 0;
    seqs[i].seq = strdup(ks_buffer[i].seq.s);
    seqs[i].qual = ks_buffer[i].qual.l? strdup(ks_buffer[i].qual.s) : 0;
    seqs[i].l_seq = strlen(seqs[i].seq);
  }

  // return the ks_buffer to the queue
  kseq_queue.push(ks_buffer);

  // Finish the remain steps
	if (!aux->copy_comment) {
		for (int i = 0; i < batch_num; i++) {
			free(seqs[i].comment);
			seqs[i].comment = 0;
		}
    DLOG_IF(INFO, FLAGS_v >= 2 && input.start_idx == 0) << 
      "Do not append seq comment";
  }

  SeqsRecord record;
  record.start_idx = input.start_idx;
  record.batch_num = batch_num;
  record.seqs = seqs;
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished KseqsToBseqs() in " << 
                                  getUs() - start_ts << " us";
  return record;
}

void SeqsRead::compute() {

  int num_seqs_produced = 0;

  while (true) {
    uint64_t start_ts = getUs();

    // Read from file input, get mem_chains
    int batch_num = 0;
    bseq1_t *seqs = bseq_read(10000000, 
        &batch_num, aux->ks, aux->ks2);

    if (!seqs) break;

    if (!aux->copy_comment) {
	    for (int i = 0; i < batch_num; i++) {
		    free(seqs[i].comment);
		    seqs[i].comment = 0;
	    }
	    VLOG_IF(2, num_seqs_produced == 0) << "Do not append seq comment";
    }

    VLOG(1) << "Read " << batch_num << " seqs in "
	    << getUs() - start_ts << " us";

    // Construct output record
    SeqsRecord record;
    record.start_idx = num_seqs_produced;
    record.batch_num = batch_num;
    record.seqs = seqs;

    pushOutput(record);
    num_seqs_produced += batch_num;
  }
}

SeqsRecord SeqsToSams::compute(SeqsRecord const & input) {

  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started SeqsToSams() for one input";
  uint64_t start_ts = getUs();

  bseq1_t* seqs = input.seqs;
  uint64_t start_idx = input.start_idx;
  int batch_num = input.batch_num;

  mem_alnreg_v* alnreg = new mem_alnreg_v[batch_num];
  if (NULL == alnreg) {
    LOG(ERROR) << strerror(errno) << " due to "
               << ((errno==12) ? "out-of-memory" : "internal failure") ;
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < batch_num; i++) {
    mem_chain_v chains = seq2chain(aux, &seqs[i]);
    kv_init(alnreg[i]);
    for (int j = 0; j < chains.n; j++) {
      mem_chain2aln(
          aux->opt, 
          aux->idx->bns, 
          aux->idx->pac,
          seqs[i].l_seq,
          (uint8_t*)seqs[i].seq,
          &chains.a[j],
          &alnreg[i]);

      free(chains.a[j].seeds);
    }
    free(chains.a);

    // Post-process each chain before output
    alnreg[i].n = mem_sort_dedup_patch(
        aux->opt, 
        aux->idx->bns, 
        aux->idx->pac, 
        (uint8_t*)seqs[i].seq, 
        alnreg[i].n, 
        alnreg[i].a);

    for (int j = 0; j < alnreg[i].n; j++) {
      mem_alnreg_t *p = &alnreg[i].a[j];
      if (p->rid >= 0 && aux->idx->bns->anns[p->rid].is_alt)
        p->is_alt = 1;
    }
  }
  mem_pestat_t pes[4];
  mem_pestat(aux->opt, aux->idx->bns->l_pac, batch_num, alnreg, pes);

#ifdef USE_HTSLIB
  for (int i =0; i< batch_num/2; i++) {
    seqs[i<<1].bams = bams_init();
    seqs[1+(i<<1)].bams = bams_init();
    mem_sam_pe(
        aux->opt,
        aux->idx->bns,
        aux->idx->pac,
        pes,
        (start_idx>>1)+i,
        &seqs[i<<1],
        &alnreg[i<<1],
        aux->h);
     }
#else
  for (int i = 0; i < batch_num/2; i++) {
    mem_sam_pe(
        aux->opt,
        aux->idx->bns,
        aux->idx->pac,
        pes,
        (start_idx>>1)+i,
        &seqs[i<<1],
        &alnreg[i<<1]);
  }
#endif
  freeAligns(alnreg, batch_num);

  // Free fields in seq except sam
  for (int i = 0; i < batch_num; i++) {
    free(seqs[i].name);
    free(seqs[i].comment);
    free(seqs[i].seq);
    free(seqs[i].qual);
  }

  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished SeqsToSams() for one batch in "
    << getUs() - start_ts << " us";

  return input;
}

ChainsRecord SeqsToChains::compute(SeqsRecord const & seqs_record) {

  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started SeqsToChains() for one input";

  uint64_t start_ts = getUs();
  uint64_t ref_time = 0;

  bseq1_t* seqs = seqs_record.seqs;
  uint64_t start_idx = seqs_record.start_idx;
  int batch_num = seqs_record.batch_num;

  mem_chain_v* chains = (mem_chain_v*)malloc(batch_num*sizeof(mem_chain_v));

  for (int i = 0; i < batch_num; i++) {
    chains[i] = seq2chain(aux, &seqs[i]);
  }

  ChainsRecord ret;
  ret.start_idx    = seqs_record.start_idx;
  ret.batch_num    = batch_num;
  ret.seqs         = seqs;
  ret.chains       = chains;
 
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished SeqsToChains() in " 
    << getUs() - start_ts << " us";

  return ret;
}

RegionsRecord ChainsToRegions::compute(ChainsRecord const & record) {

  uint64_t start_ts = getUs();
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started ChainsToRegions() on CPU";

  bseq1_t* seqs       = record.seqs;
  mem_chain_v* chains = record.chains;
  int batch_num       = record.batch_num;

  mem_alnreg_v* alnreg = (mem_alnreg_v*)malloc(batch_num*sizeof(mem_alnreg_v));

  for (int i = 0; i < batch_num; i++) {
    kv_init(alnreg[i]);
    for (int j = 0; j < chains[i].n; j++) {
      mem_chain2aln(
          aux->opt, 
          aux->idx->bns, 
          aux->idx->pac,
          seqs[i].l_seq,
          (uint8_t*)seqs[i].seq,
          &chains[i].a[j],
          alnreg+i);
    }
  }

  RegionsRecord output;
  output.start_idx = record.start_idx;
  output.batch_num = batch_num;
  output.seqs = seqs;
  output.alnreg = alnreg;

  freeChains(chains, batch_num);

  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished ChainsToRegions() on CPU for "
    << getUs() - start_ts << " us";

  return output;
}

SeqsRecord RegionsToSam::compute(RegionsRecord const & record) {

  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started RegionsToSam() for one input";

  uint64_t start_ts = getUs();
  uint64_t seedcov_time = 0;

  uint64_t start_idx   = record.start_idx;
  int batch_num        = record.batch_num;
  mem_alnreg_v* alnreg = record.alnreg;
  bseq1_t* seqs        = record.seqs;

  for (int i = 0; i < batch_num; i++) {
    // Post-process each chain before output
    alnreg[i].n = mem_sort_dedup_patch(
        aux->opt, 
        aux->idx->bns, 
        aux->idx->pac, 
        (uint8_t*)seqs[i].seq, 
        alnreg[i].n, 
        alnreg[i].a);

    for (int j = 0; j < alnreg[i].n; j++) {
      mem_alnreg_t *p = &alnreg[i].a[j];
      if (p->rid >= 0 && aux->idx->bns->anns[p->rid].is_alt)
        p->is_alt = 1;
    }
  }

  if(aux->opt->flag&MEM_F_PE) {
    mem_pestat_t pes[4];
    mem_pestat(aux->opt, aux->idx->bns->l_pac, batch_num, alnreg, pes);
#ifdef USE_HTSLIB
    for (int i =0; i< batch_num/2; i++) {
      seqs[i<<1].bams = bams_init();
      seqs[1+(i<<1)].bams = bams_init();
      mem_sam_pe(
          aux->opt,
          aux->idx->bns,
          aux->idx->pac,
          pes,
          (start_idx>>1)+i,
          &seqs[i<<1],
          &alnreg[i<<1],
          aux->h);
       }
#else
    for (int i = 0; i < batch_num/2; i++) {
      mem_sam_pe(
          aux->opt,
          aux->idx->bns,
          aux->idx->pac,
          pes,
          (start_idx>>1)+i,
          &seqs[i<<1],
          &alnreg[i<<1]);
    }
#endif
  }
  else {
    for (int i=0; i<batch_num; i++) {
#ifdef USE_HTSLIB
      seqs[i].bams = bams_init();
#endif
      mem_mark_primary_se(
          aux->opt,
          alnreg[i].n,
          alnreg[i].a,
          start_idx+i
          );
      mem_reg2sam(
          aux->opt,
          aux->idx->bns,
          aux->idx->pac,
          &seqs[i],
          &alnreg[i],
          0,
          0,
          aux->h
          );
    }
  }

  freeAligns(alnreg, batch_num);

  // Free fields in seq except sam
  for (int i = 0; i < batch_num; i++) {
    free(seqs[i].name);
    free(seqs[i].comment);
    free(seqs[i].seq);
    free(seqs[i].qual);
  }

  SeqsRecord output;
  output.start_idx = start_idx;
  output.batch_num = batch_num;
  output.seqs = seqs;

  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished RegionsToSam() in " << getUs() - start_ts << " us";
  return output;
}

#ifdef USE_HTSLIB
typedef bam1_t *bam1_p;
KSORT_INIT(sort, bam1_p, bam1_lt)
#endif

void SamsReorder::compute(int wid) {

  uint64_t n_processed = 0;
  std::unordered_map<uint64_t, SeqsRecord> record_buf;
#ifdef USE_HTSLIB
  int max_batch_records = FLAGS_max_batch_records;
  bam1_t** bam_buffer = NULL;
  int bam_buffer_idx = 0;
  int batch_records = 0;
  int bam_buffer_order = 0;

  // bams->l may be larger than 1
  // the bam_buffer needs to be realloc
  int max_buffer_records = 0;

  BamsRecord output;
#else
  SeqsRecord output;
#endif
  while (true) {
    SeqsRecord input;
    SeqsRecord record;
    bool ready = this->getInput(input);
    while (!this->isFinal() && !ready) {
      boost::this_thread::sleep_for(boost::chrono::microseconds(100));
      ready = this->getInput(input);
    }
    if (!ready) { 
      // this means isFinal() is true and input queue is empty
      break; 
    }
    if (FLAGS_inorder_output) {
      record_buf[input.start_idx] = input;

      while (record_buf.count(n_processed)) {
        record = record_buf[n_processed];
        record_buf.erase(n_processed);
        n_processed += record.batch_num;
#ifdef USE_HTSLIB
        int batch_num = record.batch_num;
        bseq1_t* seqs = record.seqs;
    
        if (!bam_buffer) {
          max_buffer_records = max_batch_records*batch_num*2;
          bam_buffer = (bam1_t**)malloc(max_buffer_records*sizeof(bam1_t*));
        }
        for (int i = 0; i < batch_num; i++) {
          if (seqs[i].bams) {
            for (int j =0; j < seqs[i].bams->l; j++) {
              bam1_t* bam_record = seqs[i].bams->bams[j];
              if (FLAGS_filter == 0 || 
                  (bam_record->core.flag & FLAGS_filter) == 0) 
              {
                if (bam_buffer_idx >= max_buffer_records) {
                  max_buffer_records *= 2;
                  bam_buffer = (bam1_t**)realloc(bam_buffer, max_buffer_records*sizeof(bam1_t*)); 
                }
                bam_buffer[bam_buffer_idx++] = bam_record;
              }
            }
          }
          free(seqs[i].bams->bams);
          free(seqs[i].bams); seqs[i].bams = NULL;    
        }
        free(seqs);
 
        if (++batch_records >= max_batch_records) {
          if(FLAGS_sort) {
            uint64_t start_ts_st = getUs();
            //std::sort(bam_buffer, bam_buffer+bam_buffer_idx, bam1_lt);
            ks_mergesort(sort, bam_buffer_idx, (bam1_p *)bam_buffer, 0);
            DLOG_IF(INFO, VLOG_IS_ON(3)) << "Sorted " << bam_buffer_idx << " records in " <<
              getUs() - start_ts_st << " us";
          }
          output.bam_buffer = bam_buffer;
          output.bam_buffer_idx = bam_buffer_idx;
          output.bam_buffer_order = bam_buffer_order;
          bam_buffer_idx = 0;
          batch_records = 0;
          bam_buffer_order = bam_buffer_order + 1;
          pushOutput(output);
          bam_buffer = NULL;
        }
#else
        output = record;
        pushOutput(output);
#endif
      }
    }
    else {
      record = input;
#ifdef USE_HTSLIB
      int batch_num = record.batch_num;
      bseq1_t* seqs = record.seqs;
  
      if (!bam_buffer) {
        // in case bams->l is not 1 
        max_buffer_records = max_batch_records*batch_num*2;
        bam_buffer = (bam1_t**)malloc(max_buffer_records*sizeof(bam1_t*));
      }
      for(int i = 0; i < batch_num; i++) {
        if (seqs[i].bams) {
          for (int j =0; j < seqs[i].bams->l; j++) {
            bam1_t* bam_record = seqs[i].bams->bams[j];
            if (FLAGS_filter == 0 || 
                (bam_record->core.flag & FLAGS_filter) == 0) 
            {
              if (bam_buffer_idx >= max_buffer_records) {
                max_buffer_records *= 2;
                bam_buffer = (bam1_t**)realloc(bam_buffer, max_buffer_records*sizeof(bam1_t*));
                DLOG(INFO) << "realloc() bam_buffer to " << max_buffer_records;
              }
              bam_buffer[bam_buffer_idx++] = bam_record;
            }
          }
        }
        free(seqs[i].bams->bams);
        free(seqs[i].bams); seqs[i].bams = NULL;    
      }
      free(seqs);
      if (++ batch_records >=  max_batch_records) {
        if(FLAGS_sort) {
          uint64_t start_ts_st = getUs();
          //std::sort(bam_buffer, bam_buffer+bam_buffer_idx, bam1_lt);
          ks_mergesort(sort, bam_buffer_idx, (bam1_p *)bam_buffer, 0);
          DLOG_IF(INFO, VLOG_IS_ON(3)) << "Sorted " << bam_buffer_idx << " records in " <<
            getUs() - start_ts_st << " us";
        }
        output.bam_buffer = bam_buffer;
        output.bam_buffer_idx = bam_buffer_idx;
        output.bam_buffer_order = bam_buffer_order;
        bam_buffer_idx = 0;
        batch_records = 0;
        bam_buffer_order = bam_buffer_order + 1;
        pushOutput(output);
        bam_buffer = NULL;
      }
#else
      output = record;
      pushOutput(output);
#endif
    }
  }
#ifdef USE_HTSLIB
  //finish the remaining sam
  if(batch_records > 0){
    if(FLAGS_sort) {
      uint64_t start_ts_st = getUs();
      //std::sort(bam_buffer, bam_buffer+bam_buffer_idx, bam1_lt);
      ks_mergesort(sort, bam_buffer_idx, (bam1_p *)bam_buffer, 0);
      DLOG_IF(INFO, VLOG_IS_ON(3)) << "Sorted " << bam_buffer_idx << " records in " <<
        getUs() - start_ts_st << " us";
    }
    output.bam_buffer = bam_buffer;
    output.bam_buffer_idx = bam_buffer_idx;
    output.bam_buffer_order = bam_buffer_order;
    bam_buffer_idx = 0;
    batch_records = 0;
    bam_buffer_order = bam_buffer_order + 1;
    pushOutput(output);
    bam_buffer = NULL; 
  }
#endif
}

#ifdef USE_HTSLIB
int WriteOutput::compute(BamsRecord const &input)
#else
int WriteOutput::compute(SeqsRecord const &input)
#endif
{
  boost::any var = this->getConst("sam_dir");
  std::string out_dir = boost::any_cast<std::string>(var);
  bool use_file = !out_dir.empty();
#ifdef USE_HTSLIB
  // write bam output
  if (!use_file) {
    LOG(ERROR) << "Bams output only works with file output,"
      << " please specify --output_dir";
    exit(1);
  }
  const char *modes[] = {"wb", "wb0", "w"};
  samFile *fout = NULL;

  uint64_t start_ts = getUs();
  int file_id = input.bam_buffer_order;
  int bam_buffer_idx = input.bam_buffer_idx;
  bam1_t** bam_buffer = input.bam_buffer;

  std::stringstream ss;
  ss << out_dir << "/part-" << std::setw(6) << std::setfill('0') << file_id;
  DLOG_IF(INFO, VLOG_IS_ON(2)) << "Writting to " << ss.str();
  fout = sam_open(ss.str().c_str(), modes[FLAGS_output_flag]); 
  if (!fout) {
    throw std::runtime_error("Cannot open bam output file");
  }
  int status = sam_hdr_write(fout, aux->h);
  if (status) {
    LOG(ERROR) << "sam_hdr_write error: " << status;
  }
  // start writing to file
  for (int i = 0; i < bam_buffer_idx; ++i){
    size_t n = sam_write1(fout, aux->h, bam_buffer[i]); 
    bam_destroy1(bam_buffer[i]);
  }
  sam_close(fout);
  free(bam_buffer);
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Written " << bam_buffer_idx
          << " records in "
          << getUs() - start_ts << " us";
  return 0;
#else
  // write sam output
  FILE* fout;
  fout = stdout;
  uint64_t start_ts = getUs();
  int batch_num = input.batch_num;
  bseq1_t* seqs = input.seqs;
  for (int i = 0; i < batch_num; ++i) {
    if (seqs[i].sam) fputs(seqs[i].sam, fout);
    free(seqs[i].sam);
  }
  free(seqs);
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Written batch " << input.start_idx 
    << " to file in " << getUs() - start_ts << " us";
#endif
}
