#include "config.h"
#include "SWRead.h"
#include "util.h"


SWRead::SWRead(int start_idx, int idx,
    ktp_aux_t* aux,
    const bseq1_t* seq, 
    const mem_chain_v* chains,
    mem_alnreg_v* alnregs,
    mem_chainref_t* ref,
    std::vector<int>* chain_idxes):
  is_pend_(false), 
  start_idx_(start_idx),
  read_idx_(idx), chain_idx_(0),
  aux_(aux), seq_(seq), chains_(chains),
  alnregs_(alnregs), ref_(ref),
  chain_idxes_(chain_idxes)
{
  if (chains && chains->n > 0) {
    seed_idx_ = chains->a[0].n - 1;
    is_finished_ = false;
  }
  else {
    seed_idx_ = -1;
    is_finished_ = true;
  }
}

SWRead::~SWRead() {
  for (int j = 0; j < chains_->n; j++) {
      free(ref_[j].srt);
      free(ref_[j].rseq);
  }
  free(ref_);
}

void SWRead::finish() {
  is_pend_ = false;
}

inline void SWRead::getChainRef(
      const ktp_aux_t* aux,
      const bseq1_t* seq,
      const mem_chain_v* chain,
      mem_chainref_t* &ref
) {
  int64_t l_pac = aux->idx->bns->l_pac;

  // input for SmithWaterman for each chain
  ref = (mem_chainref_t*)malloc(chain->n*sizeof(mem_chainref_t));

  for (int j = 0; j < chain->n; j++) {
    // Prepare the maxspan and rseq for each seed
    mem_chain_t *c = &chain->a[j]; 

    int64_t rmax[2], tmp, max = 0;
    uint8_t *rseq = 0;
    uint64_t *srt;

    if (c->n == 0) {
      continue;
    }
    // get the max possible span
    rmax[0] = l_pac<<1; rmax[1] = 0;
    for (int k = 0; k < c->n; ++k) {
      int64_t b, e;
      const mem_seed_t *t = &c->seeds[k];
      b = t->rbeg - (t->qbeg + cal_max_gap(aux->opt, t->qbeg));
      e = t->rbeg + t->len + ((seq->l_seq - t->qbeg - t->len)
          + cal_max_gap(aux->opt, seq->l_seq - t->qbeg - t->len));
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
    int rid;
    rseq = bns_fetch_seq(aux->idx->bns, aux->idx->pac, &rmax[0], c->seeds[0].rbeg, &rmax[1], &rid);
    assert(c->rid == rid);

    srt = (uint64_t *)malloc(c->n * 8);
    for (int l = 0; l < c->n; ++l)
      srt[l] = (uint64_t)c->seeds[l].score<<32 | l;
    ks_introsort_64(c->n, srt);

    ref[j].rmax[0] = rmax[0];
    ref[j].rmax[1] = rmax[1];
    ref[j].rseq = rseq;
    ref[j].srt = srt;
  }
}


//void SWRead::packRead(char* buffer, int &buffer_idx, int &task_num, int &read_num, mem_alnreg_t** &region_batch) {
void SWRead::packRead(char* buffer, int &buffer_idx, int &task_num, int &read_num,
    mem_alnreg_t** &region_batch, mem_chainref_t* &chainref_p) {

  getChainRef(aux, seq_, chains_, ref_);
  chainref_p = ref_;
  if (is_finished_) {
    return ;
  }
  //TODO put this part in the later loop, should be verified  
  for ( ; chain_idx_ < chains_->n; 
          seed_idx_ > 0 ? seed_idx_-- : seed_idx_ = chains_->a[++chain_idx_].n-1)
  {
      uint32_t sorted_idx = (uint32_t)(ref_[chain_idx_].srt[seed_idx_]);

    // get next available seed in the current read
      mem_seed_t* seed_array = &chains_->a[chain_idx_].seeds[sorted_idx];
      // initialize the newreg
      mem_alnreg_t* newreg = kv_pushp(mem_alnreg_t, *alnregs_);
      memset(newreg, 0, sizeof(mem_alnreg_t));

      // Push chain_idx to added newreg
      chain_idxes_->push_back(chain_idx_);

      newreg->score  = seed_array->len * aux_->opt->a;
      newreg->truesc = seed_array->len * aux_->opt->a;
      newreg->qb = seed_array->qbeg;
      newreg->rb = seed_array->rbeg;
      newreg->qe = seed_array->qbeg + seed_array->len;
      newreg->re = seed_array->rbeg + seed_array->len; 
      newreg->rid = chains_->a[chain_idx_].rid;
      newreg->seedlen0 = seed_array->len; 
      newreg->frac_rep = chains_->a[chain_idx_].frac_rep;
      newreg->w = aux_->opt->w;
  }
  chain_idx_ = 0;
  seed_idx_ = 0;
  int counter8 = 0;
  int tmp_int = 0;
  int chain_num = 0;
  int seed_num = 0;
  read_num += 1;
  int chain_num_addr = 0;
  int seed_num_addr = 0;
  int idx_end_addr = buffer_idx;
  buffer_idx += 4;
  // pack the read sequence
  *((int*)(&buffer[buffer_idx])) = seq_->l_seq;
  buffer_idx += 4;
  for ( int i = 0; i < seq_->l_seq; i++ ){
    counter8 += 1;
    tmp_int = tmp_int << 4 | seq_->seq[i];
    if ( counter8 % 8 ==0 ) {
      *(uint32_t*) (&buffer[buffer_idx]) = tmp_int ;
      buffer_idx += 4 ;
    }
  }
  if ( counter8 % 8 !=0 ) {
    tmp_int = tmp_int << (4*(8 - counter8 % 8));
    *(uint32_t*) (&buffer[buffer_idx]) = tmp_int ;
    buffer_idx += 4 ;
  }
  counter8 = 0;
  tmp_int = 0;
  chain_num_addr = buffer_idx;
  buffer_idx += 4;
  int region_num = 0; // used to keep track of the current region
  for ( chain_idx_ = 0; chain_idx_ < chains_->n; chain_idx_++) {
    // Pack the maxspan and rseq
    *((int64_t*)(&buffer[buffer_idx]))= ref_[chain_idx_].rmax[0];
    buffer_idx += 8;
    *((int64_t*)(&buffer[buffer_idx]))= ref_[chain_idx_].rmax[1];
    buffer_idx += 8;
    for ( int i = 0; i < ref_[chain_idx_].rmax[1] - ref_[chain_idx_].rmax[0]; i++) {
      counter8 = counter8 + 1;
      tmp_int = tmp_int << 4 | ref_[chain_idx_].rseq[i];
      if ( counter8 % 8 ==0 ) {
        *(uint32_t*) (&buffer[buffer_idx]) = tmp_int ;
        buffer_idx += 4 ;
      }
    } 
    if ( counter8 % 8 !=0 ) {
      tmp_int = tmp_int << (4*(8 - counter8 % 8));
      *(uint32_t*) (&buffer[buffer_idx]) = tmp_int ;
      buffer_idx += 4 ;
    }
    counter8 = 0;
    // record the address of seed number
    seed_num_addr = buffer_idx ;
    seed_num = 0;
    buffer_idx += 4 ; 
    // pack the seed information
    for ( seed_idx_ = chains_->a[chain_idx_].n -1 ; seed_idx_ >= 0 ; seed_idx_--) {
      uint32_t sorted_idx = (uint32_t)(ref_[chain_idx_].srt[seed_idx_]);
      // get next available seed in the current read
      mem_seed_t* seed_array = &chains_->a[chain_idx_].seeds[sorted_idx];
     
      if (seed_array->qbeg > 0 ||
          seed_array->qbeg + seed_array->len != seq_->l_seq)
      {
        region_batch[task_num] = &(alnregs_->a[region_num]) ;
        *((int*)(&buffer[buffer_idx])) = task_num ; 
        buffer_idx += 4;
        task_num += 1;
        seed_num += 1;
        *((int64_t*)(&buffer[buffer_idx])) = seed_array->rbeg ;
        buffer_idx += 8;
        *((int32_t*)(&buffer[buffer_idx])) = seed_array->qbeg ;
        buffer_idx += 4;
        *((int32_t*)(&buffer[buffer_idx])) = seed_array->len ;
        buffer_idx += 4;
      }
      region_num++;
    }
     *((int*)(&buffer[seed_num_addr])) = seed_num ;
     chain_num += 1;
    // if (chain_num >= 100){
    //   *((int*)(&buffer[chain_num_addr])) = chain_num ;
    //   *((int*)(&buffer[idx_end_addr])) = buffer_idx/4 ;
    //   chain_num = 0;
    //   read_num += 1;
    //   idx_end_addr = buffer_idx;
    //   buffer_idx += 4;
    //   // pack the read sequence
    //   *((int*)(&buffer[buffer_idx])) = seq_->l_seq;
    //   buffer_idx += 4;
    //   for ( int i = 0; i < seq_->l_seq; i++ ){
    //     counter8 += 1;
    //     tmp_int = tmp_int << 4 | seq_->seq[i];
    //     if ( counter8 % 8 ==0 ) {
    //       *(uint32_t*) (&buffer[buffer_idx]) = tmp_int ;
    //       buffer_idx += 4 ;
    //     }
    //   }
    //   if ( counter8 % 8 !=0 ) {
    //     tmp_int = tmp_int << (4*(8 - counter8 % 8));
    //     *(uint32_t*) (&buffer[buffer_idx]) = tmp_int ;
    //     buffer_idx += 4 ;
    //   }
    //   counter8 = 0;
    //   tmp_int = 0;
    //   chain_num_addr = buffer_idx;
    //   buffer_idx += 4;
    // }
  }
  *((int*)(&buffer[chain_num_addr])) = chain_num ;
  *((int*)(&buffer[idx_end_addr])) = buffer_idx/4 ;

}

void SWRead::nextTask(ExtParam** &task_batch,int &task_num) {
  if (is_finished_) {
    return ;
  }
  for ( ; chain_idx_ < chains_->n; 
          seed_idx_ > 0 ? seed_idx_-- : seed_idx_ = chains_->a[++chain_idx_].n-1)
  {
      uint32_t sorted_idx = (uint32_t)(ref_[chain_idx_].srt[seed_idx_]);

    // get next available seed in the current read
      mem_seed_t* seed_array = &chains_->a[chain_idx_].seeds[sorted_idx];
      // initialize the newreg
      mem_alnreg_t* newreg = kv_pushp(mem_alnreg_t, *alnregs_);
      memset(newreg, 0, sizeof(mem_alnreg_t));

      // Push chain_idx to added newreg
      chain_idxes_->push_back(chain_idx_);

      newreg->score  = seed_array->len * aux_->opt->a;
      newreg->truesc = seed_array->len * aux_->opt->a;
      newreg->qb = 0;
      newreg->rb = seed_array->rbeg;
      newreg->qe = seq_->l_seq;
      newreg->re = seed_array->rbeg + seed_array->len; 
      newreg->rid = chains_->a[chain_idx_].rid;
      newreg->seedlen0 = seed_array->len; 
      newreg->frac_rep = chains_->a[chain_idx_].frac_rep;
      newreg->w = aux_->opt->w;
   }
  int region_num = 0;
  int read_task_num = 0;
  chain_idx_ = 0;
  seed_idx_ = chains_->a[0].n-1;
 //TODO 
  for ( ; chain_idx_ < chains_->n; 
          seed_idx_ > 0 ? seed_idx_-- : seed_idx_ = chains_->a[++chain_idx_].n-1)
  {
      uint32_t sorted_idx = (uint32_t)(ref_[chain_idx_].srt[seed_idx_]);

      // get next available seed in the current read
      mem_seed_t* seed_array = &chains_->a[chain_idx_].seeds[sorted_idx];

      if (seed_array->qbeg > 0 ||
          seed_array->qbeg + seed_array->len != seq_->l_seq)
      { 
        // need to do extension
        task_batch[task_num + read_task_num] = getTask(
            aux_->opt,
            seed_array,
            (const uint8_t*)seq_->seq,
            seq_->l_seq,
            ref_[chain_idx_].rmax[0],
            ref_[chain_idx_].rmax[1],
            ref_[chain_idx_].rseq,
            read_idx_);

        task_batch[task_num + read_task_num]->read_obj = this;
        task_batch[task_num + read_task_num]->newreg = &(alnregs_->a[region_num]);
        task_batch[task_num + read_task_num]->chain_idx = chain_idx_;
        task_batch[task_num + read_task_num]->chain = &chains_->a[chain_idx_];
        region_num++;
        read_task_num++;
      }
      else {
        region_num++;
      }
  }
  task_num += read_task_num;
  return ;
}

inline ExtParam* SWRead::getTask(
    mem_opt_t *opt,
    const mem_seed_t *seed, 
    const uint8_t *query,
    int l_query,
    int64_t rmax_0, 
    int64_t rmax_1,
    uint8_t *rseq,
    int idx)
{
  uint64_t start_ts = getNs();
  ExtParam *SwTask = new ExtParam;
  int i = 0;
  SwTask->leftQlen = seed->qbeg;
  if(SwTask->leftQlen > 0)
  {
    SwTask->leftQs = new uint8_t[SwTask->leftQlen];
    for(i = 0;i < SwTask->leftQlen; i++)
      SwTask->leftQs[i] = query[SwTask->leftQlen-1-i];
    SwTask->leftRlen =(int)( seed->rbeg - rmax_0) ;
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
    SwTask->rightQs =(uint8_t *) query +qe;
    int64_t re = seed->rbeg + seed->len - rmax_0;
    SwTask->rightRlen =(int) (rmax_1 - rmax_0 -re) ;
    SwTask->rightRs = (uint8_t *)rseq +re;
  }
  else
  {
    SwTask->rightQs = NULL;
    SwTask->rightRlen = 0;
    SwTask->rightRs = NULL;
  }

  SwTask->h0 = seed->len*opt->a ;
  SwTask->regScore = seed->len*opt->a ;
  SwTask->qBeg = seed->qbeg ;
  SwTask->rBeg = seed->rbeg ;     // for testing
  SwTask->seedLength = seed->len ;
  SwTask->idx = idx ;
  SwTask->l_query = l_query ;
  get_task_time += ( getNs() - start_ts );
  return SwTask;
}

inline int SWRead::testExtension(
    mem_opt_t *opt,
    mem_seed_t& seed,
    mem_alnreg_v& alnregv, 
    int l_query) 
{
  long rdist = -1;
  int qdist = -1;  
  int maxgap = -1;
  int mindist = -1;
  int w = -1;
  int breakidx = 0;
  bool isbreak = false;

  int i;
  for (i = 0; i < alnregv.n; i++) {
    if (seed.rbeg < alnregv.a[i].rb ||
        seed.rbeg + seed.len> alnregv.a[i].re||
        seed.qbeg < alnregv.a[i].qb ||
        seed.qbeg + seed.len > alnregv.a[i].qe)
    {
      continue; 
    }
    if (seed.len - alnregv.a[i].seedlen0 > .1 * l_query){
      continue;
    }
    qdist   = seed.qbeg - alnregv.a[i].qb;
    rdist   = seed.rbeg - alnregv.a[i].rb;
    mindist = (qdist < rdist) ? qdist : (int)rdist;
    maxgap  = cal_max_gap(opt,mindist);
    w = (maxgap < alnregv.a[i].w) ? maxgap : alnregv.a[i].w;

    if (qdist - rdist < w &&
        rdist - qdist < w) 
    {
      break;
    }

    qdist   = alnregv.a[i].qe -(seed.qbeg + seed.len);
    rdist   = alnregv.a[i].re - (seed.rbeg + seed.len);
    mindist = (qdist < rdist) ? qdist : (int)rdist;
    maxgap  = cal_max_gap(opt, mindist);
    w = (maxgap < alnregv.a[i].w) ? maxgap : alnregv.a[i].w;
    if (qdist - rdist < w &&
        rdist - qdist < w) 
    {
      break;
    }
  }
  return i;
}

inline int SWRead::checkOverlap(
    int startidx,
    mem_seed_t& seed,
    mem_chain_t& chain,
    uint64_t *srt)
{
  int i;
  for (i = startidx; i < chain.n; i++) {
    const mem_seed_t* targetseed;
    if(srt[i]==0) {
      continue;
    }
    targetseed = &chain.seeds[(uint32_t)srt[i]];
    if (targetseed->len < seed.len* 0.95) {
      continue;
    }
    if (seed.qbeg <= targetseed->qbeg && 
        seed.qbeg + seed.len - targetseed->qbeg >= seed.len>>2 &&
        targetseed->qbeg - seed.qbeg != targetseed->rbeg-seed.rbeg) {
      break;
    }
    if (targetseed->qbeg <= seed.qbeg &&
        targetseed->qbeg + targetseed->len - seed.qbeg >= seed.len>>2 &&
        seed.qbeg - targetseed->qbeg != seed.rbeg - targetseed->rbeg) 
    {
      break;
    }
  }
  return i;
}

