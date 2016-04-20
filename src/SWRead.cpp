#include "SWRead.h"

SWRead::SWRead(int idx,
    ktp_aux_t* aux,
    const bseq1_t* seq, 
    mem_chainref_t* ref,
    const mem_chain_v* chains,
    mem_alnreg_v* alnregs):
  read_idx_(idx),
  aux_(aux), seq_(seq), ref_(ref), 
  chains_(chains), alnregs_(alnregs),
  is_pend_(false), chain_idx_(0)
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

void SWRead::finish() {
  is_pend_ = false;
}

enum SWRead::TaskStatus SWRead::nextTask(ExtParam* &task, mem_alnreg_t* newreg) {

  if (is_pend_) {
    return TaskStatus::Pending;
  }
  else if (is_finished_) {
    return TaskStatus::Finished;
  }

  for ( ; chain_idx_ < chains_->n; 
          seed_idx_ > 0 ? seed_idx_-- : seed_idx_ = chains_->a[++chain_idx_].n-1)
  {
    uint32_t sorted_idx = (uint32_t)(ref_[chain_idx_].srt[seed_idx_]);

    // get next available seed in the current read
    mem_seed_t* seed_array = &chains_->a[chain_idx_].seeds[sorted_idx];

    if (alnregs_->n > testExtension(
            aux_->opt, 
            *seed_array,
            *alnregs_, 
            seq_->l_seq) 
        && 
        chains_->a[chain_idx_].n == checkOverlap(
            seed_idx_+1,
            *seed_array,
            chains_->a[chain_idx_],
            ref_[chain_idx_].srt)) 
    {                                               
      // Mark this seed as uncomputed
      ref_[chain_idx_].srt[seed_idx_] = 0;         
    }
    else {
      // initialize the newreg
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

      if (seed_array->qbeg > 0 ||
          seed_array->qbeg + seed_array->len != seq_->l_seq)
      { 
        // need to do extension
        task = getTask(
            aux_->opt,
            seed_array,
            (const uint8_t*)seq_->seq,
            seq_->l_seq,
            ref_[chain_idx_].rmax[0],
            ref_[chain_idx_].rmax[1],
            ref_[chain_idx_].rseq,
            read_idx_);

        task->read_obj = this;
        task->newreg = newreg;
        task->chain_idx = chain_idx_;
        task->chain = &chains_->a[chain_idx_];

        //fprintf(stderr, "SMRead got task for %d: %d, %d\n", read_idx_, chain_idx_, seed_idx_);

        is_pend_ = true;
        // update chain and seed indexes before return
        if (seed_idx_ > 0) {
          seed_idx_ --;
        }
        else {
          seed_idx_ = chains_->a[++chain_idx_].n - 1;
        }

        return TaskStatus::Successful;
      }
      else {
        // no need to extend, just push to alnreg_v 
        newreg->seedcov=0;

        for (int j = 0; j < chains_->a[chain_idx_].n; j++) {
          const mem_seed_t *t = &chains_->a[chain_idx_].seeds[j];
          if (t->qbeg >= newreg->qb && 
              t->qbeg + t->len <= newreg->qe && 
              t->rbeg >= newreg->rb && 
              t->rbeg + t->len <= newreg->re) {
            newreg->seedcov += t->len; 
          }
        }
        kv_push(mem_alnreg_t, *alnregs_, *newreg);
      }
    }
  }
  is_finished_ = true;
  return TaskStatus::Finished;
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
  if (seed->qbeg > 0 ||
      seed->qbeg + seed->len != l_query) {
    
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
      SwTask->rightQs   = new uint8_t[SwTask->rightQlen];
      for(int i = 0;i<SwTask->rightQlen;i++)
        SwTask->rightQs[i] = query[i+qe];
      int64_t re = seed->rbeg + seed->len - rmax_0;
      SwTask->rightRlen =(int) (rmax_1 - rmax_0 -re) ;
      SwTask->rightRs = new uint8_t[SwTask->rightRlen];
      for(int i = 0;i<SwTask->rightRlen; i++)
        SwTask->rightRs[i] = rseq[i+(int)re];
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
    SwTask->l_query = l_query ;
    return SwTask;
  }
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

inline int SWRead::cal_max_gap(
    const mem_opt_t *opt, 
    int qlen
) {
	int l_del = (int)((double)(qlen * opt->a - opt->o_del) / opt->e_del + 1.);
	int l_ins = (int)((double)(qlen * opt->a - opt->o_ins) / opt->e_ins + 1.);
	int l = l_del > l_ins? l_del : l_ins;
	l = l > 1? l : 1;
	return l < opt->w<<1? l : opt->w<<1;
}


