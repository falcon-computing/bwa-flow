#include <stdio.h>
#include <stdlib.h>
#include <sstream>

#include "bwa/bwamem.h"
#include "bwa/bntseq.h"
#include "bwa/utils.h"

#include "util.h"
#include "allocation_wrapper.h"

template <typename T>
void assertEQ(int i, int j, const char* msg, T const &a, T const &b) {
  if (a != b) {
    std::stringstream ss;
    ss << "[" << i << "," << j << "] " << msg << ":"
       << a << "!=" << b << "\n";
    throw resultsError(ss.str().c_str());
  }
}

// Encode a scalar value to serialized data
template <typename T>
static inline void putT(std::stringstream &ss, T value) {
  ss.write(reinterpret_cast<char*>(&value), sizeof(value));
}

// Decode a scalar value from serialized data
template <typename T>
static inline void getT(std::stringstream &ss, T &value) {
  ss.read(reinterpret_cast<char*>(&value), sizeof(value));
}

// Store a string with its length to serialized data
static inline void putStr(std::stringstream &ss, const char* str) {
  if (str) {
    size_t length = strlen(str);
    putT(ss, length);
    ss.write(str, length);
  }
  else {
    size_t length = 0;
    putT(ss, length);
  }
}

// Retrieve a string from serialized data
static inline size_t getStr(std::stringstream &ss, char* &dst) {
  size_t length = 0;
  getT(ss, length);

  if (length > 0) {
    dst = (char*)malloc(length+1);
    ss.read(dst, length);
    dst[length] = '\0';
  }
  return length;
}

// for chain_record ser
void serialize(std::stringstream &ss, bseq1_t& seq) {
  putStr(ss, seq.name);
  putStr(ss, seq.comment);

  // special treatment because seq.seq may use 2-bit encoding
  putT(ss, seq.l_seq);
  ss.write(seq.seq, seq.l_seq);

  putStr(ss, seq.qual);
}

void deserialize(std::stringstream &ss, bseq1_t& seq) {
  if (!getStr(ss, seq.name)) seq.name = 0;
  if (!getStr(ss, seq.comment)) seq.comment = 0;

  // special treatment because seq.seq may use 2-bit encoding
  getT(ss, seq.l_seq);
  seq.seq = (char*)malloc(seq.l_seq);
  ss.read(seq.seq, seq.l_seq);

  if (!getStr(ss, seq.qual)) seq.qual = 0;
}

void serialize(std::stringstream &ss, mem_chain_v& chains) {
  putT(ss, chains.n);
  putT(ss, chains.m);
  for (int i = 0; i < chains.n; i++) {
    serialize(ss, chains.a[i]);
  }
}

void deserialize(std::stringstream &ss, mem_chain_v& chains) {
  getT(ss, chains.n);
  getT(ss, chains.m);
  chains.a = (mem_chain_t*)malloc(chains.m*sizeof(mem_chain_t));
  for (int i = 0; i < chains.n; i++) {
     deserialize(ss, chains.a[i]);
  }
}

typedef struct {
  uint32_t w:29, kept:2, is_alt:1;
} chain_opts_t;

typedef union {
  chain_opts_t opts;
  uint32_t     value;
} chain_opts_u;

typedef struct {
	int n_comp:30, is_alt:2; // number of sub-alignments chained together
} alnreg_opts_t;

typedef union {
  alnreg_opts_t opts;
  int           value;
} alnreg_opts_u;

void serialize(std::stringstream &ss, mem_chain_t& chain) {
  putT(ss, chain.n);
  putT(ss, chain.m);
  putT(ss, chain.first);
  putT(ss, chain.rid);

  // take care of the bitfield in original struct
  chain_opts_u v;
  v.opts.w = chain.w;
  v.opts.kept = chain.kept;
  v.opts.is_alt = chain.is_alt;

  putT(ss, v.value);

  putT(ss, chain.frac_rep);
  putT(ss, chain.pos);
  for (int i = 0; i < chain.n; i++) {
    //serialize(ss, &chain->seeds[i]);
    putT(ss, chain.seeds[i]);
  }
}

void deserialize(std::stringstream &ss, mem_chain_t& chain) {

  getT(ss, chain.n);
  getT(ss, chain.m);
  getT(ss, chain.first);
  getT(ss, chain.rid);

  // here need to convert uint32_t to struct
  chain_opts_u v;
  getT(ss, v.value);

  chain.w = v.opts.w;
  chain.kept = v.opts.kept;
  chain.is_alt = v.opts.is_alt;

  getT(ss, chain.frac_rep);
  getT(ss, chain.pos);

  chain.seeds = (mem_seed_t*)malloc(chain.m*sizeof(mem_seed_t));
  for (int i = 0; i < chain.n; i++) {
    //serialize(ss, &chain->seeds[i]);
    getT(ss, chain.seeds[i]);
  }
}

void serialize(std::stringstream &ss, mem_seed_t& seed) {
  putT(ss, seed.rbeg);
  putT(ss, seed.qbeg);
  putT(ss, seed.len);
  putT(ss, seed.score);
}

void serialize(std::stringstream &ss, mem_alnreg_t& alnreg) {
  putT(ss, alnreg.rb); 
  putT(ss, alnreg.re); 
  putT(ss, alnreg.qb); 
  putT(ss, alnreg.qe); 
	putT(ss, alnreg.rid);
	putT(ss, alnreg.score);
	putT(ss, alnreg.truesc);
	putT(ss, alnreg.sub);
	putT(ss, alnreg.alt_sc);
	putT(ss, alnreg.csub);
	putT(ss, alnreg.sub_n);
	putT(ss, alnreg.w);
	putT(ss, alnreg.seedcov);
	putT(ss, alnreg.secondary);
	putT(ss, alnreg.secondary_all);
	putT(ss, alnreg.seedlen0);

  alnreg_opts_u v;
  v.opts.n_comp = alnreg.n_comp;
  v.opts.is_alt = alnreg.is_alt;
  putT(ss, v.value);

	putT(ss, alnreg.frac_rep);
	putT(ss, alnreg.hash);
}

void deserialize(std::stringstream &ss, mem_alnreg_t& alnreg) {
  getT(ss, alnreg.rb); 
  getT(ss, alnreg.re); 
  getT(ss, alnreg.qb); 
  getT(ss, alnreg.qe); 
	getT(ss, alnreg.rid);
	getT(ss, alnreg.score);
	getT(ss, alnreg.truesc);
	getT(ss, alnreg.sub);
	getT(ss, alnreg.alt_sc);
	getT(ss, alnreg.csub);
	getT(ss, alnreg.sub_n);
	getT(ss, alnreg.w);
	getT(ss, alnreg.seedcov);
	getT(ss, alnreg.secondary);
	getT(ss, alnreg.secondary_all);
	getT(ss, alnreg.seedlen0);

  alnreg_opts_u v;
  getT(ss, v.value);

  alnreg.n_comp = v.opts.n_comp;
  alnreg.is_alt = v.opts.is_alt;

	getT(ss, alnreg.frac_rep);
	getT(ss, alnreg.hash);
}

void serialize(std::stringstream &ss, mem_alnreg_v& alnregs) {
  putT(ss, alnregs.n);
  putT(ss, alnregs.m);
  for (int i = 0; i < alnregs.n; i++) {
    putT(ss, alnregs.a[i]);
  }
}

void deserialize(std::stringstream &ss, mem_alnreg_v& alnregs) {
  getT(ss, alnregs.n);
  getT(ss, alnregs.m);
  alnregs.a = (mem_alnreg_t*)malloc(alnregs.m*sizeof(mem_alnreg_t));
  for (int i = 0; i < alnregs.n; i++) {
    getT(ss, alnregs.a[i]);
  }
}

std::string serialize(SeqsRecord& record) {

  uint64_t start_idx = record.start_idx;
  int      batch_num = record.batch_num;

  std::stringstream ss;

  putT(ss, start_idx);
  putT(ss, batch_num);

  for (int i = 0; i < batch_num; i++) {
    serialize(ss, record.seqs[i]);
  }

  return ss.str();
}

void deserialize(const char* data, size_t length, SeqsRecord &output) {

  uint64_t start_idx = 0;
  int      batch_num = 0;

  std::stringstream ss;
  ss.write(data, length);

  getT(ss, start_idx);
  getT(ss, batch_num);

  bseq1_t* seqs = (bseq1_t*)calloc(batch_num, sizeof(bseq1_t));

  for (int i = 0; i < batch_num; i++) {
    deserialize(ss, seqs[i]);
  }

  output.start_idx = start_idx;
  output.batch_num = batch_num;
  output.seqs = seqs;
}

std::string serialize(ChainsRecord& record) {
  uint64_t start_idx = record.start_idx;
  int      batch_num = record.batch_num;

  std::stringstream ss;

  putT(ss, start_idx);
  putT(ss, batch_num);

  for (int i = 0; i < batch_num; i++) {
    serialize(ss, record.seqs[i]);
    serialize(ss, record.chains[i]);
  }

  return ss.str();
}

void deserialize(const char* data, size_t length, ChainsRecord& output) {
  uint64_t start_idx = 0;
  int      batch_num = 0;

  std::stringstream ss;
  ss.write(data, length);

  getT(ss, start_idx);
  getT(ss, batch_num);

  bseq1_t*     seqs   = (bseq1_t*)calloc(batch_num, sizeof(bseq1_t));
  mem_chain_v* chains = (mem_chain_v*)calloc(batch_num, sizeof(mem_chain_v));

  for (int i = 0; i < batch_num; i++) {
    deserialize(ss, seqs[i]);
    deserialize(ss, chains[i]);
  }

  output.start_idx = start_idx;
  output.batch_num = batch_num;
  output.seqs = seqs;
  output.chains = chains;
}

std::string serialize(RegionsRecord& record) {
  uint64_t start_idx = record.start_idx;
  int      batch_num = record.batch_num;

  std::stringstream ss;

  putT(ss, start_idx);
  putT(ss, batch_num);

  for (int i = 0; i < batch_num; i++) {
    serialize(ss, record.seqs[i]);
    serialize(ss, record.alnreg[i]);
  }

  return ss.str();
}

void deserialize(const char* data, size_t length, RegionsRecord& output) {
  uint64_t start_idx = 0;
  int      batch_num = 0;

  std::stringstream ss;
  ss.write(data, length);

  getT(ss, start_idx);
  getT(ss, batch_num);

  bseq1_t*     seqs     = (bseq1_t*)calloc(batch_num, sizeof(bseq1_t));
  mem_alnreg_v* alnregs = (mem_alnreg_v*)calloc(batch_num, sizeof(mem_alnreg_v));

  for (int i = 0; i < batch_num; i++) {
    deserialize(ss, seqs[i]);
    deserialize(ss, alnregs[i]);
  }

  output.start_idx = start_idx;
  output.batch_num = batch_num;
  output.seqs = seqs;
  output.alnreg = alnregs;
}
