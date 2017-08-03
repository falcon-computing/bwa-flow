#include <stdio.h>
#include <stdlib.h>
#include <sstream>

#include "bwa/bwamem.h"
#include "bwa/bntseq.h"
#include "bwa/utils.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

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
static inline void getStr(std::stringstream &ss, char* &dst) {
  size_t length = 0;
  getT(ss, length);

  if (length > 0) {
    dst = (char*)malloc(length+1);
    ss.read(dst, length);
    dst[length] = '\0';
  }
}

void serialize(std::stringstream &ss, bseq1_t& seq) {
  putStr(ss, seq.name);
  putStr(ss, seq.comment);
  putStr(ss, seq.seq);
  putStr(ss, seq.qual);
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
