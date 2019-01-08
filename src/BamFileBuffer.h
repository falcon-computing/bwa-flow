#ifndef BAMFILEBUFFER_H
#define BAMFILEBUFFER_H

#include <iostream>
#include <string.h>

//#include "hts_internal.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "htslib/hfile.h"

/* Wrapper class for hFile and bam_write1
 */
class BamFileBuffer {
 public:
  BamFileBuffer(int is_be, int size = 4*1024*1024):
    size_(size), pos_(0), is_be_(is_be)
  {
    data_ = (char*)malloc(size);
  }

  ~BamFileBuffer() {
    free(data_);
  }

  size_t write(bam1_t* b) {
    const bam1_core_t *c = &b->core;
    uint32_t x[8], block_len = b->l_data + 32, y;
    int i, ok;
    x[0] = c->tid;
    x[1] = c->pos;
    x[2] = (uint32_t)c->bin<<16 | c->qual<<8 | c->l_qname;
    x[3] = (uint32_t)c->flag<<16 | c->n_cigar;
    x[4] = c->l_qseq;
    x[5] = c->mtid;
    x[6] = c->mpos;
    x[7] = c->isize;
    ok = check_space(4 + block_len);
    if (is_be_) {
        for (i = 0; i < 8; ++i) ed_swap_4p(x + i);
        y = block_len;
        if (ok) ok = (write_buf((char*)ed_swap_4p(&y), 4) >= 0);
        swap_data(c, b->l_data, b->data, 1);
    } else {
        if (ok) ok = (write_buf((char*)&block_len, 4) >= 0);
    }
    if (ok) ok = (write_buf((char*)x, 32) >= 0);
    if (ok) ok = (write_buf((char*)b->data, b->l_data) >= 0);
    if (is_be_) swap_data(c, b->l_data, b->data, 0);
    return ok? 4 + block_len : -1;
  }

  const char* get_data() const {
    return data_;
  }

  size_t get_size() const {
    return pos_;
  }

 private:

  int check_space(int nbytes) {
    while (pos_ + nbytes >= size_) {
      size_ *= 2;
      data_ = (char*)realloc(data_, size_);
      if (!data_) {
        LOG(ERROR) << "cannot allocate memory for BAM buffer (" 
          << size_ << " bytes)";
        return 0;
      }
    }
    return 1;
  }

  size_t write_buf(char* in, size_t nbytes) {
    memcpy(data_ + pos_, in, nbytes);
    pos_ += nbytes;

    return nbytes;
  }

  static void swap_data(const bam1_core_t *c, int l_data, uint8_t *data, int is_host)
  {
    uint8_t *s;
    uint32_t *cigar = (uint32_t*)(data + c->l_qname);
    uint32_t i, n;
    s = data + c->n_cigar*4 + c->l_qname + c->l_qseq + (c->l_qseq + 1)/2;
    for (i = 0; i < c->n_cigar; ++i) ed_swap_4p(&cigar[i]);
    while (s < data + l_data) {
        int size;
        s += 2; // skip key
        size = aux_type2size(*s); ++s; // skip type
        switch (size) {
        case 1: ++s; break;
        case 2: ed_swap_2p(s); s += 2; break;
        case 4: ed_swap_4p(s); s += 4; break;
        case 8: ed_swap_8p(s); s += 8; break;
        case 'Z':
        case 'H':
            while (*s) ++s;
            ++s;
            break;
        case 'B':
            size = aux_type2size(*s); ++s;
            if (is_host) memcpy(&n, s, 4), ed_swap_4p(s);
            else ed_swap_4p(s), memcpy(&n, s, 4);
            s += 4;
            switch (size) {
            case 1: s += n; break;
            case 2: for (i = 0; i < n; ++i, s += 2) ed_swap_2p(s); break;
            case 4: for (i = 0; i < n; ++i, s += 4) ed_swap_4p(s); break;
            case 8: for (i = 0; i < n; ++i, s += 8) ed_swap_8p(s); break;
            }
            break;
        }
    }
  }

  static inline int aux_type2size(uint8_t type)
  {
    switch (type) {
      case 'A': case 'c': case 'C':
        return 1;
      case 's': case 'S':
        return 2;
      case 'i': case 'I': case 'f':
        return 4;
      case 'd':
        return 8;
      case 'Z': case 'H': case 'B':
        return type;
      default:
        return 0;
    }
  }

  size_t size_;
  size_t pos_;
  int    is_be_;
  char*  data_;
};

#endif
