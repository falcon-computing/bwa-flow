#include <ctype.h>
#include <glog/logging.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

#include "bwa/bwa.h"
#include "bwa/bwamem.h"
#include "bwa/kvec.h"
#include "bwa/utils.h"
#include "bwa/bntseq.h"
#include "bwa/kseq.h"
#include "bwa_wrapper.h"
#include "config.h"
#include "allocation_wrapper.h"

extern unsigned char nst_nt4_table[256];

extern gzFile fp_idx, fp2_read2 ;
extern void *ko_read1 , *ko_read2 ;

int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: bwa (alignment via Burrows-Wheeler transformation)\n");
	fprintf(stderr, "Version: %s\n", VERSION);
	fprintf(stderr, "Contact: Heng Li <lh3@sanger.ac.uk>\n\n");
	fprintf(stderr, "Usage:   bwa <command> [options]\n\n");
	fprintf(stderr, "Command: index         index sequences in the FASTA format\n");
	fprintf(stderr, "         mem           BWA-MEM algorithm\n");
	fprintf(stderr, "         fastmap       identify super-maximal exact matches\n");
	fprintf(stderr, "         pemerge       merge overlapping paired ends (EXPERIMENTAL)\n");
	fprintf(stderr, "         aln           gapped/ungapped alignment\n");
	fprintf(stderr, "         samse         generate alignment (single ended)\n");
	fprintf(stderr, "         sampe         generate alignment (paired ended)\n");
	fprintf(stderr, "         bwasw         BWA-SW for long queries\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "         shm           manage indices in shared memory\n");
	fprintf(stderr, "         fa2pac        convert FASTA to PAC format\n");
	fprintf(stderr, "         pac2bwt       generate BWT from PAC\n");
	fprintf(stderr, "         pac2bwtgen    alternative algorithm for generating BWT\n");
	fprintf(stderr, "         bwtupdate     update .bwt to the new format\n");
	fprintf(stderr, "         bwt2sa        generate SA from BWT and Occ\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Note: To use BWA, you need to first index the genome with `bwa index'.\n"
                  "      There are three alignment algorithms in BWA: `mem', `bwasw', and\n"
                  "      `aln/samse/sampe'. If you are not sure which to use, try `bwa mem'\n"
                  "      first. Please `man ./bwa.1' for the manual.\n\n");
	return 1;
}

static void update_a(mem_opt_t *opt, const mem_opt_t *opt0) {
	if (opt0->a) { // matching score is changed
		if (!opt0->b) opt->b *= opt->a;
		if (!opt0->T) opt->T *= opt->a;
		if (!opt0->o_del) opt->o_del *= opt->a;
		if (!opt0->e_del) opt->e_del *= opt->a;
		if (!opt0->o_ins) opt->o_ins *= opt->a;
		if (!opt0->e_ins) opt->e_ins *= opt->a;
		if (!opt0->zdrop) opt->zdrop *= opt->a;
		if (!opt0->pen_clip5) opt->pen_clip5 *= opt->a;
		if (!opt0->pen_clip3) opt->pen_clip3 *= opt->a;
		if (!opt0->pen_unpaired) opt->pen_unpaired *= opt->a;
	}
}

int pack_bwa_mem_args(std::vector<const char*> & bwa_mem_args) {
  if (FLAGS_k) bwa_mem_args.push_back("-k"), bwa_mem_args.push_back(std::to_string(FLAGS_k).c_str());
  if (FLAGS_1) bwa_mem_args.push_back("-1");
  if (!FLAGS_x.empty()) bwa_mem_args.push_back("-x"), bwa_mem_args.push_back(FLAGS_x.c_str());
  if (FLAGS_w) bwa_mem_args.push_back("-w"), bwa_mem_args.push_back(std::to_string(FLAGS_w).c_str());
  if (FLAGS_A) bwa_mem_args.push_back("-A"), bwa_mem_args.push_back(std::to_string(FLAGS_A).c_str());
  if (FLAGS_B) bwa_mem_args.push_back("-B"), bwa_mem_args.push_back(std::to_string(FLAGS_B).c_str());
  if (FLAGS_T) bwa_mem_args.push_back("-T"), bwa_mem_args.push_back(std::to_string(FLAGS_T).c_str());
  if (FLAGS_U) bwa_mem_args.push_back("-U"), bwa_mem_args.push_back(std::to_string(FLAGS_U).c_str());
  if (FLAGS_P) bwa_mem_args.push_back("-P");
  if (FLAGS_a) bwa_mem_args.push_back("-a");
  if (FLAGS_p) bwa_mem_args.push_back("-p");
  if (FLAGS_S) bwa_mem_args.push_back("-S");
  if (FLAGS_Y) bwa_mem_args.push_back("-Y");
  if (FLAGS_V) bwa_mem_args.push_back("-V");
  if (FLAGS_c) bwa_mem_args.push_back("-c"), bwa_mem_args.push_back(std::to_string(FLAGS_c).c_str()); 
  if (FLAGS_d) bwa_mem_args.push_back("-d"), bwa_mem_args.push_back(std::to_string(FLAGS_d).c_str()); 
  if (FLAGS_j) bwa_mem_args.push_back("-j");
  if (FLAGS_r) bwa_mem_args.push_back("-r"), bwa_mem_args.push_back(std::to_string(FLAGS_r).c_str());
  if (FLAGS_D) bwa_mem_args.push_back("-D"), bwa_mem_args.push_back(std::to_string(FLAGS_D).c_str());
  if (FLAGS_m) bwa_mem_args.push_back("-m"), bwa_mem_args.push_back(std::to_string(FLAGS_m).c_str());
  if (FLAGS_s) bwa_mem_args.push_back("-s"), bwa_mem_args.push_back(std::to_string(FLAGS_s).c_str());
  if (FLAGS_G) bwa_mem_args.push_back("-G"), bwa_mem_args.push_back(std::to_string(FLAGS_G).c_str());
  if (FLAGS_N) bwa_mem_args.push_back("-N"), bwa_mem_args.push_back(std::to_string(FLAGS_N).c_str());
  if (FLAGS_W) bwa_mem_args.push_back("-W"), bwa_mem_args.push_back(std::to_string(FLAGS_W).c_str());
  if (FLAGS_y) bwa_mem_args.push_back("-y"), bwa_mem_args.push_back(std::to_string(FLAGS_y).c_str());
  if (FLAGS_C) bwa_mem_args.push_back("-C");
  if (FLAGS_K) bwa_mem_args.push_back("-K"), bwa_mem_args.push_back(std::to_string(FLAGS_K).c_str());
  if (FLAGS_X) bwa_mem_args.push_back("-X"), bwa_mem_args.push_back(std::to_string(FLAGS_X).c_str());
#ifdef USE_HTSLIB
  if (FLAGS_o) bwa_mem_args.push_back("-o"), bwa_mem_args.push_back(std::to_string(FLAGS_o).c_str());
#endif
  if (!FLAGS_h.empty()) bwa_mem_args.push_back("-h"), bwa_mem_args.push_back(FLAGS_h.c_str());
  if (FLAGS_Q) bwa_mem_args.push_back("-Q"), bwa_mem_args.push_back(std::to_string(FLAGS_Q).c_str());
  if (!FLAGS_O.empty()) bwa_mem_args.push_back("-O"), bwa_mem_args.push_back(FLAGS_O.c_str());
  if (!FLAGS_E.empty()) bwa_mem_args.push_back("-E"), bwa_mem_args.push_back(FLAGS_E.c_str());
  if (!FLAGS_L.empty()) bwa_mem_args.push_back("-L"), bwa_mem_args.push_back(FLAGS_L.c_str());
  if (!FLAGS_H.empty()) bwa_mem_args.push_back("-H"), bwa_mem_args.push_back(FLAGS_H.c_str());
  if (!FLAGS_I.empty()) bwa_mem_args.push_back("-I"), bwa_mem_args.push_back(FLAGS_I.c_str());

  return 0;
}

int pre_process(int argc,
    char *argv[],
    ktp_aux_t* aux,
    bool is_master
) {
	mem_opt_t *opt, opt0;
	int fd, fd2, i, c, ignore_alt = 0, no_mt_io = 0;
	int fixed_chunk_size = -1;
	char *p, *rg_line = 0, *hdr_line = 0;
	const char *mode = 0;
	mem_pestat_t pes[4];
	memset(pes, 0, 4 * sizeof(mem_pestat_t));
	for (i = 0; i < 4; ++i) pes[i].failed = 1;

	aux->opt = opt = mem_opt_init();
	memset(&opt0, 0, sizeof(mem_opt_t));
#ifdef USE_HTSLIB
  while ((c = getopt(argc, argv, "1paMCSPVYjk:c:v:s:r:t:R:A:B:O:E:U:w:L:d:T:Q:D:m:I:N:W:x:G:h:y:K:X:H:o:")) >= 0) {
#else
    while ((c = getopt(argc, argv, "1paMCSPVYjk:c:v:s:r:t:R:A:B:O:E:U:w:L:d:T:Q:D:m:I:N:W:x:G:h:y:K:X:H:")) >= 0) {
#endif
      if (c == 'k') opt->min_seed_len = atoi(optarg), opt0.min_seed_len = 1;
      else if (c == '1') no_mt_io = 1;
      else if (c == 'x') mode = optarg;
      else if (c == 'w') opt->w = atoi(optarg), opt0.w = 1;
      else if (c == 'A') opt->a = atoi(optarg), opt0.a = 1;
      else if (c == 'B') opt->b = atoi(optarg), opt0.b = 1;
      else if (c == 'T') opt->T = atoi(optarg), opt0.T = 1;
      else if (c == 'U') opt->pen_unpaired = atoi(optarg), opt0.pen_unpaired = 1;
      else if (c == 't') opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1;
      else if (c == 'P') opt->flag |= MEM_F_NOPAIRING;
      else if (c == 'a') opt->flag |= MEM_F_ALL;
      else if (c == 'p') opt->flag |= MEM_F_PE | MEM_F_SMARTPE;
      else if (c == 'M') opt->flag |= MEM_F_NO_MULTI;
      else if (c == 'S') opt->flag |= MEM_F_NO_RESCUE;
      else if (c == 'Y') opt->flag |= MEM_F_SOFTCLIP;
      else if (c == 'V') opt->flag |= MEM_F_REF_HDR;
      else if (c == 'c') opt->max_occ = atoi(optarg), opt0.max_occ = 1;
      else if (c == 'd') opt->zdrop = atoi(optarg), opt0.zdrop = 1;
      else if (c == 'v') bwa_verbose = atoi(optarg);
      else if (c == 'j') ignore_alt = 1;
      else if (c == 'r') opt->split_factor = atof(optarg), opt0.split_factor = 1.;
      else if (c == 'D') opt->drop_ratio = atof(optarg), opt0.drop_ratio = 1.;
      else if (c == 'm') opt->max_matesw = atoi(optarg), opt0.max_matesw = 1;
      else if (c == 's') opt->split_width = atoi(optarg), opt0.split_width = 1;
      else if (c == 'G') opt->max_chain_gap = atoi(optarg), opt0.max_chain_gap = 1;
      else if (c == 'N') opt->max_chain_extend = atoi(optarg), opt0.max_chain_extend = 1;
      else if (c == 'W') opt->min_chain_weight = atoi(optarg), opt0.min_chain_weight = 1;
      else if (c == 'y') opt->max_mem_intv = atol(optarg), opt0.max_mem_intv = 1;
      else if (c == 'C') aux->copy_comment = 1;
      else if (c == 'K') fixed_chunk_size = atoi(optarg);
      else if (c == 'X') opt->mask_level = atof(optarg);
#ifdef USE_HTSLIB
      else if (c == 'o') ; //opt->bam_output = atoi(optarg), opt0.bam_output = 1;
#endif
      else if (c == 'h') {
        opt0.max_XA_hits = opt0.max_XA_hits_alt = 1;
        opt->max_XA_hits = opt->max_XA_hits_alt = strtol(optarg, &p, 10);
        if (*p != 0 && ispunct(*p) && isdigit(p[1]))
          opt->max_XA_hits_alt = strtol(p+1, &p, 10);
      }
      else if (c == 'Q') {
        opt0.mapQ_coef_len = 1;
        opt->mapQ_coef_len = atoi(optarg);
        opt->mapQ_coef_fac = opt->mapQ_coef_len > 0? log(opt->mapQ_coef_len) : 0;
      } else if (c == 'O') {
        opt0.o_del = opt0.o_ins = 1;
        opt->o_del = opt->o_ins = strtol(optarg, &p, 10);
        if (*p != 0 && ispunct(*p) && isdigit(p[1]))
          opt->o_ins = strtol(p+1, &p, 10);
      } else if (c == 'E') {
        opt0.e_del = opt0.e_ins = 1;
        opt->e_del = opt->e_ins = strtol(optarg, &p, 10);
        if (*p != 0 && ispunct(*p) && isdigit(p[1]))
          opt->e_ins = strtol(p+1, &p, 10);
      } else if (c == 'L') {
        opt0.pen_clip5 = opt0.pen_clip3 = 1;
        opt->pen_clip5 = opt->pen_clip3 = strtol(optarg, &p, 10);
        if (*p != 0 && ispunct(*p) && isdigit(p[1]))
          opt->pen_clip3 = strtol(p+1, &p, 10);
      } else if (c == 'R') {
        if ((rg_line = bwa_set_rg(optarg)) == 0) return 1; // FIXME: memory leak
      } else if (c == 'H') {
        if (optarg[0] != '@') {
          FILE *fp;
          if ((fp = fopen(optarg, "r")) != 0) {
            char *buf;
            buf = (char*)calloc(1, 0x10000);
            while (fgets(buf, 0xffff, fp)) {
              i = strlen(buf);
              assert(buf[i-1] == '\n'); // a long line
              buf[i-1] = 0;
              hdr_line = bwa_insert_header(buf, hdr_line);
            }
            free(buf);
            fclose(fp);
          }
        } else hdr_line = bwa_insert_header(optarg, hdr_line);
      } else if (c == 'I') { // specify the insert size distribution
        aux->pes0 = pes;
        pes[1].failed = 0;
        pes[1].avg = strtod(optarg, &p);
        pes[1].std = pes[1].avg * .1;
        if (*p != 0 && ispunct(*p) && isdigit(p[1]))
          pes[1].std = strtod(p+1, &p);
        pes[1].high = (int)(pes[1].avg + 4. * pes[1].std + .499);
        pes[1].low  = (int)(pes[1].avg - 4. * pes[1].std + .499);
        if (pes[1].low < 1) pes[1].low = 1;
        if (*p != 0 && ispunct(*p) && isdigit(p[1]))
          pes[1].high = (int)(strtod(p+1, &p) + .499);
        if (*p != 0 && ispunct(*p) && isdigit(p[1]))
          pes[1].low  = (int)(strtod(p+1, &p) + .499);
        if (bwa_verbose >= 3)
          fprintf(stderr, "[M::%s] mean insert size: %.3f, stddev: %.3f, max: %d, min: %d\n",
              __func__, pes[1].avg, pes[1].std, pes[1].high, pes[1].low);
      }
      else return 1;
    }

    if (rg_line) {
      hdr_line = bwa_insert_header(rg_line, hdr_line);
      free(rg_line);
    }

    if (opt->n_threads < 1) opt->n_threads = 1;
    if (optind + 1 >= argc || optind + 3 < argc) {
      fprintf(stderr, "\n");
      fprintf(stderr, "Usage: bwa mem [options] <idxbase> <in1.fq> [in2.fq]\n\n");
      fprintf(stderr, "Algorithm options:\n\n");
      fprintf(stderr, "       -t INT        number of threads [%d]\n", opt->n_threads);
      fprintf(stderr, "       -k INT        minimum seed length [%d]\n", opt->min_seed_len);
      fprintf(stderr, "       -w INT        band width for banded alignment [%d]\n", opt->w);
      fprintf(stderr, "       -d INT        off-diagonal X-dropoff [%d]\n", opt->zdrop);
      fprintf(stderr, "       -r FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT [%g]\n", opt->split_factor);
      fprintf(stderr, "       -y INT        seed occurrence for the 3rd round seeding [%ld]\n", (long)opt->max_mem_intv);
      //		fprintf(stderr, "       -s INT        look for internal seeds inside a seed with less than INT occ [%d]\n", opt->split_width);
      fprintf(stderr, "       -c INT        skip seeds with more than INT occurrences [%d]\n", opt->max_occ);
      fprintf(stderr, "       -D FLOAT      drop chains shorter than FLOAT fraction of the longest overlapping chain [%.2f]\n", opt->drop_ratio);
      fprintf(stderr, "       -W INT        discard a chain if seeded bases shorter than INT [0]\n");
      fprintf(stderr, "       -m INT        perform at most INT rounds of mate rescues for each read [%d]\n", opt->max_matesw);
      fprintf(stderr, "       -S            skip mate rescue\n");
      fprintf(stderr, "       -P            skip pairing; mate rescue performed unless -S also in use\n");
      fprintf(stderr, "\nScoring options:\n\n");
      fprintf(stderr, "       -A INT        score for a sequence match, which scales options -TdBOELU unless overridden [%d]\n", opt->a);
      fprintf(stderr, "       -B INT        penalty for a mismatch [%d]\n", opt->b);
      fprintf(stderr, "       -O INT[,INT]  gap open penalties for deletions and insertions [%d,%d]\n", opt->o_del, opt->o_ins);
      fprintf(stderr, "       -E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [%d,%d]\n", opt->e_del, opt->e_ins);
      fprintf(stderr, "       -L INT[,INT]  penalty for 5'- and 3'-end clipping [%d,%d]\n", opt->pen_clip5, opt->pen_clip3);
      fprintf(stderr, "       -U INT        penalty for an unpaired read pair [%d]\n\n", opt->pen_unpaired);
      fprintf(stderr, "       -x STR        read type. Setting -x changes multiple parameters unless overriden [null]\n");
      fprintf(stderr, "                     pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref)\n");
      fprintf(stderr, "                     ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref)\n");
      fprintf(stderr, "                     intractg: -B9 -O16 -L5  (intra-species contigs to ref)\n");
      fprintf(stderr, "\nInput/output options:\n\n");
      fprintf(stderr, "       -p            smart pairing (ignoring in2.fq)\n");
      fprintf(stderr, "       -R STR        read group header line such as '@RG\\tID:foo\\tSM:bar' [null]\n");
      fprintf(stderr, "       -H STR/FILE   insert STR to header if it starts with @; or insert lines in FILE [null]\n");
      fprintf(stderr, "       -j            treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)\n");
      fprintf(stderr, "\n");
      fprintf(stderr, "       -v INT        verbose level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", bwa_verbose);
      fprintf(stderr, "       -T INT        minimum score to output [%d]\n", opt->T);
      fprintf(stderr, "       -h INT[,INT]  if there are <INT hits with score >80%% of the max score, output all in XA [%d,%d]\n", opt->max_XA_hits, opt->max_XA_hits_alt);
      fprintf(stderr, "       -a            output all alignments for SE or unpaired PE\n");
      fprintf(stderr, "       -C            append FASTA/FASTQ comment to SAM output\n");
      fprintf(stderr, "       -V            output the reference FASTA header in the XR tag\n");
      fprintf(stderr, "       -Y            use soft clipping for supplementary alignments\n");
      fprintf(stderr, "       -M            mark shorter split hits as secondary\n\n");
      fprintf(stderr, "       -I FLOAT[,FLOAT[,INT[,INT]]]\n");
      fprintf(stderr, "                     specify the mean, standard deviation (10%% of the mean if absent), max\n");
      fprintf(stderr, "                     (4 sigma from the mean if absent) and min of the insert size distribution.\n");
      fprintf(stderr, "                     FR orientation only. [inferred]\n");
      fprintf(stderr, "\n");
      fprintf(stderr, "Note: Please read the man page for detailed description of the command line and options.\n");
      fprintf(stderr, "\n");
      free(opt);
      return 1;
    }

    if (mode) {
      if (strcmp(mode, "intractg") == 0) {
        if (!opt0.o_del) opt->o_del = 16;
        if (!opt0.o_ins) opt->o_ins = 16;
        if (!opt0.b) opt->b = 9;
        if (!opt0.pen_clip5) opt->pen_clip5 = 5;
        if (!opt0.pen_clip3) opt->pen_clip3 = 5;
      } else if (strcmp(mode, "pacbio") == 0 || strcmp(mode, "pbref") == 0 || strcmp(mode, "ont2d") == 0) {
        if (!opt0.o_del) opt->o_del = 1;
        if (!opt0.e_del) opt->e_del = 1;
        if (!opt0.o_ins) opt->o_ins = 1;
        if (!opt0.e_ins) opt->e_ins = 1;
        if (!opt0.b) opt->b = 1;
        if (opt0.split_factor == 0.) opt->split_factor = 10.;
        if (strcmp(mode, "ont2d") == 0) {
          if (!opt0.min_chain_weight) opt->min_chain_weight = 20;
          if (!opt0.min_seed_len) opt->min_seed_len = 14;
          if (!opt0.pen_clip5) opt->pen_clip5 = 0;
          if (!opt0.pen_clip3) opt->pen_clip3 = 0;
        } else {
          if (!opt0.min_chain_weight) opt->min_chain_weight = 40;
          if (!opt0.min_seed_len) opt->min_seed_len = 17;
          if (!opt0.pen_clip5) opt->pen_clip5 = 0;
          if (!opt0.pen_clip3) opt->pen_clip3 = 0;
        }
      } else {
        fprintf(stderr, "[E::%s] unknown read type '%s'\n", __func__, mode);
        return 1; // FIXME memory leak
      }
    } else update_a(opt, &opt0);
    bwa_fill_scmat(opt->a, opt->b, opt->mat);

    aux->idx = 0;
    if (aux->idx == 0) {
      if ((aux->idx = bwa_idx_load(argv[optind], BWA_IDX_ALL)) == 0) return 1; // FIXME: memory leak
    } else if (bwa_verbose >= 3)
      fprintf(stderr, "[M::%s] load the bwa index from shared memory\n", __func__);
    if (ignore_alt)
      for (i = 0; i < aux->idx->bns->n_seqs; ++i)
        aux->idx->bns->anns[i].is_alt = 0;

    if (is_master) {
      ko_read1 = kopen(argv[optind + 1], &fd);
      if (ko_read1 == 0) {
        if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind + 1]);
        return 1;
      }
      fp_idx = gzdopen(fd, "r");
      aux->ks = kseq_init(fp_idx);

#if 0 // Deprecated: do not turn on/off FPGA based on read size
      // Decide FPGA usage based on read size
      int read_length = kseq_read(aux->ks);
      if (read_length >= 250 || read_length < 100){
        if (FLAGS_use_fpga) {
          DLOG(WARNING) << "Disabling FPGA for read length = " << read_length;
        }
        FLAGS_use_fpga = false;
      }

      kseq_destroy(aux->ks);
      err_gzclose(fp_idx);
      kclose(ko_read1);
      ko_read1 = kopen(argv[optind + 1], &fd);
      fp_idx = gzdopen(fd, "r");
      aux->ks = kseq_init(fp_idx);
#endif
      if (optind + 2 < argc) {
        if (opt->flag&MEM_F_PE) {
          if (bwa_verbose >= 2)
            fprintf(stderr, "[W::%s] when '-p' is in use, the second query file is ignored.\n", __func__);
        } else {
          ko_read2 = kopen(argv[optind + 2], &fd2);
          if (ko_read2 == 0) {
            if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind + 2]);
            return 1;
          }
          fp2_read2 = gzdopen(fd2, "r");
          aux->ks2 = kseq_init(fp2_read2);
          opt->flag |= MEM_F_PE;
        }
      }
      aux->actual_chunk_size = fixed_chunk_size > 0? fixed_chunk_size : opt->chunk_size * opt->n_threads;
    }
#ifdef USE_HTSLIB
    bam_hdr_t *h = NULL; // TODO
    kstring_t str;
    str.l = str.m = 0; str.s = 0;
    bwa_format_sam_hdr(aux->idx->bns, hdr_line, &str);
    h = sam_hdr_parse(str.l, str.s);
    h->l_text = str.l; h->text = str.s;
    //sam_hdr_write(out, h);
    aux->h = h;
#else
    bwa_print_sam_hdr(aux->idx->bns, hdr_line);
#endif

    return 0;
  }


