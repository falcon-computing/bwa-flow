//#include <autopilot_tech.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ap_int.h>
#include <string.h>
#include <ap_utils.h>
#include <hls_stream.h>
using namespace hls;

typedef ap_uint<2> uint2_t;
typedef ap_uint<4> uint4_t;
typedef ap_uint<11> uint11_t;
typedef ap_uint<32> uuint32_t;
//typedef ap_uint<64> uint64_t;
typedef unsigned char uint8_t;
typedef short int16_t;

//#define PE_NUMS 20
//#define PEARRAY_NUM 2 
//#define PART_TASK_NUMS 40
#define READ_PE_NUM 20
#define o_del 6
#define e_del 1
#define o_ins 6
#define e_ins 1
#define pen_clip 5
#define w_in 100



void sw_extend(int qs_baddr, int ts_baddr, uint4_t *read_seq, uint4_t *chain_rseqs, 
    uint8_t qlen, uint11_t tlen, int16_t h0, int16_t *regScore, int16_t qBeg, int16_t max_ins,
    int16_t max_del, int16_t *w_ret, int16_t *qle_ret, int16_t *tle_ret, int16_t *gtle_ret,
    int16_t *gscore_ret, int16_t *maxoff_ret, uint4_t dir);

void data_parse(int *totalinp, int readNum, stream<int> readTask[READ_PE_NUM],
    stream<uint2_t> readTask_ctrl[READ_PE_NUM]);
void seed_proc( int64_t rmax_0,
                int64_t rmax_1,
                int read_seq_length,
                int seed_index,
                int64_t seed_rbeg,
                int seed_qbeg,
                int seed_len,
                uint4_t* read_seq,
                uint4_t* chain_rseqs,
                stream<int>& read_match
    );
char get_readMatch(stream<int>& read_match, int *results);
void results_assemble(stream<int> read_match[READ_PE_NUM], stream<int>& results);
char fill_resulBuf(stream<int>& results, int *resultsBuf, short *locAddr);
void upload_results(stream<int>& results, int *output_a);

void sw_extend(int qs_baddr, int ts_baddr, uint4_t *read_seq, uint4_t *chain_rseqs, 
    uint8_t qlen, uint11_t tlen, int16_t h0, int16_t *regScore, int16_t qBeg, int16_t max_ins,
    int16_t max_del, int16_t *w_ret, int16_t *qle_ret, int16_t *tle_ret, int16_t *gtle_ret,
    int16_t *gscore_ret, int16_t *maxoff_ret, uint4_t dir)
{
	ap_int<12> i;
	int j;
	ap_uint<2> k;
	//	ap_uint<8> l=0;
	ap_int<12> max_i, max_ie, max_off;
	ap_int<12> gscore;
	int qs_baddr_t;
  int ts_baddr_t;
	short max_j;
	char oe_del = o_del + e_del;
	char oe_ins = o_ins + e_ins;
	uint8_t beg, end;
	uint8_t backw_tmp;
	char forw_update;
	uint8_t forw_tmp;
	ap_int<10> abs_mj_m_i;
	char tmp_ehh_m_eins;
  int16_t tmp_eme;
	int16_t h1_init_val;
	//	char h1_init_tmp=0;
	int16_t max;
	int16_t h, e,M;
	int16_t e_tmp;
	int16_t h_tmp;
  int16_t h1_reg;
	int16_t t, f, h1, m;
	short mj;
	ap_uint<3> q_i, q_j;
	ap_int<10> prev;
	char isBreak;
	uint8_t aw1;
	uint8_t aw_tmp;

	char h0_arr[2];
#pragma HLS ARRAY_PARTITION variable=h0_arr complete dim=0

	const char my_mat[5][5]={{1, -4, -4, -4, -1}, {-4, 1, -4, -4, -1}, {-4, -4, 1, -4, -1}, {-4, -4, -4, 1, -1}, {-1, -1, -1, -1, -1}};
#pragma HLS ARRAY_PARTITION variable=my_mat complete dim=0
	int16_t eh_h [512];
#pragma HLS ARRAY_MAP variable=eh_h instance=eh_arr vertical
#pragma HLS RESOURCE variable=eh_h core=RAM_2P_BRAM
	int16_t eh_e [512];
#pragma HLS ARRAY_MAP variable=eh_e instance=eh_arr vertical
#pragma HLS RESOURCE variable=eh_e core=RAM_2P_BRAM

	max = h0;
	max_i = max_j = -1;
	max_ie = -1;
	gscore = -1;
	max_off = 0;

	k = 0;
	isBreak = 0;
for(j=0;j<=qlen;j++)
 { 
#pragma HLS PIPELINE II=1
  eh_e[j]= 0;
  eh_h[j]= 0;
}

ext_while_loop : while ((k < 2) && (!isBreak))
				 {
#pragma HLS LOOP_TRIPCOUNT min=2 max=2
					 prev = *regScore;
					 aw_tmp = w_in << k;
					 aw1 = aw_tmp < max_ins ? aw_tmp : max_ins;
					 aw1 = aw1 < max_del ? aw1 : max_del;
					 beg = 0;
					 end = qlen;
					 tmp_eme = h0 - oe_ins;
					 tmp_eme = (tmp_eme > 0) ? tmp_eme : 0;
					 h1_init_val = h0 - o_del;
target_loop : for (i = 0; i < tlen; i++) {
#pragma HLS LOOP_TRIPCOUNT min=10 max=10
					 f = 0; m = 0; mj = -1;
					 // ts_baddr_t = ts_baddr + i;
           if (dir == 0)
             ts_baddr_t = ts_baddr - i;
           else
             ts_baddr_t = ts_baddr + i;
					 q_i = chain_rseqs[ts_baddr_t];

					 if (beg < i - aw1) beg = i - aw1;
					 //#pragma HLS resource variable=beg core=AddSub_DSP
					 if (end > i + aw1 + 1) end = i + aw1 + 1;
					 //#pragma HLS resource variable=end core=AddSub_DSP
					 if (end > qlen) end = qlen;
					 if(beg ==0){
             h1_init_val -= e_del;
					   h1 = h1_init_val;
					   if (h1 < 0) h1 = 0;
           }
           else h1 = 0;
					 backw_tmp = 0; 
					 forw_tmp = 0; 
					 forw_update = 0;
query_loop : for (j = beg; j < end; ++j) {
#pragma HLS LOOP_TRIPCOUNT min=10 max=10
#pragma AP pipeline II=6
//#pragma AP dependence variable=eh_e array inter false
//#pragma AP dependence variable=eh_h array inter false
           if (dir == 0)
             qs_baddr_t = qs_baddr - j;
           else
             qs_baddr_t = qs_baddr + j;
					 q_j = read_seq[qs_baddr_t];
					 h_tmp = eh_h[j];// get H(i-1,j-1) and E(i-1,j)
					 e_tmp = eh_e[j];
					 if (i == 0) {
						 e = 0;
						 if (j == 0) {
							 h = h0;
               M = h0;
						 }
						 else if (j == 1) {
							 h = tmp_eme;
               M = tmp_eme;
						 }
						 else {
							 tmp_eme -= e_ins;
							 h = (tmp_eme > 0) ? tmp_eme : 0;
							 M = (tmp_eme > 0) ? tmp_eme : 0;
						 }
					 }
					 else {
						 e = e_tmp;
						 h = h_tmp;
             M = h_tmp;
					 }
					 h1_reg = h1;
					 //h += my_mat[q_i][q_j];
                                         M = M?M+my_mat[q_i][q_j]:0;
					 h = M > e? M : e;
					 h = h > f? h : f;
					 h1 = h;             // save H(i,j) to h1 for the next column
					 
					 t = M - oe_del;
					 t = (t > 0) ? t : 0;
					 e -= e_del;
					 e = (e > t) ? e : t;   // computed E(i+1,j)
					 t = M - oe_ins;
					 t = t > 0? t : 0;
					 f -= e_ins;
					 f = (f > t) ? f : t;   // computed F(i,j+1)
					 eh_e[j] = e; // save E(i+1,j) for the next row
					 eh_h[j] = h1_reg;          // set H(i,j-1) for the next row
                                        /* if (h1_reg == 0 && e ==0) {
						 backw_tmp = 0;
					 }
					 else {
						 backw_tmp++;
					 }*/
					 if (m <= h)
					 {
						 mj = j;
						 m = h;
					 }
           if (forw_update == 0) { //((h1_reg == 0) &&
             if (h1_reg == 0 && e ==0) {
               forw_tmp++;
             }
             else {
               forw_update = 1;
             }
           }
           if(h1_reg ==0 && e ==0){
             backw_tmp++;
           } 
           else backw_tmp = 0;
                                      
                                         
			 }
			 eh_h[end] = h1;
			 eh_e[end] = 0; 
       if( h1 == 0) 
         backw_tmp++; 
       else 
         backw_tmp = 0; 
			/* if ((forw_update == 0) && (h1 != 0)) {
				 if ((j >= mj+2) || (forw_tmp != 0)) {
					 forw_tmp++;
				 }
			 }*/
			 if (j == qlen) {
				 if (gscore <= h1) {
					 max_ie = i;
					 gscore = h1;
				 }
			 }
			 if (m == 0) break;
			 if (m > max) {
				 max = m;
				 max_i = i;
				 max_j = mj;
				 abs_mj_m_i = abs(mj - i);
				 //#pragma HLS resource variable=abs_mj_m_i core=AddSub_DSP
				 if (max_off < abs_mj_m_i) max_off = abs_mj_m_i;
			 }
			 beg = beg + forw_tmp;
			 end = end - backw_tmp + 2 <qlen ? end-backw_tmp +2:qlen;
		}
		*qle_ret = max_j + 1;
		*tle_ret = max_i + 1;
		*gtle_ret = max_ie + 1;
		*gscore_ret = gscore;
		*maxoff_ret = max_off;
		*regScore = max;
		if (max == prev || ( max_off < (aw_tmp >> 1) + (aw_tmp >> 2))) isBreak = 1;
		k++;
	}
	*w_ret = aw_tmp;
}

// distribute the data by read
void data_parse(int *totalinp, int readNum, stream<int> readTask[READ_PE_NUM],
    stream<uint2_t> readTask_ctrl[READ_PE_NUM])
{
  int readStartIndex = 0;
  int readEndIndex ;
  int readDataLength = 0 ;
  int i = 0;
  int j = 0;
  int k = 0;
  data_parse_loop :
  for (i=0; i<readNum; i++) {
#pragma HLS LOOP_TRIPCOUNT min=100 max=200
    readEndIndex = totalinp[readStartIndex];
    while (1){
#pragma HLS LOOP_TRIPCOUNT min=1 max=2
      if (!readTask_ctrl[j].full()){
        for (k= readStartIndex + 1 ; k < readEndIndex; k++){
#pragma HLS LOOP_TRIPCOUNT min=60 max=100
#pragma HLS pipeline II=1
          readTask[j].write(totalinp[k]);
        }
        readStartIndex = readEndIndex;
        readTask_ctrl[j].write(1);
        j++;
        if (j >= READ_PE_NUM){
          j = 0;
        }
        break;
      }
      else{
        j++;
        if (j >= READ_PE_NUM){
          j = 0;
        }
      }
    }
  }
  // write the end signal
  for (i=0; i<READ_PE_NUM; i++){
#pragma HLS UNROLL
    readTask_ctrl[i].write(3);
  }
}

void seed_proc( int64_t rmax_0,
                int64_t rmax_1,
                int read_seq_length,
                int seed_index,
                int64_t seed_rbeg,
                int seed_qbeg,
                int seed_len,
                uint4_t* read_seq,
                uint4_t* chain_rseqs,
                stream<int>& read_match
    )
{
// for now I dont use flexible opts 
#pragma HLS INLINE  
  uint8_t qlen[2];
  uint11_t tlen[2];
  int16_t max_ins[2];
  int16_t max_del[2];
  int16_t h0 = 0;
  int16_t regScore = 0;

  int16_t qle =0;
  int16_t tle =0;
  int16_t gtle =0;
  int16_t gscore =0;
  int16_t maxoff =0;
  int16_t qBeg =0;
	int16_t rBeg=0;
	int16_t qEnd=0;
	int16_t rEnd=0;
	int16_t score=0;
	int16_t trueScore=0;
	int16_t width=0;
  int16_t aw[2];
  int16_t sc0 =0;
  int16_t h0_arr[2];
  int16_t qs_baddr;
  int16_t ts_baddr;

  uint4_t i = 0;

  qlen[0] = seed_qbeg;
  tlen[0] = seed_rbeg - rmax_0;
  qlen[1] = read_seq_length - seed_len - seed_qbeg;
  tlen[1] = rmax_1 - seed_len - seed_rbeg;
  h0 = seed_len;
  regScore = seed_len;
  max_ins[0] = qlen[0]; // for static opt
  max_del[0] = qlen[0];
  max_ins[1] = qlen[1];
  max_del[1] = qlen[1];
  aw[0] = w_in; 
  aw[1] = w_in;
	qBeg = 0;
	qEnd = qlen[1];
	rBeg = 0;
	rEnd = 0;
	trueScore = regScore;
	qle = -1;
	tle = -1;
	gtle = -1;
	gscore = -1;
	maxoff = -1;
  
  qs_baddr = qlen[0] -1; 
  ts_baddr = tlen[0] -1;

  for (i=0; i<2; i++){
    sc0 = regScore;
    h0_arr[0] = h0;
    h0_arr[1] = sc0;
    sw_extend(
        qs_baddr,
        ts_baddr,
        read_seq,
        chain_rseqs,
        qlen[i],
        tlen[i],
        h0_arr[i],
        &regScore,
        seed_qbeg,
        max_ins[i],
        max_del[i],
        &aw[i],
        &qle,
        &tle,
        &gtle,
        &gscore,
        &maxoff,
        i
      );
    score = regScore;
		if (gscore <= 0 || gscore <= (regScore - pen_clip)) {
		  if (i == 0) {
			  qBeg = seed_qbeg - qle;
			  rBeg = -tle;
			  trueScore = regScore;
		  }
		  else {
			  qEnd = qle;
			  rEnd = tle;
			  trueScore += regScore - sc0;
		  }
		}
		else {
		  if (i == 0) {
			  qBeg = 0;
			  rBeg = -gtle;
			  trueScore = gscore;
		  }
		  else {
			  qEnd = qlen[1];
			  rEnd = gtle;
			  trueScore += gscore - sc0;
		  }
		}
    qs_baddr = qlen[0] + seed_len;
    ts_baddr = tlen[0] + seed_len;
  }
  if (aw[0] > aw[1]) width = aw[0];
  else width = aw[1];
  // get the results
  read_match.write(seed_index);
  read_match.write((qBeg & 0xFFFF) | ((qEnd<<16) & 0xFFFF0000));
  read_match.write((rBeg & 0xFFFF) | ((rEnd<<16) & 0xFFFF0000));
  read_match.write((score & 0xFFFF) | ((trueScore<<16) & 0xFFFF0000));
  read_match.write(width & 0xFFFF);
}

void read_proc(stream<int>& readTask, stream<uint2_t>& readTask_ctrl, stream<int>& read_match)
{
  uint4_t read_seq[1024];
  uint4_t chain_rseqs[1024];
  
  int read_seq_length;
  int read_seq_length_div8;
  int rseq_length;
  int rseq_length_div8;
  int i =0;
  int j =0;
  char k =0;
  int tmp_int_data;
  int tmp_int_data_r;
  int chain_num;
  int seed_num;
  int64_t rmax_0;
  int64_t rmax_1;
  int seed_index;
  int seed_qbeg;
  int64_t seed_rbeg;
  int seed_len;
  uint2_t endFlag0=0; 
  int debug_read_num = 0;

  while(1){
    while(1){
      if (!readTask_ctrl.empty())
        break;
    }
    readTask_ctrl.read(endFlag0);
    if (endFlag0 == 3) {
      read_match.write(0xFFFFFFFF);
      break;
    }
    else{
      debug_read_num += 1;
      readTask.read(read_seq_length);
      if ((read_seq_length & 0x007) != 0) {
        read_seq_length_div8 = (read_seq_length >> 3) + 1;
      }
      else {
        read_seq_length_div8 = read_seq_length >> 3;
      }
      // get the read sequence
      for (i=0; i<read_seq_length_div8; i++){
        readTask.read(tmp_int_data);
        for (j=0; j<8; j++) {
          read_seq[i*8 + j] = (tmp_int_data & 0xF0000000) >> 28;
          tmp_int_data <<=4;
        }
      }
      readTask.read(chain_num);
      // get the chain information
      for (i=0; i< chain_num; i++){
        readTask.read(tmp_int_data);
        readTask.read(tmp_int_data_r);
        // lower bits come out first
        rmax_0 = (uuint32_t)tmp_int_data_r;
        rmax_0 <<= 32;
        rmax_0 |= (uuint32_t)tmp_int_data;
        readTask.read(tmp_int_data);
        readTask.read(tmp_int_data_r);
        rmax_1 = (uuint32_t)tmp_int_data_r;
        rmax_1 <<= 32;
        rmax_1 |= (uuint32_t)tmp_int_data;
        rseq_length = rmax_1 - rmax_0 ;
        if ((rseq_length & 0x007) != 0) {
          rseq_length_div8 = (rseq_length >> 3) + 1;
        }
        else {
          rseq_length_div8 = rseq_length >> 3;
        }
        // get the chain rseq
        for (j=0; j<rseq_length_div8; j++){
          readTask.read(tmp_int_data);
          for (k=0; k<8; k++){
            chain_rseqs[j*8 +k] = (tmp_int_data & 0xF0000000) >> 28;
            tmp_int_data <<=4;
          }
        }
        // proc the seeds
        readTask.read(seed_num);
        for (j=0; j<seed_num; j++){
          readTask.read(seed_index);
          readTask.read(tmp_int_data);
          readTask.read(tmp_int_data_r);
          seed_rbeg = (uuint32_t)tmp_int_data_r;
          seed_rbeg <<= 32;
          seed_rbeg |= (uuint32_t)tmp_int_data;
          readTask.read(seed_qbeg);
          readTask.read(seed_len);
          seed_proc(
                rmax_0,
                rmax_1,
                read_seq_length,
                seed_index,
                seed_rbeg,
                seed_qbeg,
                seed_len,
                read_seq,
                chain_rseqs,
                read_match
              );
        }
      }
    }
  }
}

char get_peMatch(stream<int>& pe_matchs, int *oneResultBuf)
{
#pragma HLS INLINE
	char i=0;

	if (!pe_matchs.empty()) {
		pe_matchs.read(oneResultBuf[0]);
		if (oneResultBuf[0] == 0xFFFFFFFF) {
			return 2;
		}
		else {
			for (i = 0; i < 4; i++) {
#pragma HLS PIPELINE II=1
				pe_matchs.read(oneResultBuf[1 + i]);
			}
			return 1;
		}
	}
	else {
		return 0;
	}
}

char get_readMatch(stream<int>& read_match, int *results)
{
#pragma HLS INLINE
  char i = 0;
  int seed_index = 0;
  
  if (!read_match.empty()) {
    read_match.read(seed_index);
    if (seed_index == 0xFFFFFFFF) {
      return 2;
    }
    else {
      results[0] = seed_index;
      for (i=1; i < 5; i++) {
#pragma HLS PIPELINE II=1
        read_match.read(results[i]);
      }
      return 1;
    }
  }
  else {
    return 0;   
  }
}

void results_assemble(stream<int> read_match[READ_PE_NUM], stream<int>& results)
{
  char j = 0;
  char k = 0;
  char getReadMatchFlag = 0;
  char readOver_cnt = 0;
  int oneMatchBuf[5];

  while(1) {
    getReadMatchFlag = get_readMatch(read_match[j], oneMatchBuf);
    if (getReadMatchFlag == 0) {
      j++;
      if (j >= READ_PE_NUM) j = 0;
    }
    else if (getReadMatchFlag == 1){
      for (k = 0; k < 5; k++) {
#pragma HLS PIPELINE II=1
        results.write(oneMatchBuf[k]);
      }
    }
    else if (getReadMatchFlag == 2){
      readOver_cnt++;
    } 
    if (readOver_cnt >= READ_PE_NUM) {
      results.write(0xFFFFFFFF);
      break;
    }
  }
}

char fill_resulBuf(stream<int>& results, int *resultsBuf, short *locAddr)
{
	char i=0;
	int seed_index=0;
	int locAddr_t;

	locAddr_t = *locAddr;
	while(1) {
		if (!results.empty()) {
			results.read(seed_index);
			if (seed_index == 0xFFFFFFFF) {
				*locAddr = locAddr_t;
//printf("here!!!\n");
				if (locAddr_t != 0) {
					return 1;
				}
				else {
				    return 2;
				}
			}
			else {
				resultsBuf[locAddr_t] = seed_index;
				locAddr_t = locAddr_t + 1;
				for (i=0; i<4; i++) {
#pragma HLS PIPELINE II=1
					results.read(resultsBuf[locAddr_t]);
					locAddr_t = locAddr_t + 1;
				}
			    if (locAddr_t >= 1000) {
					*locAddr = locAddr_t;
					return 0;
			    } 
			}
		}
	}
}

void upload_results(stream<int>& results, int *output_a)
{
	char flag_over;
	short locAddr;
	int outAddr;
	int resultsBuf[1024];
  int i;
  int tmp_data;


	outAddr = 0;
	while(1) {
		locAddr = 0;
		flag_over = fill_resulBuf(results, resultsBuf, &locAddr);
		if (flag_over == 2) {
			break;
		}
		else if (flag_over == 1) {
			//memcpy(&output_a[outAddr], resultsBuf, locAddr*4);
      for ( i=0; i < locAddr; i++){
        output_a[outAddr + i] = resultsBuf[i];
      }
			break;
		}
		else if (flag_over == 0) {
			//memcpy(&output_a[outAddr], resultsBuf, locAddr*4);
      for ( i=0; i < locAddr; i++){
        output_a[outAddr + i] = resultsBuf[i];
      }
			outAddr = outAddr + locAddr;
		}
	}
}

void sw_core(int *a, int *output_a, int __inc)
{
#pragma HLS inline off
  int local_inc = 0;
  const int EXTRA = 0;
  local_inc = __inc;
  char i = 0;

#pragma HLS DATAFLOW
  stream<int> readTask[READ_PE_NUM];
  #pragma HLS STREAM variable=readTask depth=10240+EXTRA
  stream<uint2_t> readTask_ctrl[READ_PE_NUM];
  #pragma HLS STREAM variable=readTask_ctrl depth=1+EXTRA
  stream<int> read_match[READ_PE_NUM];
  #pragma HLS STREAM variable=read_match depth=512+EXTRA
  stream<int> result_out;
  #pragma HLS STREAM variable=result_out depth=2048+EXTRA

  data_parse(a, local_inc, readTask, readTask_ctrl);
  
  for (i=0; i<READ_PE_NUM; i++) {
    #pragma HLS UNROLL
    read_proc(readTask[i], readTask_ctrl[i], read_match[i] );
  }

  results_assemble(read_match, result_out);
  upload_results(result_out, output_a);
  return;
}


#ifndef HLS_
extern "C" {
#endif
void sw_top(int *a, int *output_a, int __inc)
{
#ifndef HLS_
#pragma HLS INTERFACE m_axi port=a offset=slave bundle=gmem
#pragma HLS INTERFACE m_axi port=output_a offset=slave bundle=gmem
#else
#pragma HLS INTERFACE m_axi port=a offset=slave bundle=gmem depth=20851
#pragma HLS INTERFACE m_axi port=output_a offset=slave bundle=gmem depth=1024*5
#endif
#pragma HLS INTERFACE s_axilite port=a bundle=control
#pragma HLS INTERFACE s_axilite port=output_a bundle=control
#pragma HLS INTERFACE s_axilite port=__inc bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control

    sw_core(a, output_a, __inc);
   return;
}
#ifndef HLS_
}
#endif
