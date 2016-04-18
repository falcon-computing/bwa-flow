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
typedef unsigned char uint8_t;
typedef short int16_t;

//#define PE_NUMS 1
#define PE_NUMS 20
#define PEARRAY_NUM 1
#define PART_TASK_NUMS 40
//#define LOG2_PART_TASK_NUMS 5
#define BLKBUF_DEPTH 3072

//#define TOTAL_TASK_NUMS 4096 
//#define RESULT_SIZE TOTAL_TASK_NUMS*4

void sw_extend(uint11_t qs_baddr, uint4_t *qs, uint11_t ts_baddr, uint8_t qlen, uint11_t tlen, char o_ins,
			   char e_ins, char o_del, char e_del, char penClip, char w_in, int16_t h0, int16_t *regScore, int16_t qBeg, int16_t max_ins, int16_t max_del,
			   int16_t *w_ret, int16_t *qle_ret, int16_t *tle_ret, int16_t *gtle_ret, int16_t *gscore_ret, int16_t *maxoff_ret);
void leftright_ext(stream<int>& pe_seeds, stream<int>& pe_matchs);
void load_data(int *totalinp, int *partInp, int startTaskIndex, int partTaskNums);
void load_task(int *partInp, int startTaskIndex, int partTaskNums, int staOffset, int i, stream<int>& pe_seeds, stream<uint2_t>& pe_seeds_ctrl);
bool feed_peArray(int *blockTaskBuf, stream<int>& peArr_blockTask, stream<uint2_t>& peArr_blockTask_ctrl);
//bool feed_seed(int *oneTaskBuf, stream<int>& pe_seeds, stream<uint2_t>& pe_seeds_ctrl);
void proc_block(stream<int>& peArr_blockTask, stream<uint2_t>& peArr_blockTask_ctrl, stream<int> seeds[PE_NUMS], stream<uint2_t> seeds_ctrl[PE_NUMS]);
void proc_element(stream<int>& pe_seeds, stream<uint2_t>& pe_seeds_ctrl, stream<int>& pe_matchs);
//void write_match(int *oneMatchBuf, stream<int>& pe_matchs);
void task_parse(int *partInp, int startTaskIndex, int partTaskNums, stream<int> seeds[PE_NUMS], stream<uint2_t> seeds_ctrl[PE_NUMS]);
void feed_data(char mode, int *blockTaskMem, stream<int> blockTask[PEARRAY_NUM], stream<uint2_t> blockTask_ctrl[PEARRAY_NUM]);
void data_parse(int *totalinp, int __inc, stream<int> blockTask[PEARRAY_NUM], stream<uint2_t> blockTask_ctrl[PEARRAY_NUM]);
char get_peMatch(stream<int>& pe_matchs, int *oneResultBuf);
void receive_match(stream<int> matchs[PE_NUMS], stream<int>& peArr_blockMatchs);
char get_peArrMatch(stream<int>& peArr_blockMatchs, int *results);
void results_assemble(stream<int> blockMatchs[PEARRAY_NUM], stream<int>& results);
char fill_resulBuf(stream<int>& results, int *resultsBuf, short *locAddr);
void upload_results(stream<int>& results, int *output_a);

void sw_extend(uint11_t qs_baddr, uint4_t *qs, uint11_t ts_baddr, uint8_t qlen, uint11_t tlen, char o_ins,
			   char e_ins, char o_del, char e_del, char penClip, char w_in, int16_t h0, int16_t *regScore, int16_t qBeg, int16_t max_ins, int16_t max_del,
			   int16_t *w_ret, int16_t *qle_ret, int16_t *tle_ret, int16_t *gtle_ret, int16_t *gscore_ret, int16_t *maxoff_ret)
{
	ap_int<12> i;
	uint8_t j;
	ap_uint<2> k;
	//	ap_uint<8> l=0;
	ap_int<12> max_i, max_ie, max_off;
	ap_int<12> gscore;
	uint11_t ts_baddr_t;
	char max_j;
	char oe_del = o_del + e_del;
	char oe_ins = o_ins + e_ins;
	uint8_t beg, end;
	char backw_tmp;
	char backw_reg;
	char forw_update;
	char forw_tmp;
	char forw_reg;
	ap_int<10> abs_mj_m_i;
	char tmp_ehh_m_eins;
        int16_t tmp_eme;
	int16_t h1_init_val;
	//	char h1_init_tmp=0;
	int16_t max;
	int16_t h, e;
	int16_t e_tmp;
	int16_t h_tmp;
        int16_t h1_reg;
	int16_t t, f, h1, m;
	char mj;
	ap_uint<3> q_i, q_j;
	ap_int<10> prev;
	char isBreak;
	char aw1;
	char aw_tmp;
	//	char h_i0=0;

	char h0_arr[2];
#pragma HLS ARRAY_PARTITION variable=h0_arr complete dim=0

	const char my_mat[5][5]={{1, -4, -4, -4, -1}, {-4, 1, -4, -4, -1}, {-4, -4, 1, -4, -1}, {-4, -4, -4, 1, -1}, {-1, -1, -1, -1, -1}};
#pragma HLS ARRAY_PARTITION variable=my_mat complete dim=0
	int16_t eh_h [256];
#pragma HLS ARRAY_MAP variable=eh_h instance=eh_arr vertical
#pragma HLS RESOURCE variable=eh_h core=RAM_2P_BRAM
	int16_t eh_e [256];
#pragma HLS ARRAY_MAP variable=eh_e instance=eh_arr vertical
#pragma HLS RESOURCE variable=eh_e core=RAM_2P_BRAM

	max = h0;
	max_i = max_j = -1;
	max_ie = -1;
	gscore = -1;
	max_off = 0;

	k = 0;
	isBreak = 0;
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
					 ts_baddr_t = ts_baddr + i;
					 //#pragma HLS resource variable=ts_baddr_t core=AddSub_DSP
					 q_i = qs[ts_baddr_t];
					 h1_init_val -= e_del;
					 h1 = h1_init_val;
					 if (h1 < 0) h1 = 0;
					 if (beg < i - aw1) beg = i - aw1;
					 //#pragma HLS resource variable=beg core=AddSub_DSP
					 if (end > i + aw1 + 1) end = i + aw1 + 1;
					 //#pragma HLS resource variable=end core=AddSub_DSP
					 if (end > qlen) end = qlen;
					 backw_tmp = 0; backw_reg = 0;
					 forw_tmp = 0; forw_reg = 0;
					 forw_update = 0;
query_loop : for (j = beg; j < end; ++j) {
#pragma HLS LOOP_TRIPCOUNT min=10 max=10
#pragma AP pipeline II=6
//#pragma AP dependence variable=eh_e array inter false
//#pragma AP dependence variable=eh_h array inter false
					 q_j = qs[qs_baddr + j];
					 h_tmp = eh_h[j];// get H(i-1,j-1) and E(i-1,j)
					 e_tmp = eh_e[j];
					 if (i == 0) {
						 e = 0;
						 if (j == 0) {
							 h = h0;
						 }
						 else if (j == 1) {
							 h = tmp_eme;
						 }
						 else {
							 tmp_eme -= e_ins;
							 h = (tmp_eme > 0) ? tmp_eme : 0;
						 }
					 }
					 else {
						 e = e_tmp;
						 h = h_tmp;
					 }
					 h1_reg = h1;
					 h += my_mat[q_i][q_j];
					 h = h > e? h : e;
					 h = h > f? h : f;
					 h1 = h;             // save H(i,j) to h1 for the next column
					 if (h1_reg == 0) {
						 backw_tmp = 0;
					 }
					 else {
						 backw_tmp++;
					 }
					 if (m <= h)
					 {
						 mj = j;
						 m = h;
						 backw_reg = backw_tmp;
					 }
					 if (j >= mj+2) {
						 if (forw_update == 0) { //((h1_reg == 0) &&
							 if (h1_reg == 0) {
								 forw_update = 1;
							 }
							 else {
								 forw_tmp++;
							 }
						 }
					 }
					 else {
						 forw_tmp = 0;
						 forw_update = 0;
					 }
					 t = h - oe_del;
					 t = (t > 0) ? t : 0;
					 e -= e_del;
					 e = (e > t) ? e : t;   // computed E(i+1,j)
					 t = h - oe_ins;
					 t = t > 0? t : 0;
					 f -= e_ins;
					 f = (f > t) ? f : t;   // computed F(i,j+1)
					 eh_e[j] = e; // save E(i+1,j) for the next row
					 eh_h[j] = h1_reg;          // set H(i,j-1) for the next row
			 }
			 eh_h[end] = h1;
			 eh_e[end] = 0;
			 if ((forw_update == 0) && (h1 != 0)) {
				 if ((j >= mj+2) || (forw_tmp != 0)) {
					 forw_tmp++;
				 }
			 }
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
			 //j = mj - backw_reg;
			 //beg = j + 1;
			 beg = mj - backw_reg + 1;
			 //#pragma HLS resource variable=beg core=AddSub_DSP
			 //j = mj + 2 + forw_tmp;
			 //end = j;
			 end = mj + 2 + forw_tmp;
			 //#pragma HLS resource variable=end core=AddSub_DSP
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


void leftright_ext(stream<int>& pe_seeds, stream<int>& pe_matchs)
{
#pragma HLS INLINE
	char o_del = 0;
	char e_del = 0;
	char o_ins = 0;
	char e_ins = 0;
	char penClip[2];
#pragma HLS ARRAY_PARTITION variable=penClip complete dim=0
	char w_in = 0;
	int result_data = 0;
	int tmp_compar = 0;
	int seed_index = 0;
	uint8_t qlen[2];
#pragma HLS ARRAY_PARTITION variable=qlen complete dim=0
	uint11_t tlen[2];
#pragma HLS ARRAY_PARTITION variable=tlen complete dim=0
	int16_t max_ins[2];
#pragma HLS ARRAY_PARTITION variable=max_ins complete dim=0
	int16_t max_del[2];
#pragma HLS ARRAY_PARTITION variable=max_del complete dim=0
	int16_t h0=0;
	int16_t regScore=0;
	int16_t qBeg_ori=0;

	int16_t qle=0;
	int16_t tle=0;
	int16_t gtle=0;
	int16_t gscore=0;
	int16_t maxoff=0;
	int16_t qBeg=0;
	int16_t rBeg=0;
	int16_t qEnd=0;
	int16_t rEnd=0;
	int16_t score=0;
	int16_t trueScore=0;
	int16_t width=0;

	uint4_t i = 0;
	uint11_t k = 0;
	uint4_t l = 0;
	uint11_t qlen2 = 0;
	uint11_t tlen2 = 0;
	uint11_t qs_baddr = 0;
	uint11_t ts_baddr = 0;
	int16_t aw[2];
#pragma HLS ARRAY_PARTITION variable=aw complete dim=0
	int tmp_parame = 0;
	int qr_offset = 0;
	int16_t sc0 = 0;
	int16_t h0_arr[2];
#pragma HLS ARRAY_PARTITION variable=h0_arr complete dim=0

	uint4_t query_mem[2048];
#pragma HLS RESOURCE variable=query_mem core=RAM_2P_BRAM
	int tmp_qr_data = 0;
//	int param_mem[8];
	int param_data = 0;
	int query_data = 0;
	int ctrl_flag = 0;
	uint11_t qrLen_div8 = 0;

//	pe_seeds.read();
	pe_seeds.read(tmp_compar);
//	tmp_compar = oneSeedBuf[1];
	o_del = tmp_compar & 0xFF; //6
	e_del = (tmp_compar >> 8) & 0xFF; //1
	o_ins = (tmp_compar >> 16) & 0xFF; //6
	e_ins = (tmp_compar >> 24) & 0xFF; //1
	pe_seeds.read(tmp_compar);
//	tmp_compar = oneSeedBuf[2];
	penClip[0] = tmp_compar & 0xFF; //100
	penClip[1] = (tmp_compar >> 8) & 0xFF;
	w_in = (tmp_compar >> 16) & 0xFF;

//load_param :
//	for (k = 0; k < 8; k++) {
//		param_mem[k] = oneSeedBuf[3+k];
//	}

	pe_seeds.read(tmp_parame);
//	 tmp_parame = oneSeedBuf[0+3];
	 qlen[0] = tmp_parame & 0xFFFF; //55
	 tlen[0] = (tmp_parame >> 16) & 0xFFFF; //105
	pe_seeds.read(tmp_parame);
//	 tmp_parame = oneSeedBuf[1+3];
	 qlen[1] = tmp_parame & 0xFFFF; //55
	 tlen[1] = (tmp_parame >> 16) & 0xFFFF; //105
	pe_seeds.read(tmp_parame);
//	 tmp_parame = oneSeedBuf[2+3];
	 qr_offset = tmp_parame;
	pe_seeds.read(tmp_parame);
//	 tmp_parame = oneSeedBuf[3+3];
	 regScore = tmp_parame & 0xFFFF; //100
	 qBeg_ori = (tmp_parame >> 16) & 0xFFFF;
	pe_seeds.read(tmp_parame);
//	 tmp_parame = oneSeedBuf[4+3];
	 h0 = tmp_parame & 0xFFFF;
	pe_seeds.read(tmp_parame);
//	 tmp_parame = oneSeedBuf[5+3];
	 max_ins[0] = tmp_parame & 0xFFFF;
	 max_del[0] = (tmp_parame >> 16) & 0xFFFF;
	pe_seeds.read(tmp_parame);
//	 tmp_parame = oneSeedBuf[6+3];
	 max_ins[1] = tmp_parame & 0xFFFF;
	 max_del[1] = (tmp_parame >> 16) & 0xFFFF;
	pe_seeds.read(seed_index);
//	 seed_index = oneSeedBuf[7+3];

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
	 qlen2 = qlen[0] + qlen[1];
	 tlen2 = tlen[0] + tlen[1];
	 qrLen_div8 = qlen2+tlen2;
//#pragma HLS resource variable=qrLen_div8 core=AddSub_DSP
	 if ((qrLen_div8 & 0x007) != 0) {
		 qrLen_div8 = (qrLen_div8 >> 3) + 1;
	 }
	 else {
		 qrLen_div8 = (qrLen_div8 >> 3);
	 }
load_query :
	 for (k=0; k<qrLen_div8; k++) {
		pe_seeds.read(tmp_qr_data);
//		 tmp_qr_data = oneSeedBuf[11+k];
		 for (l=0; l<8; l++) {
			 query_mem[k*8 + l] = (tmp_qr_data & 0xF0000000) >> 28;
			 tmp_qr_data <<= 4;
		 }
	 }
	 qs_baddr = 0;
	 ts_baddr = qlen2;
left_right_loop :
         score = h0;
	 for (i=0; i<2; i++)
	  {
#pragma HLS LOOP_TRIPCOUNT min=2 max=2
		  sc0 = regScore;
		  h0_arr[0] = h0;
		  h0_arr[1] = sc0;
		  if (qlen[i] > 0) {
			  sw_extend(qs_baddr, query_mem, ts_baddr, qlen[i], tlen[i], o_ins, e_ins, o_del, e_del, penClip[i],
				  w_in, h0_arr[i], &regScore, qBeg_ori, max_ins[i], max_del[i], &aw[i], &qle, &tle, &gtle, &gscore, &maxoff);
			  score = regScore;
			  if (gscore <= 0 || gscore <= (regScore - penClip[i])) {
				  if (i == 0) {
					  qBeg = qBeg_ori - qle;
					//#pragma HLS resource variable=qBeg core=AddSub_DSP
					  rBeg = -tle;
					//#pragma HLS resource variable=rBeg core=AddSub_DSP
					  trueScore = regScore;
				  }
				  else {
					  qEnd = qle;
					  rEnd = tle;
					  trueScore += regScore - sc0;
					//#pragma HLS resource variable=trueScore core=AddSub_DSP
				  }
			  }
			  else {
				  if (i == 0) {
					  qBeg = 0;
					  rBeg = -gtle;
					//#pragma HLS resource variable=rBeg core=AddSub_DSP
					  trueScore = gscore;
				  }
				  else {
					  qEnd = qlen[1];
					  rEnd = gtle;
					  trueScore += gscore - sc0;
					//#pragma HLS resource variable=trueScore core=AddSub_DSP
				  }
			  }
		  }
		  qs_baddr += qlen[i];
		//#pragma HLS resource variable=qs_baddr core=AddSub_DSP
		  ts_baddr += tlen[i];
		//#pragma HLS resource variable=ts_baddr core=AddSub_DSP
	  }
	  if (aw[0] > aw[1]) width = aw[0];
	  else width = aw[1];

	  pe_matchs.write(seed_index);
	  //result_data = (qBeg & 0xFFFF) | ((qEnd<<16) & 0xFFFF0000);
	  pe_matchs.write((qBeg & 0xFFFF) | ((qEnd<<16) & 0xFFFF0000));
	  //result_data = (rBeg & 0xFFFF) | ((rEnd<<16) & 0xFFFF0000);
	  pe_matchs.write((rBeg & 0xFFFF) | ((rEnd<<16) & 0xFFFF0000));
	  //result_data = (score & 0xFFFF) | ((trueScore<<16) & 0xFFFF0000);
	  pe_matchs.write((score & 0xFFFF) | ((trueScore<<16) & 0xFFFF0000));
	  //result_data = width & 0xFFFF;
	  pe_matchs.write(width & 0xFFFF);
}

void load_data(int *totalinp, int *partInp, int startTaskIndex, int partTaskNums)
{
	//int i=0;
	int staOffset = 0;
	int endOffset = 0;
	int firstTaskPos = 0;
	int lastTaskPos = 0;
	uint8_t qlen[2];
#pragma HLS ARRAY_PARTITION variable=qlen complete dim=0
	uint11_t tlen[2];
#pragma HLS ARRAY_PARTITION variable=tlen complete dim=0
	uint11_t qlen2 = 0;
	uint11_t tlen2 = 0;
	int dataLength = 0;
	int tmp_parame = 0;
	uint11_t qrLen_div8 = 0;
	int headerLength = 0;

	firstTaskPos = (startTaskIndex + 1)*8;
//#pragma HLS resource variable=firstTaskPos core=AddSub_DSP
	memcpy(&partInp[0+3], (const void*)(&totalinp[0]), (8)*4);
	memcpy(&partInp[8+3], (const void*)(&totalinp[firstTaskPos]), (partTaskNums*8)*4);

	headerLength = (partTaskNums+1)*8+3;
	staOffset = totalinp[firstTaskPos + 2];
	lastTaskPos = (startTaskIndex + partTaskNums)*8;
	tmp_parame = totalinp[lastTaskPos + 0];
	qlen[0] = tmp_parame & 0xFFFF; //55
	tlen[0] = (tmp_parame >> 16) & 0xFFFF; //105
	tmp_parame = totalinp[lastTaskPos + 1];
	qlen[1] = tmp_parame & 0xFFFF; //55
	tlen[1] = (tmp_parame >> 16) & 0xFFFF; //105
	endOffset = totalinp[lastTaskPos + 2];
	qlen2 = qlen[0] + qlen[1];
	tlen2 = tlen[0] + tlen[1];
	qrLen_div8 = qlen2+tlen2;
	if ((qrLen_div8 & 0x00000007) != 0) {
		qrLen_div8 = (qrLen_div8 >> 3) + 1;
	}
	else {
		qrLen_div8 = (qrLen_div8 >> 3);
	}

	dataLength = (endOffset - staOffset + qrLen_div8);
//#pragma HLS resource variable=dataLength core=AddSub_DSP
	memcpy(&partInp[headerLength], (const void*)(&totalinp[staOffset]), dataLength*4);
	partInp[0] = headerLength + dataLength;
	partInp[1] = startTaskIndex;
	partInp[2] = partTaskNums;
}

void load_task(int *partInp, int startTaskIndex, int partTaskNums, int staOffset, int i, stream<int>& pe_seeds, stream<uint2_t>& pe_seeds_ctrl)
{
#pragma HLS INLINE
	char k = 0;
	uint11_t l = 0;
	uint8_t qlen[2];
#pragma HLS ARRAY_PARTITION variable=qlen complete dim=0
	uint11_t tlen[2];
#pragma HLS ARRAY_PARTITION variable=tlen complete dim=0
	uint11_t qlen2 = 0;
	uint11_t tlen2 = 0;
	int tmp_parame = 0;
	int endOffset = 0;
	uint11_t qrLen_div8 = 0;
	int taskParaPos = 0;
	int taskDataPos = 0;

	pe_seeds.write(partInp[0]);
	pe_seeds.write(partInp[1]);

	taskParaPos = (i+1)*8;
	tmp_parame = partInp[taskParaPos + 0];
	qlen[0] = tmp_parame & 0xFFFF; //55
	tlen[0] = (tmp_parame >> 16) & 0xFFFF; //105
	tmp_parame = partInp[taskParaPos + 1];
	qlen[1] = tmp_parame & 0xFFFF; //55
	tlen[1] = (tmp_parame >> 16) & 0xFFFF; //105
	tmp_parame = partInp[taskParaPos + 2];
	endOffset = tmp_parame;
	qlen2 = qlen[0] + qlen[1];
	tlen2 = tlen[0] + tlen[1];
	qrLen_div8 = qlen2+tlen2;
	if ((qrLen_div8 & 0x00000007) != 0) {
		qrLen_div8 = (qrLen_div8 >> 3) + 1;
	}
	else {
		qrLen_div8 = (qrLen_div8 >> 3);
	}
	taskDataPos = endOffset - staOffset + (partTaskNums+1)*8;
//#pragma HLS resource variable=taskDataPos core=AddSub_DSP
//	pe_seeds.write(qrLen_div8 + 11);
//	pe_seeds.write(startTaskIndex + i);

	for (k = 0; k < 8; k++) {
#pragma HLS PIPELINE II=1
		pe_seeds.write(partInp[taskParaPos + k]);
	}
	pe_seeds_ctrl.write(1);
	for (l = 0; l < qrLen_div8; l++) {
#pragma HLS PIPELINE II=1
		pe_seeds.write(partInp[taskDataPos + l]);
	}
}

bool feed_peArray(int *blockTaskBuf, stream<int>& peArr_blockTask, stream<uint2_t>& peArr_blockTask_ctrl)
{
#pragma HLS INLINE
	short i=0;
	int blockTaskLen=0;

	if (!peArr_blockTask_ctrl.full()) {
		blockTaskLen = blockTaskBuf[0];
		for (i=0; i<blockTaskLen; i++) {
#pragma HLS PIPELINE II=1
			peArr_blockTask.write(blockTaskBuf[i]);
		}
		peArr_blockTask_ctrl.write(1);
		//peArr_blockTask_ctrl.write(1);
		return true;
	}
	else {
		return false;
	}
}

//bool feed_seed(int *oneTaskBuf, stream<int>& pe_seeds, stream<uint2_t>& pe_seeds_ctrl)
//{
//#pragma HLS INLINE
//	short i=0;
//	int taskLen=0;
//
//	if (!pe_seeds_ctrl.full()) {
//		taskLen = oneTaskBuf[0];
//		for (i=0; i<taskLen+1; i++) {
//			pe_seeds.write(oneTaskBuf[i]);
//		}
//		pe_seeds_ctrl.write(1);
//		//pe_seeds_ctrl.write(1);
//		return true;
//	}
//	else {
//		return false;
//	}
//}

void proc_block(stream<int>& peArr_blockTask, stream<uint2_t>& peArr_blockTask_ctrl, stream<int> seeds[PE_NUMS], stream<uint2_t> seeds_ctrl[PE_NUMS])
{
	int i=0;
	int j=0;
	char m=0;
	uint2_t endFlag0=0;
	//uint2_t endFlag1=0;
	int seedLen=0;
	int partInp[BLKBUF_DEPTH];
	int startTaskIndex=0;
	int partTaskNums=0;

	while(1) {
#if 1
		while (1) {
			if (!peArr_blockTask_ctrl.empty()) {
				break;
			}
		}
#endif
		peArr_blockTask_ctrl.read(endFlag0);
		//peArr_blockTask_ctrl.read(endFlag1);
		if (endFlag0 == 3) {
			for(m=0; m<PE_NUMS; m++) {
				//seeds_ctrl[m].write(3);
				seeds_ctrl[m].write(3);
			}
			break;
		}
		else {
			peArr_blockTask.read(seedLen);
			peArr_blockTask.read(startTaskIndex);
			peArr_blockTask.read(partTaskNums);
			for (j=0; j<seedLen-3; j++) {
#pragma HLS PIPELINE II=1
				peArr_blockTask.read(partInp[j]);
			}
			task_parse(partInp, startTaskIndex, partTaskNums, seeds, seeds_ctrl);
		}
	}
}

void proc_element(stream<int>& pe_seeds, stream<uint2_t>& pe_seeds_ctrl, stream<int>& pe_matchs)
{
//	short i=0;
//	short j=0;
//	char m=0;
	uint2_t endFlag0=0;
	uint2_t endFlag1=0;
//	int seedLen=0;
//	int oneSeedBuf[256];
//	int matchsBuf[5];
//#pragma HLS RESOURCE variable=matchsBuf core=RAM_1P_LUTRAM

	while(1) {
#if 1
		while (1) {
			if (!pe_seeds_ctrl.empty()) {
				break;
			}
		}
#endif
		pe_seeds_ctrl.read(endFlag0);
		//pe_seeds_ctrl.read(endFlag1);
		if (endFlag0 == 3) {
			pe_matchs.write(0xFFFFFFFF);
			break;
		}
//		pe_seeds.read(seedLen);
//		for (j=0; j<seedLen; j++) {
//			pe_seeds.read(oneSeedBuf[j]);
//		}
//		leftright_ext(oneSeedBuf, matchsBuf);
		leftright_ext(pe_seeds, pe_matchs);
//		write_match(matchsBuf, pe_matchs);
	}
}

//void write_match(int *oneMatchBuf, stream<int>& pe_matchs)
//{
//#pragma HLS INLINE
//	char i=0;
//
//	for (i=0; i<5; i++) {
//#pragma HLS PIPELINE II=1
//		pe_matchs.write(oneMatchBuf[i]);
//	}
//}

void task_parse(int *partInp, int startTaskIndex, int partTaskNums, stream<int> seeds[PE_NUMS], stream<uint2_t> seeds_ctrl[PE_NUMS])
{
#pragma HLS INLINE
	int i = 0;
	char j = 0;
	char beg = 0;
	char m = 0;
	int staOffset = 0;
	int oneTaskBuf[256];

	staOffset = partInp[1*8 + 2];
	for(i=0; i<partTaskNums; i++) {
//		load_task(partInp, startTaskIndex, partTaskNums, staOffset, i, oneTaskBuf);
#if 1
		while(1) {
//			if (feed_seed(oneTaskBuf, seeds[j], seeds_ctrl[j])) {
			if (!seeds_ctrl[j].full()) {
				load_task(partInp, startTaskIndex, partTaskNums, staOffset, i, seeds[j], seeds_ctrl[j]);
				break;
			}
			else {
               j++;
               if (j >= PE_NUMS) j = 0;
			}
		}
#else
     static int jj = 0;

    load_task(partInp, startTaskIndex, partTaskNums, staOffset, i, seeds[jj], seeds_ctrl[jj]);
    jj++;
    if (jj >= PE_NUMS) jj=0;
#endif
	}
}

void feed_data(char mode, int *blockTaskMem, stream<int> blockTask[PEARRAY_NUM], stream<uint2_t> blockTask_ctrl[PEARRAY_NUM])
{
//#pragma HLS INLINE
	char j;
	char k;

	j = 0;
	while(1) {
		if (feed_peArray(blockTaskMem, blockTask[j], blockTask_ctrl[j])) {
			break;
		}
		else {
            j++;
            if (j >= PEARRAY_NUM) j = 0;
		}
	}

	if (mode == 1) {
	    for (k=0; k<PEARRAY_NUM; k++) {
	    	//blockTask_ctrl[k].write(3);
	    	blockTask_ctrl[k].write(3);
	    }
	}
}
void data_parse(int *totalinp, int __inc, stream<int> blockTask[PEARRAY_NUM], stream<uint2_t> blockTask_ctrl[PEARRAY_NUM])
{
//#pragma HLS INLINE
	int i=0;
	int partNums;
	int blockTaskMem_A[BLKBUF_DEPTH];
	int blockTaskMem_B[BLKBUF_DEPTH];
	int startTaskIndex;
	int partTaskNums;
	int last_partTaskNums;
	char mode;

	//if (__inc & (PART_TASK_NUMS-1)) {
	//	partNums = (__inc >> LOG2_PART_TASK_NUMS) + 1;
	//}
	//else {
	//	partNums = (__inc >> LOG2_PART_TASK_NUMS);
	//}
	partNums = (__inc + PART_TASK_NUMS -1) / PART_TASK_NUMS;
	last_partTaskNums = __inc % PART_TASK_NUMS;
	startTaskIndex = 0;
	mode = 0;
data_parse_loop :
    for (i=0; i<partNums+1; i++) {
		if (i == (partNums - 1)) {
			//if ((__inc & (PART_TASK_NUMS-1)) == 0) {
			//	partTaskNums = PART_TASK_NUMS;
			//}
			//else {
			//	partTaskNums = __inc & (PART_TASK_NUMS-1);
			//}
			if (last_partTaskNums == 0) {
				partTaskNums = PART_TASK_NUMS;
			}
			else {
				partTaskNums = last_partTaskNums;
			}
		}
		else {
			partTaskNums = PART_TASK_NUMS;
		}
		if (i%2 == 0) {
			if (i == 0) {
				load_data(totalinp, blockTaskMem_A, startTaskIndex, partTaskNums);
			}
			else if (i == partNums) {
				feed_data(1, blockTaskMem_B, blockTask, blockTask_ctrl);
			}
			else {
				load_data(totalinp, blockTaskMem_A, startTaskIndex, partTaskNums);
				feed_data(0, blockTaskMem_B, blockTask, blockTask_ctrl);
			}
		}
		else {
			if (i == partNums) {
				feed_data(1, blockTaskMem_A, blockTask, blockTask_ctrl);
			}
			else {
				load_data(totalinp, blockTaskMem_B, startTaskIndex, partTaskNums);
				feed_data(0, blockTaskMem_A, blockTask, blockTask_ctrl);
			}
		}
		startTaskIndex += partTaskNums;
		//#pragma HLS resource variable=startTaskIndex core=AddSub_DSP
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

void receive_match(stream<int> matchs[PE_NUMS], stream<int>& peArr_blockMatchs)
{
//	#pragma HLS INLINE
	char l=0;
	char j=0;
	char getPeMatchFlag = 0;
	int matchBuf[5];
#pragma HLS RESOURCE variable=matchBuf core=RAM_1P_LUTRAM
	char peOver_cnt=0;

	peOver_cnt = 0;
	j = 0;
	while(1) {
		getPeMatchFlag = get_peMatch(matchs[j], matchBuf);
		if (getPeMatchFlag == 0) {
            j++;
            if (j >= PE_NUMS) j = 0;
		}
		else if (getPeMatchFlag == 1) {
			for (l=0; l<5; l++) {
#pragma HLS PIPELINE II=1
				peArr_blockMatchs.write(matchBuf[l]);
			}
		}
		else {
			peOver_cnt = peOver_cnt + 1;
		}
		if (peOver_cnt >= PE_NUMS) {
			peArr_blockMatchs.write(0xFFFFFFFF);
			break;
		}
	}
}

char get_peArrMatch(stream<int>& peArr_blockMatchs, int *results)
{
#pragma HLS INLINE
	char i=0;
	int seed_index=0;
	int resul_addr;

	if (!peArr_blockMatchs.empty()) {
		peArr_blockMatchs.read(seed_index);
		if (seed_index == 0xFFFFFFFF) {
			return 2;
		}
		else {
			results[0] = seed_index;
			for (i = 1; i < 5; i++) {
#pragma HLS PIPELINE II=1
			//	resul_addr = seed_index*4 + i;
			//#pragma HLS resource variable=resul_addr core=AddSub_DSP
				peArr_blockMatchs.read(results[i]);
			}
			return 1;
		}
	}
	else {
		return 0;
	}
}

void results_assemble(stream<int> blockMatchs[PEARRAY_NUM], stream<int>& results)
{
	char j=0;
	char k=0;
	char getPeArrMatchFlag = 0;
	char peArrOver_cnt=0;
	int seed_index=0;
	int oneMatchBuf[5];

	peArrOver_cnt = 0;
	while(1) {
		getPeArrMatchFlag = get_peArrMatch(blockMatchs[j], oneMatchBuf);
		if (getPeArrMatchFlag == 0) {
            j++;
            if (j >= PEARRAY_NUM) j = 0;
		}
		else if (getPeArrMatchFlag == 1) {
			for (k = 0; k < 5; k++) {
#pragma HLS PIPELINE II=1
			//	resul_addr = seed_index*4 + i;
			//#pragma HLS resource variable=resul_addr core=AddSub_DSP
				results.write(oneMatchBuf[k]);
			}
		}
		else if (getPeArrMatchFlag == 2) {
			peArrOver_cnt = peArrOver_cnt + 1;
		}
		if (peArrOver_cnt >= PEARRAY_NUM) {
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

	outAddr = 0;
	while(1) {
		locAddr = 0;
		flag_over = fill_resulBuf(results, resultsBuf, &locAddr);
		if (flag_over == 2) {
			break;
		}
		else if (flag_over == 1) {
			memcpy(&output_a[outAddr], resultsBuf, locAddr*4);
			break;
		}
		else if (flag_over == 0) {
			memcpy(&output_a[outAddr], resultsBuf, locAddr*4);
			outAddr = outAddr + locAddr;
		}
	}
}

void sw_core(int *a, int *output_a, int __inc)
{
#pragma HLS inline off
	int local_inc=0;
    const int EXTRA = 0;

	local_inc = __inc;

#pragma HLS DATAFLOW
		char i, j;
		stream<int> blockTask[PEARRAY_NUM];
	#pragma HLS STREAM variable=blockTask depth=4096+EXTRA
		stream<uint2_t> blockTask_ctrl[PEARRAY_NUM];
	#pragma HLS STREAM variable=blockTask_ctrl depth=1+EXTRA

		stream<int> seeds[PEARRAY_NUM][PE_NUMS];
#pragma HLS STREAM variable=seeds depth=256+EXTRA
		stream<uint2_t> seeds_ctrl[PEARRAY_NUM][PE_NUMS];
#pragma HLS STREAM variable=seeds_ctrl depth=1+EXTRA
		stream<int> matchs[PEARRAY_NUM][PE_NUMS];
#pragma HLS STREAM variable=matchs depth=128+EXTRA

		stream<int> blockMatchs[PEARRAY_NUM];
	#pragma HLS STREAM variable=blockMatchs depth=512+EXTRA

		stream<int> resul_out("resul_out");
	#pragma HLS STREAM variable=resul_out depth=2048+EXTRA

		data_parse(a, local_inc, blockTask, blockTask_ctrl);
		for (i=0; i<PEARRAY_NUM; i++) {
#pragma HLS UNROLL
			proc_block(blockTask[i], blockTask_ctrl[i], seeds[i], seeds_ctrl[i]);
			for (j=0; j<PE_NUMS; j++) {
#pragma HLS UNROLL
				proc_element(seeds[i][j], seeds_ctrl[i][j], matchs[i][j]);
			}
			receive_match(matchs[i], blockMatchs[i]);
		}

		results_assemble(blockMatchs, resul_out);
		upload_results(resul_out, output_a);
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
