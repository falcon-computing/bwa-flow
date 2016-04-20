#ifndef SWTASK_H
#define SWTASK_H

#include "bwa_wrapper.h"

class SWRead;

class ExtParam
{
  public:
    uint8_t *leftQs;
    int     leftQlen ;
    uint8_t *leftRs;
    int     leftRlen ;
    uint8_t *rightQs;
    int     rightQlen ;
    uint8_t *rightRs;
    int     rightRlen ;
    int     w ;
    int8_t *mat;
    int     oDel ;
    int     eDel ;
    int     oIns ;
    int     eIns ;
    int     penClip5 ;
    int     penClip3 ;
    int     zdrop   ;
    int     h0 ;
    int     regScore ;
    int     qBeg    ;
    int     idx    ;
    int64_t     rBeg   ;    // just for testing on cpu
    int     seedLength ;
    int     l_query ;

    SWRead* read_obj;
    mem_alnreg_t* newreg;
    int chain_idx;
    mem_chain_t* chain;
};

#endif
