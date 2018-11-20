#ifndef SAMBLASTER
#define SAMBLASTER
#define __STDC_FORMAT_MACROS

#include <stdlib.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <map>
#include "sbhash.h"

typedef uint64_t UINT64;
typedef uint32_t UINT32;
typedef int32_t   INT32;

// mempcpy is a GNU extension and not available everywhere.
#ifndef _GNU_SOURCE
inline void *mempcpy(void *dest, const void *src, size_t n) {
    return (char*) memcpy(dest, src, n) + n; 
}
#endif
inline bool streq(char * s1, const char * s2) __attribute__((always_inline));
inline bool streq(char * s1, const char * s2) {
    return (strcmp(s1, s2) ==0);
}

void fatalError(const char* errorStr);
void fsError(const char * filename);

// Stuff needed for timings.
// // To turn timing off, set the below to 0.
#define TIMING 1
#if TIMING

void fprintTimeSeconds (FILE * out, double seconds, int precision);
void fprintTimeMicroSeconds (FILE * out, UINT64 microSeconds, int precision);

inline UINT64 diffTVs (struct timeval * startTV, struct timeval * endTV) {
    return (((endTV->tv_sec - startTV->tv_sec) * 1000000) + (endTV->tv_usec - startTV->tv_usec));
}

#include <sys/times.h>
#include <sys/resource.h>
#include <time.h>
#endif // TIMING

typedef INT32 pos_t;
#define MAX_SEQUENCE_LENGTH 250 // Current illumina paired-end reads are at most 150 + 150
inline int padLength(int length) {
    return length + (2 * MAX_SEQUENCE_LENGTH);
}
inline int padPos(int pos) {
    return pos + MAX_SEQUENCE_LENGTH;
}

typedef UINT64 sgn_t; // Type for signatures for offsets and lengths.
typedef struct splitLine splitLine_t;
extern splitLine_t * splitLineFreeList;
struct splitLine {
    // General fields for any split line.
    splitLine_t * next;
    char * buffer;
    int bufLen;
    size_t maxBufLen;
    char **fields;
    int numFields;
    int maxFields;
    bool split;
    // Special SAM fields that we need to access as other than strings.
    // It this were a class, these would be in a subclass.
    int   flag;
    pos_t pos;
    int   seqNum;
    pos_t binPos;
    int   binNum;
    int   SQO;
    int   EQO;
    int   sclip;
    int   eclip;
    int   rapos;
    int   raLen;
    int   qaLen;
    bool  CIGARprocessed;
    bool  discordant;
    bool  splitter;
    bool  unmappedClipped;
};

splitLine_t * makeSplitLine();
void deleteSplitLine(splitLine_t * line);
void cleanUpSplitLines();
void disposeSplitLines(splitLine_t * line);
splitLine_t * getSplitLine();
void splitSplitLine(splitLine_t * line, int maxSplits);
void unsplitSplitLine(splitLine_t * line);
void resizeSplitLine(splitLine_t * line, int newsize);
void changeFieldSplitLine(splitLine_t * line, int fnum, char * newValue);
void addTag(splitLine_t * line, const char * header, const char * val);
splitLine_t * readLine(FILE * input);

inline void outputString(char * str, FILE * output){
  // Do the error checking here so we don't have to do it elsewhere.
  if (fputs(str, output) < 0) {
        fatalError("samblaster: Unable to write to output file.\n");
    }
}

inline void writeLine(splitLine_t * line, FILE * output) {
    unsplitSplitLine(line);
    outputString(line->buffer, output);
}

void checkBAMfile(splitLine_t * line);

// Define SAM field offsets.
#define QNAME  0
#define FLAG   1
#define RNAME  2
#define POS    3
#define MAPQ   4
#define CIGAR  5
#define RNEXT  6
#define PNEXT  7
#define TLEN   8
#define SEQ    9
#define QUAL  10
#define TAGS  11

// Define SAM flag accessors
#define MULTI_SEGS     0x1
#define FIRST_SEG      0x40
#define SECOND_SEG     0x80
inline bool checkFlag(splitLine_t * line, int bits) { return ((line->flag & bits) != 0); }

inline void setFlag(splitLine_t * line, int bits) { line->flag |= bits; }

inline bool isPaired(splitLine_t * line) { return checkFlag(line, MULTI_SEGS); }

inline bool isConcordant(splitLine_t * line) { return checkFlag(line, 0x2); }

inline bool isDiscordant(splitLine_t * line) { return !isConcordant(line); }

inline bool isUnmapped(splitLine_t * line) { return checkFlag(line, 0x4); }

inline bool isNextUnmapped(splitLine_t * line) { return checkFlag(line, 0x8); }

inline bool isNextMapped(splitLine_t * line) { return !isNextUnmapped(line); }

inline bool isMapped(splitLine_t * line) { return !isUnmapped(line); }

inline bool isReverseStrand(splitLine_t * line) { return checkFlag(line, 0x10); }

inline bool isForwardStrand(splitLine_t * line) { return !isReverseStrand(line); }

inline bool isFirstRead(splitLine_t * line) { return checkFlag(line, FIRST_SEG); }

inline bool isSecondRead(splitLine_t * line) { return checkFlag(line, SECOND_SEG); }

// These determine alignment type.
// // Things may get more complicated than this once we have alternate contigs such as in build 38 of human genome.
inline bool isPrimaryAlignment(splitLine_t * line) { 
  return !(checkFlag(line, 0x100) || checkFlag(line, 0x800)); 
}

// We have to hande secondard and complementary alignments differently depending on compatMode.
// // So, we store which bits are being included in each.
extern int complementaryBits;
inline bool isComplementaryAlignment(splitLine_t * line) { 
  return checkFlag(line, complementaryBits);
}

extern int secondaryBits;
inline bool isSecondaryAlignment(splitLine_t * line)  {
  return checkFlag(line, secondaryBits); 
}

inline bool isDuplicate(splitLine_t * line) { return checkFlag(line, 0x400); }

inline void setDuplicate(splitLine_t * line) { setFlag(line, 0x400); }

typedef hashTable_t sigSet_t;

inline int str2int (char * str) {
    return strtol(str, NULL, 0);
}
// Need to change this if pos is unsigned.
inline pos_t str2pos (char * str) {
    return strtol(str, NULL, 0);
}
// Temp buffer to use to form new flag field when marking dups.
extern char tempBuf[10];
inline void markDup(splitLine_t * line) {
    setDuplicate(line);
    sprintf(tempBuf, "%d", line->flag);
    changeFieldSplitLine(line, FLAG, tempBuf);
}
// Special version of write line that appends an id number to the output.
// // Used to output splitters.
void writeSAMlineWithIdNum(splitLine_t * line, FILE * output);


///////////////////////////////////////////////////////////////////////////////
//// Sequence Map
/////////////////////////////////////////////////////////////////////////////////
//
//// We use a map instead of a hash map for sequence names.
//// This is because the default hash function on char * hashes the ptr values.
//// So, we would need to define our own hash on char * to get things to work properly.
//// Not worth it for a structure holding so few members.
//
//// Function needed to get char * map to work.
struct less_str {
   bool operator()(char const *a, char const *b) const
   {
      return strcmp(a, b) < 0;
   }
};
// This stores the map between sequence names and sequence numbers.
typedef std::map<const char *, int, less_str> seqMap_t;

inline void addSeq(seqMap_t * seqs, char * item, int val) {
    (*seqs)[item] = val;
}

///////////////////////////////////////////////////////////////////////////////
//// Struct for processing state
/////////////////////////////////////////////////////////////////////////////////
struct state_struct {
    char *         inputFileName;
    FILE *         inputFile;
    char *         outputFileName;
    FILE *         outputFile;
    FILE *         discordantFile;
    char *         discordantFileName;
    FILE *         splitterFile;
    char *         splitterFileName;
    FILE *         unmappedClippedFile;
    char *         unmappedClippedFileName;
    sigSet_t *     sigs;
    seqMap_t       seqs;
    UINT32 *       seqLens;
    UINT64 *       seqOffs;
    splitLine_t ** splitterArray;
    int            splitterArrayMaxSize;
    UINT32         sigArraySize;
    int            binCount;
    int            minNonOverlap;
    int            maxSplitCount;
    int            minIndelSize;
    int            maxUnmappedBases;
    int            minClip;
    int            unmappedFastq;
    bool           acceptDups;
    bool           excludeDups;
    bool           removeDups;
    bool           addMateTags;
    bool           compatMode;
    bool           ignoreUnmated;
    bool           quiet;
};
typedef struct state_struct state_t;
state_t * makeState ();
void deleteState(state_t * s);

///////////////////////////////////////////////////////////////////////////////
//// Signatures
/////////////////////////////////////////////////////////////////////////////////
//// We now calculate signatures as offsets into a super contig that includes the entire genome.
//// And partition it into equaly sized bins.
//// This performs better on genomes with large numbers of small contigs,
////  without performance degradation on more standard genomes.
//// Thanks to https://github.com/carsonhh for the suggestion.
///////////////////////////////////////////////////////////////////////////////
inline sgn_t calcSig(splitLine_t * first, splitLine_t * second) {
    // Total nonsense to get the compiler to actually work.
    UINT64 t1 = first->binPos;
    UINT64 t2 = t1 << 32;
    UINT64 final = t2 | second->binPos;
    return (sgn_t)final;
}

inline UINT32 calcSigArrOff(splitLine_t * first, splitLine_t * second, int binCount) {
    UINT32 s1 = (first->binNum * 2) + (isReverseStrand(first) ? 1 : 0);
    UINT32 s2 = (second->binNum * 2) + (isReverseStrand(second) ? 1 : 0);
    UINT32 retval = (s1 * binCount * 2) + s2;
#ifdef DEBUG
    fprintf(stderr, "1st %d %d -> %d 2nd %d %d -> %d count %d result %d read: %s\n",
            first->binNum, isReverseStrand(first), s1, second->binNum, isReverseStrand(second), s2, binCount, retval, first->fields[QNAME]);
#endif
    return retval;
}

///////////////////////////////////////////////////////////////////////////////
//// Sequences
/////////////////////////////////////////////////////////////////////////////////
inline int getSeqNum(splitLine_t * line, int field, state_t * state) __attribute__((always_inline));
inline int getSeqNum(splitLine_t * line, int field, state_t * state) {
#ifdef DEBUG
    seqMap_t::iterator findret = state->seqs.find(line->fields[field]);
    if (findret == state->seqs.end()) {
        char * temp;
        asprintf(&temp, "Unable to find seq %s for readid %s in sequence map.\n", line->fields[field], line->fields[QNAME]);
        fatalError(temp);
    }
    return findret->second;
#else
    return state->seqs.find(line->fields[field])->second;
#endif
}

///////////////////////////////////////////////////////////////////////////////
//// Helpers to process CIGAR strings
/////////////////////////////////////////////////////////////////////////////////

// This will parse a base 10 int, and change ptr to one char beyond the end of the number.
inline int parseNextInt(char **ptr) {
    int num = 0;
    for (char curChar = (*ptr)[0]; curChar != 0; curChar = (++(*ptr))[0]) {
        int digit = curChar - '0';
        if (digit >= 0 && digit <= 9) num = num*10 + digit;
        else break;
    }
    return num;
}
// This will the current char, and move the ptr ahead by one.
inline char parseNextOpCode(char **ptr) {
    return ((*ptr)++)[0];
}
// This just test for end of string.
inline bool moreCigarOps(char *ptr) {
    return (ptr[0] != 0);
}
void calcOffsets(splitLine_t * line);

inline int getStartDiag(splitLine_t * line) {
  // SRO - SQO (not strand normalized)
  // Simplify the following.
  // return (str2pos(line->fields[POS])) - line->sclip;
  return line->rapos - line->sclip;
}

inline int getEndDiag(splitLine_t * line) {
    // ERO - EQO (not strand normalized)
    // Simplify the following
    // return (line->rapos + line->raLen - 1) - (line->sclip + line->qaLen - 1)
    return (line->rapos + line->raLen) - (line->sclip + line->qaLen);
}

///////////////////////////////////////////////////////////////////////////////
//// Process SAM Blocks
/////////////////////////////////////////////////////////////////////////////////
#define BIN_SHIFT 27 //bin window is 27 bits wide
#define BIN_MASK ((1 << 27)-1) //bin window is 27 bits wide

void outputSAMBlock(splitLine_t * block, FILE * output);

inline bool needSwap(splitLine_t * first, splitLine_t * second) {
    // Sort first by ref offset.
    if (first->pos > second->pos) return true;
    if (first->pos < second->pos) return false;
    // Now by seq number.
    if (first->seqNum > second->seqNum) return true;
    if (first->seqNum < second->seqNum) return false;
    // Now by strand.
    // If they are both the same strand, it makes no difference which on is first.
    if (isReverseStrand(first) == isReverseStrand(second)) return false;
    if (isReverseStrand(first) && isForwardStrand(second)) return true;
    return false;
}

inline void swapPtrs(splitLine_t ** first, splitLine_t ** second) {
    splitLine_t * temp = *first;
    *first = *second;
    *second = temp;
}

void brokenBlock(splitLine_t *block, int count);

template <bool excludeSecondaries>
int fillSplitterArray(splitLine_t * block, state_t * state, int mask, bool flagValue);
void markDupsDiscordants(splitLine_t * block, state_t * state);

int compQOs(const void * p1, const void * p2);

template <bool excludeSecondaries>
int fillSplitterArray(splitLine_t * block, state_t * state, int mask, bool flagValue);
void markSplitterUnmappedClipped(splitLine_t * block, state_t * state, int mask, bool flagValue);

void writeUnmappedClipped(splitLine_t * line, state_t * state);
void processSAMBlock(splitLine_t * block, state_t * state);
#else
#endif
