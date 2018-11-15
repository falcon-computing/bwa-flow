// This define is needed for portable definition of PRIu64
#define __STDC_FORMAT_MACROS

#include <stdlib.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <map>
#include "samblaster.h"

int secondaryBits = 0x100;
int complementaryBits = 0x800;
splitLine_t* splitLineFreeList = NULL;

void fatalError(const char * errorStr) {
    fprintf(stderr, "%s\n", errorStr);
    exit(1);
}

void fsError(const char * filename) {
    char * temp;
    if (errno == ENOENT)
        asprintf(&temp, "samblaster: File '%s' does not exist.\n", filename);
    else
        asprintf(&temp, "samblaster: File system error on %s: %d.\n", filename, errno);
    fatalError(temp);
}

///////////////////////////////////////////////////////////////////////////////
// Runtime Statistics
///////////////////////////////////////////////////////////////////////////////

// A convenience function for outputing time is seconds in a more useful metric.
void fprintTimeSeconds (FILE * out, double seconds, int precision) {
    double totalseconds = seconds;
    int hours = seconds/3600.;
    if (hours > 0) {
        seconds -= hours * 3600;
        fprintf(out, "%dH", hours);
    }
    int minutes = seconds/60.;
    if (minutes > 0) {
        seconds -= minutes * 60;
        fprintf(out, "%dM", minutes);
    }
    if (hours + minutes > 0) {
        fprintf(out, "%.0fS", seconds);
        fprintf(out, "(%.*fS)", precision, totalseconds);
    }
    else fprintf(out, "%.*fS", precision, totalseconds);
}

void fprintTimeMicroSeconds (FILE * out, UINT64 microSeconds, int precision) {
    fprintTimeSeconds(out, ((double)microSeconds/1000000.0), precision);
}

///////////////////////////////////////////////////////////////////////////////
// Split Lines
///////////////////////////////////////////////////////////////////////////////

// The structure to store "split" input lines, especially SAM lines.
// They form a singly linked list so that we can form groups of them,
//     and also so that we can keep a freelist of them.

// Reference offsets (pos) can be negative or run off the end of a contig due to clipping.
// Therefore, we will use a padding strategy.
// The space allocated to each contig will be padded by twice the max read length.
// This leaves room for both offset underflow and overflow.
// And all offsets will be shifted higher by the max read length.
// This will eliminate negative offsets and "center" offsets within the offset range for the contig.

// Creator for splitLine
splitLine_t * makeSplitLine() {
    splitLine_t * line = (splitLine_t *)malloc(sizeof(splitLine_t));
    line->bufLen = 0;
    line->maxBufLen = 1500;
    line->buffer = (char *)malloc(line->maxBufLen);
    line->numFields = 0;
    line->maxFields = 100;
    line->fields = (char **)malloc(line->maxFields * sizeof(char *));
    return line;
}

// Destructor for split line.
void deleteSplitLine(splitLine_t * line) {
    free(line->buffer);
    free(line->fields);
    free(line);
}

// Descructor for a list of splitLines
void cleanUpSplitLines() {
    splitLine_t * l = splitLineFreeList;
    while (l != NULL)
    {
        splitLine_t * next = l->next;
        deleteSplitLine(l);
        l = next;
    }
}

// Like descructor for splitLine except don't free memory.
// Instead, put the linked list of objects back on the free list.
void disposeSplitLines(splitLine_t * line) {
    // First find the last line in the list.
    // Then get rid of them all.
    splitLine_t * last = line;
    for (splitLine_t * l = line->next; l!=NULL; l = l->next) last = l;
    last->next = splitLineFreeList;
    splitLineFreeList = line;
}

// Like constuctor, except take struct off free list if there is one.
splitLine_t * getSplitLine() {
    splitLine_t * line;
    if (splitLineFreeList ==  NULL) {
        line = makeSplitLine();
    }
    else {
        line = splitLineFreeList;
        splitLineFreeList = splitLineFreeList->next;
    }
    line->next = NULL;
    line->CIGARprocessed = false;
    // Mark these here so that the code for these doesn't have to stand on its head to do it.
    line->discordant = false;
    line->splitter = false;
    line->unmappedClipped = false;
    return line;
}

// Split the line into fields.
void splitSplitLine(splitLine_t * line, int maxSplits) {
    line->numFields = 0;
    int fieldStart = 0;
    // replace the newline with a tab so that it works like the rest of the fields.
    line->buffer[line->bufLen-1] = '\t';
    for (int i=0; i<line->bufLen; ++i) {
        if (line->buffer[i] == '\t') {
            line->fields[line->numFields] = line->buffer + fieldStart;
            line->numFields += 1;
            if (line->numFields == maxSplits) break;
            line->buffer[i] = 0;
            // Get ready for the next iteration.
            fieldStart = i+1;
        }
    }
    // replace the tab at the end of the line with a null char to terminate the final string.
    line->buffer[line->bufLen-1] = 0;
    line->split = true;
}

// Unsplit the fields back into a single string.
// This will mung the strings, so only call this when all processing on the line is done.
void unsplitSplitLine(splitLine_t * line) {
    // First make sure we are still split.
    if (!line->split) return;
    // First undo the splits.
    // We will undo the splits backwards from the next field to avoid having to calculate strlen each time.
    for (int i=1; i<line->numFields; ++i) {
        line->fields[i][-1] = '\t';
    }
    // Now put the newline back in.
    line->buffer[line->bufLen-1] = '\n';
    // Mark as no longer split.
    line->split = false;
}

// Resize the buffer of a splitLine.
// Since the new buffer may not be in the same place, we need to first unsplit, resize, then resplit.
void resizeSplitLine(splitLine_t * line, int newsize) {
    // First unsplit it.
    unsplitSplitLine(line);
    // Resize the buffer, giving a little extra room.
    line->maxBufLen = newsize + 50;
    line->buffer = (char *)realloc(line->buffer, line->maxBufLen);
    if (line->buffer == NULL) {
        fatalError("samblaster: Failed to reallocate to a larger read buffer size.\n");
    }
    // Now resplit the line.
    splitSplitLine(line, line->numFields);
}


// Change a field into a value.
// This will be tough given how we output lines.
// So, we might have to try a few things.
// Start with simply expanding/contracting the string to put in the new value.
void changeFieldSplitLine(splitLine_t * line, int fnum, char * newValue) {
    // What we do will depend on the lengths of the two strings.
    // So, start by calculaing these only once.
    char * fp = line->fields[fnum];
    int oldLen = strlen(fp);
    int newLen = strlen(newValue);
    // Now see if we need to first change the length of the whole line.
    int move = newLen - oldLen;
    if (move != 0) {
        // This should never happen, but to be robust we need to check.
        // It is messy to fix it, as all the field ptrs will now be wrong.
        // For now, punt.
        if ((size_t)(line->bufLen + move) >= line->maxBufLen) {
            resizeSplitLine(line, line->bufLen + move);
            fp = line->fields[fnum];
        }
        // Calculate the size of the tail that is still needed.
        int distance = 1 + line->bufLen - (fp - line->buffer) - oldLen;
        // Do the copy.
        memmove(fp+newLen, fp+oldLen, distance);
        // Correct the total length of the buffer.
        line->bufLen += move;
        // We need to correct the other ptrs as well.
        for (int i=fnum+1; i<line->numFields; i++) line->fields[i] += move;
    }
    // Copy in the new value.
    memcpy(fp, newValue, newLen);
}

void addTag(splitLine_t * line, const char * header, const char * val) {
    int hl = strlen(header);
    int vl = strlen(val);
    // Make sure everything will fit.
    int newlen = line->bufLen + hl + vl;
    if ((size_t)newlen >= line->maxBufLen) {
        resizeSplitLine(line, newlen);
    }
    // Copy over the header and the value.
    char * ptr = line->buffer + line->bufLen - 1;
    ptr = (char *)mempcpy(ptr, header, hl);
    ptr = (char *)mempcpy(ptr, val, vl);
    // Add the null terminator for the field, and for the record.
    ptr[0] = 0;
    ptr[1] = 0;
    // Fix the buffer length.
    line->bufLen = newlen;
}


// Read a line from the file and split it.
splitLine_t * readLine(FILE * input) {
    splitLine_t * sline = getSplitLine();
    sline->bufLen = getline(&sline->buffer, &sline->maxBufLen, input);
    if (sline->bufLen < 1) {
        disposeSplitLines(sline);
        return NULL;
    }
    splitSplitLine(sline, 12);
    return sline;
}

// Check the first line of a file (e.g. input) for bam signature.
void checkBAMfile(splitLine_t * line) {
    // If the file is a bam file, we can't rely on fields.
    // So, look at the underlying buffer for the line.

    // First define the signature to look for.
    int values[] = {31, -117, 8, 4, 66, 67,  2};
    int offsets[] = {0,    1, 2, 3, 12, 13, 14};
    int count = 7;

    // Check for empty file or likely a sam file with a header.
    // This is necessary, as the @HD line may not be long enough to check for the BAM signature.
    if (line == NULL) fatalError("samblaster: Input file is empty. Exiting.\n");
    if (line->buffer[0] == '@') return;

    // If a SAM file has no header, an alignment row should easily be long enough.
    if (line->bufLen <= offsets[count-1]) fatalError("samblaster: Input file is empty. Exiting.\n");
    // Check for the BAM signature.
    for (int i=0; i<count; i++) {
        if ((int)line->buffer[offsets[i]] != values[i]) return;
    }

    // If we are here, we almost certainly have a bamfile.
    fatalError("samblaster: Input file appears to be in BAM format. SAM input is required. Exiting.\n");
}

///////////////////////////////////////////////////////////////////////////////
// SAM and signature set related structures.
///////////////////////////////////////////////////////////////////////////////
char tempBuf[10];

void writeSAMlineWithIdNum(splitLine_t * line, FILE * output) {
    // Unsplit the line.
    unsplitSplitLine(line);
    // Split it two ways to isolate the id field.
    splitSplitLine(line, 2);
    outputString(line->fields[0], output);
    if (isPaired(line)) fprintf(output, "_%d\t", isFirstRead(line) ? 1 : 2);
    else                fprintf(output, "\t");
    outputString(line->fields[1], output);
    fprintf(output, "\n");
}

//// Struct for processing state

state_t * makeState () {
    state_t * s = new state_t();
    s->inputFile = stdin;
    s->inputFileName = (char *)"stdin";
    s->outputFile = stdout;
    s->outputFileName = (char *)"stdout";
    s->discordantFile = NULL;
    s->discordantFileName = (char *)"";
    s->splitterFile = NULL;
    s->splitterFileName = (char *)"";
    s->unmappedClippedFile = NULL;
    s->unmappedClippedFileName = (char *)"";
    s->sigs = NULL;
    s->minNonOverlap = 20;
    s->maxSplitCount = 2;
    s->minIndelSize = 50;
    s->maxUnmappedBases = 50;
    s->minClip = 20;
    s->acceptDups = false;
    s->excludeDups = false;
    s->removeDups = false;
    s->addMateTags = false;
    s->compatMode = false;
    s->ignoreUnmated = false;
    s->quiet = false;
    // Start this as -1 to indicate we don't know yet.
    // Once we are outputting our first line, we will decide.
    s->unmappedFastq = -1;
    // Used as a temporary location for ptrs to splitter for sort routine.
    s->splitterArrayMaxSize = 1000;
    s->splitterArray = (splitLine_t **)(malloc(s->splitterArrayMaxSize * sizeof(splitLine_t *)));
    return s;
}

void deleteState(state_t * s) {
    free(s->splitterArray);
    if (s->sigs != NULL) {
        // delete[] s->sigs;
        for (UINT32 i=0; i<s->sigArraySize; i++) s->sigs[i].deleteHashTable();
        free (s->sigs);
    }
    for (seqMap_t::iterator iter = s->seqs.begin(); iter != s->seqs.end(); ++iter) {
        free((char *)(iter->first));
    }
    if (s->seqLens != NULL) free(s->seqLens);
    if (s->seqOffs != NULL) free(s->seqOffs);
    delete s;
}

///////////////////////////////////////////////////////////////////////////////
// Helpers to process CIGAR strings
///////////////////////////////////////////////////////////////////////////////

void calcOffsets(splitLine_t * line) {
    if (line->CIGARprocessed) return;
    char * cigar = line->fields[CIGAR];
    line->raLen = 0;
    line->qaLen = 0;
    line->sclip = 0;
    line->eclip = 0;
    bool first = true;
    while (moreCigarOps(cigar)) {
        int opLen = parseNextInt(&cigar);
        char opCode = parseNextOpCode(&cigar);
        if      (opCode == 'M' || opCode == '=' || opCode == 'X') {
            line->raLen += opLen;
            line->qaLen += opLen;
            first = false;
        }
        else if (opCode == 'S' || opCode == 'H') {
            if (first) line->sclip += opLen;
            else       line->eclip += opLen;
        }
        else if (opCode == 'D' || opCode == 'N') {
            line->raLen += opLen;
        }
        else if (opCode == 'I') {
            line->qaLen += opLen;
        }
        else {
            fprintf(stderr, "Unknown opcode '%c' in CIGAR string: '%s'\n", opCode, line->fields[CIGAR]);
        }
    }
    line->rapos = str2pos(line->fields[POS]);
    if (isForwardStrand(line)) {
        line->pos = line->rapos - line->sclip;
        line->SQO = line->sclip;
        line->EQO = line->sclip + line->qaLen - 1;
    }
    else {
        line->pos = line->rapos + line->raLen + line->eclip - 1;
        line->SQO = line->eclip;
        line->EQO = line->eclip + line->qaLen - 1;
    }
    // Need to pad the pos in case it is negative
    line->pos = padPos(line->pos);
    // Let's not calculate these again for this line.
    line->CIGARprocessed = true;
}

///////////////////////////////////////////////////////////////////////////////
// Process SAM Blocks
///////////////////////////////////////////////////////////////////////////////

// This is apparently no longer called.
void outputSAMBlock(splitLine_t * block, FILE * output) {
    for (splitLine_t * line = block; line != NULL; line = line->next) {
        writeLine(line, output);
    }
    disposeSplitLines(block);
}

void brokenBlock(splitLine_t *block, int count) {
    char * temp;
    asprintf(&temp, "samblaster: Can't find first and/or second of pair in sam block of length %d for id: %s\n%s%s:%s\n%s",
             count, block->fields[QNAME], "samblaster:    At location: ", block->fields[RNAME], block->fields[POS],
             "samblaster:    Are you sure the input is sorted by read ids?");
    fatalError(temp);
}

// Some fields for statistics.
UINT64 idCount = 0;
UINT64 dupCount = 0;
UINT64 discCount = 0;
UINT64 splitCount = 0;
UINT64 unmapClipCount = 0;
UINT64 unmatedCount = 0;

// This is the main workhorse that determines if lines are dups or not.
void markDupsDiscordants(splitLine_t * block, state_t * state) {
    splitLine_t * first = NULL;
    splitLine_t * second = NULL;
    int count = 0;
    for (splitLine_t * line = block; line != NULL; line = line->next) {
        count += 1;
        // Do this conversion once and store the result.
        line->flag = str2int(line->fields[FLAG]);
        // We make our duplicate decisions based solely on primary alignments.
        if (!isPrimaryAlignment(line)) continue;
        // Allow unpaired reads to go through (as the second so that signature is correct).
        // According to the SAM spec, this must be checked first.
        if (!isPaired(line)) second = line;
        // Figure out if this is the first half or second half of a pair.
        else if (isFirstRead(line))  first  = line;
        else if (isSecondRead(line)) second = line;
    }

    // Figure out what type of "pair" we have.
    bool orphan = false;
    bool dummyFirst = false;
    // First get rid of the useless case of having no first AND no second.
    if (first == NULL && second == NULL) goto outOfHere;
    // Now see if we have orphan with the unmapped read missing.
    if (first == NULL || second == NULL) {
        // Get the NULL one in the first slot.
        if (second == NULL) swapPtrs(&first, &second);
        // If the only read says its paired, and it is unmapped or its mate is mapped, something is wrong.
        if (isPaired(second) && (isUnmapped(second) || isNextMapped(second))) goto outOfHere;
        // If the only read we have is unmapped, then it can't be a dup.
        if (isUnmapped(second)) return;
        // Now MAKE a dummy record for the first read, but don't put it into the block.
        // That way we won't have to worry about it getting output, or when processing splitters etc.
        // But, we will have to remember to dispose of it.
        first = getSplitLine();
        // Set the flag field to what it would have been.
        // What if this "pair" is a singleton read?
        first->flag = isFirstRead(second) ? 0x85 : 0x45;
        orphan = true;
        dummyFirst = true;
    }
    else {
        // Handle the addition of MC and MQ tags if requested.
        if (state->addMateTags) {
            int mask = (FIRST_SEG | SECOND_SEG);
            // Process the first of the pair.
            // Get the list of reads that match the second of the pair.
            if (isMapped(first)) {
                int count = fillSplitterArray<false>(block, state, second->flag & mask, true);
                for (int i=0; i<count; ++i) {
                    splitLine_t * line = state->splitterArray[i];
                    addTag(line, "	MC:Z:", first->fields[CIGAR]);
                    addTag(line, "	MQ:i:", first->fields[MAPQ]);
                }
            }
            // Process the second of the pair.
            // Get the list of reads that match the first of the pair.
            if (isMapped(second)) {
                count = fillSplitterArray<false>(block, state, first->flag & mask, true);
                for (int i=0; i<count; ++i) {
                    splitLine_t * line = state->splitterArray[i];
                    addTag(line, "	MC:Z:", second->fields[CIGAR]);
                    addTag(line, "	MQ:i:", second->fields[MAPQ]);
                }
            }
        }

        // Never mark pairs as dups if both sides are unmapped.
        if (isUnmapped(first) && isUnmapped(second)) return;

        // We need to properly handle orphans to get the correct reference offsets and sequence numbers.
        orphan = (isUnmapped(first) || isUnmapped(second));
        // Orphan that needs to be swapped.
        // We need the unmapped one in the first slot so that they won't all collide in the hash table.
        if (isMapped(first) && isUnmapped(second)) {
            swapPtrs(&first, &second);
        }
    }

    // Now look for duplicates.
    if (!state->acceptDups) {
        // Calculate and store the second position and sequence name.
        calcOffsets(second);
        second->seqNum = getSeqNum(second, RNAME, state);
        UINT64 seqOff = state->seqOffs[second->seqNum]; //genome relative position
        second->binNum = (seqOff + second->pos) >> BIN_SHIFT;
        second->binPos = (seqOff + second->pos) &  BIN_MASK;

        if (orphan) {
            // We have an orphan, so we just zero out the pos and seqnum
            first->pos = 0;
            first->seqNum = 0;
            first->binNum = 0;
            first->binPos = 0;
        }
        else {
            // Not an orphan, so handle first on its own.
            calcOffsets(first);
            first->seqNum = getSeqNum(first, RNAME, state);
            seqOff = state->seqOffs[first->seqNum]; //genome relative position
            first->binNum = (seqOff + first->pos) >> BIN_SHIFT;
            first->binPos = (seqOff + first->pos) &  BIN_MASK;
        }

        // The fact of which alignment is first or second in the template is not relevant for determining dups.
        // Therefore, we normalize the pairs based on their characteristics.
        // We have already swapped orphans.
        // Otherwise, sort by pos, and if equal, sort by sequence num, then by strand.
        if (!orphan && needSwap(first, second)) swapPtrs(&first, &second);

        // Now find the signature of the pair.
        sgn_t sig = calcSig(first, second);
        // Calculate the offset into the signatures array.
        UINT32 off = calcSigArrOff(first, second, state->binCount);
        // Attempt insert into the sigs structure.
        // The return value will tell us if it was already there.
        bool insert = state->sigs[off].hashTableInsertLocked(sig);
        // Check if the insertion actually happened.
        if (!insert) {
            dupCount += 1;
            // We always mark all or none of a block as dup.
            for (splitLine_t * line = block; line != NULL; line = line->next) {
                markDup(line);
            }
        }
    }

    // If we have a dummy first, we can't have a discordant pair.
    if (dummyFirst) {
        disposeSplitLines(first);
        return;
    }

    // The first and second help us mark the discordants.
    // Both sides mapped, but pair not properly aligned.
    if (!orphan && isDiscordant(first)) {
        first->discordant = true;
        second->discordant = true;
    }
    return;

outOfHere:
    if (state->ignoreUnmated) { unmatedCount += 1; return; }
    else                       brokenBlock(block, count);
}

// Sort ascending in SQO.
int compQOs(const void * p1, const void * p2) {
    splitLine_t * l1 = (*(splitLine_t **)p1);
    splitLine_t * l2 = (*(splitLine_t **)p2);
    return (l1->SQO - l2->SQO);
}

template <bool excludeSecondaries>
int fillSplitterArray(splitLine_t * block, state_t * state, int mask, bool flagValue) {
    // Count the secondaries we have for this read (if any), and store their ptrs into an array.
    int count = 0;
    for (splitLine_t * line = block; line != NULL; line = line->next) {
        // For all the ones that are the current read of interest....
        // Check if they are a primary or complementary alignment.
        if (checkFlag(line, mask) == flagValue && !(excludeSecondaries && isSecondaryAlignment(line))) {
            // Add the ptr to this line to the sort array.
            // If it won't fit, double the array size.
            if (count >= state->splitterArrayMaxSize) {
                state->splitterArrayMaxSize *= 2;
                state->splitterArray = (splitLine_t **)(realloc(state->splitterArray,
                                                                state->splitterArrayMaxSize * sizeof(splitLine_t *)));
            }
            state->splitterArray[count] = line;
            count += 1;
        }
    }
    return count;
}

void markSplitterUnmappedClipped(splitLine_t * block, state_t * state, int mask, bool flagValue) {
    // Count the secondaries we have for this read (if any), and store their ptrs into an array.
    int count = fillSplitterArray<true>(block, state, mask, flagValue);

    // We have the lines of interest in an array.
    // Decide what to do next based on the number of reads.
    if (count == 0) return;
    if (count == 1) {
        if (state->unmappedClippedFile == NULL) return;
        // Process unmapped or clipped.
        splitLine_t * line = state->splitterArray[0];
        // Unmapped or clipped alignments should be primary.
        if (!isPrimaryAlignment(line)) return;
        if (isUnmapped(line)) {
            line->unmappedClipped = true;
        }
        else {
            // Process the CIGAR string.
            // As this is expensive, we delay as long as possible.
            calcOffsets(line);
            if (line->sclip >= state->minClip || line->eclip >= state->minClip) {
                line->unmappedClipped = true;
            }
        }
        return;
    }

    // See if we need to process for splitters.
    if (state->splitterFile == NULL || count > state->maxSplitCount) return;

    // Calculate the query positions (for sorting) and do other preprocessing.
    for (int i=0; i<count; i++) {
        splitLine_t * line = state->splitterArray[i];
        // Make sure the primary is mapped!
        if (isPrimaryAlignment(line) && isUnmapped(line)) return;
        calcOffsets(line);
    }

    // We need to sort it by strand normalized query offsets.
    qsort(state->splitterArray, count, sizeof(splitLine_t *), compQOs);

    // Now check for pairs that match the desired parameters.
    splitLine_t * left = state->splitterArray[0];
    splitLine_t * right;
    for (int i=1; i<count; i++) {
        right = state->splitterArray[i];

        // First check for minNonOverlap.
        // We don't allow negative overlap, as that will lead to wrong non-overlap calculation.
        int overlap = std::max(1 + std::min(left->EQO, right->EQO) - std::max(left->SQO, right->SQO), 0);
        int alen1 = 1 + left->EQO - left->SQO;
        int alen2 = 1 + right->EQO - right->SQO;
        int mno = std::min(alen1-overlap, alen2-overlap);
        if (mno < state->minNonOverlap) goto nextIter;

        // Now check for the deserts and diagonal difference.
        // If they are on different chroms or strands, they pass without the other checks.
        // Since we only care if the sequences are the same, we don't need seqNums.
        // Instead just compare the strings!
        if (streq(left->fields[RNAME], right->fields[RNAME]) && (isReverseStrand(left) == isReverseStrand(right))) {
            // The start and end diags might be different if there is an indel in the alignments.
            // So, we use the end on the left, and the start on the right.
            // This will give us the net diag difference between them.
            int leftDiag, rightDiag, insSize;
            if (isReverseStrand(left)) {
                leftDiag = getStartDiag(left);
                rightDiag = getEndDiag(right);
                insSize = rightDiag - leftDiag;
            }
            else {
                leftDiag = getEndDiag(left);
                rightDiag = getStartDiag(right);
                insSize = leftDiag - rightDiag;
            }
            int desert = right->SQO - left->EQO - 1;
            // The absolute value will handle both inserts and deletes.
            // So check indel is big enough, and that there are not too many unmapped bases.
            // Subtract the inadvertant desert gap of inserts.
            if ((abs(insSize) < state->minIndelSize)
                || ((desert > 0) && ((desert - (int)std::max(0, insSize)) > state->maxUnmappedBases)))
                goto nextIter;
        }
        // We made it through the gamet.
        // Mark this pair for output.
        left->splitter = true;
        right->splitter = true;

    nextIter:
        // Get ready for the next iteration.
        left = right;
    }
}

void writeUnmappedClipped(splitLine_t * line, state_t * state) {
    // Check if we are outputting fasta or fastq.
    if (state->unmappedFastq == -1)
        state->unmappedFastq = (streq(line->fields[QUAL], "*") ? 0 : 1);

    // Print the first line.
    char firstChar = (state->unmappedFastq) ? '@' : '>';
    if (isPaired(line)) fprintf(state->unmappedClippedFile, "%c%s_%d\n", firstChar, line->fields[QNAME], (isFirstRead(line) ? 1 : 2));
    else                fprintf(state->unmappedClippedFile, "%c%s\n", firstChar, line->fields[QNAME]);

    if (state->unmappedFastq) fprintf(state->unmappedClippedFile, "%s\n+\n%s\n", line->fields[SEQ], line->fields[QUAL]);
    else                       fprintf(state->unmappedClippedFile, "%s\n", line->fields[SEQ]);
}

void processSAMBlock(splitLine_t * block, state_t * state) {
    idCount += 1;
    // First mark dups and find the discordants.
    // These share a lot of looking at flag bits, so make sense to put together.
    markDupsDiscordants(block, state);
    // Look for splitters.
    // Since this is expensive, don't do it unless the user asked for them.
    if (state->splitterFile != NULL || state->unmappedClippedFile != NULL) {
        // Check the first read for splitter.
        markSplitterUnmappedClipped(block, state, FIRST_SEG, true);
        // Check the second read for splitter.
        markSplitterUnmappedClipped(block, state, SECOND_SEG, true);
        // Check for a singleton read
        markSplitterUnmappedClipped(block, state, MULTI_SEGS, false);
    }

    // Now do the output.
    for (splitLine_t * line = block; line != NULL; line = line->next) {
        // Do the unmapped file first, as it is not sam, and doesn't sew the line back together.
        if (state->unmappedClippedFile != NULL && line->unmappedClipped && !(state->excludeDups && isDuplicate(line))) {
            writeUnmappedClipped(line, state);
            unmapClipCount += 1;
        }

        // Write to the output file.
        if (!(state->removeDups && isDuplicate(line))) {
            writeLine(line, state->outputFile);
        }

        // Write to discordant file if appropriate.
        if (state->discordantFile != NULL && line->discordant && !(state->excludeDups && isDuplicate(line))) {
            writeLine(line, state->discordantFile);
            discCount += 1;
        }
        // Write to splitter file if appropriate.
        if (state->splitterFile != NULL && line->splitter && !(state->excludeDups && isDuplicate(line))) {
            writeSAMlineWithIdNum(line, state->splitterFile);
            splitCount += 1;
        }
    }

    disposeSplitLines(block);
}


