#ifndef 	__TN93SHARED__
#define 	__TN93SHARED__

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <climits>
#include <cstdint>
#include "stringBuffer.h"

using namespace std;

#define MIN(a,b) (a) < (b) ? (a) : (b)
#define MAX(a,b) (a) > (b) ? (a) : (b)

/**
 * @brief Read the next sequence from a FASTA file.
 *
 * @param F The file to read from.
 * @param automatonState The state of the FASTA reading automaton (0 - between sequences, 1 - reading a name, 2 - reading a sequence); carried across calls.
 * @param names A string buffer that receives the sequence name.
 * @param sequences A string buffer that receives the sequence.
 * @param firstSequenceLength Set to (sequence length - 1) if it is 0 on entry.
 * @return int 1 on a parse error, 2 if a sequence was read and more may follow, 3 if the last sequence was read, 0 if there was nothing left to read.
 */
int readFASTA (FILE* F, char& automatonState, StringBuffer &names, StringBuffer& sequences, long& firstSequenceLength);
/**
 * @brief Initialize the alphabets (validFlags and reverse complement tables).
 *
 * @param doAminoAcid Set up the amino-acid alphabet instead of the nucleotide one.
 */
void initAlphabets(bool doAminoAcid = false);
/**
 * @brief Reverse complement a nucleotide sequence in place.
 *
 * @param sequence The string buffer holding the sequence.
 * @param from The index of the first character to reverse complement.
 * @param to The index of the last character to reverse complement (inclusive).
 * @return int 0 on success, 1 if a character could not be complemented.
 */
int    reverseComplement (StringBuffer& sequence, unsigned long from, unsigned long to);

extern int8_t validFlags[];

#endif
