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
#include "stringBuffer.h"

using namespace std;

#define  RESOLVE_A      0x01
#define  RESOLVE_C      0x02
#define  RESOLVE_G      0x04
#define  RESOLVE_T      0x08


#define  RESOLVE        0
#define  AVERAGE        1
#define  SKIP           2
#define  GAPMM          3
#define  SUBSET         4
#define  MISMATCH       5
#define  INFORMATIVE    6

#define RAND_RANGE 0xffffffffUL /* Maximum value returned by genrand_int32 */

#define MIN(a,b) (a) < (b) ? (a) : (b)
#define MAX(a,b) (a) > (b) ? (a) : (b)

/**
 * @brief A struct to hold the first and last non-gap characters in a sequence.
 * 
 */
struct sequence_gap_structure {
  
  long first_nongap,
       last_nongap,
       resolved_start,
       resolved_end;
  
  /**
   * @brief Construct a new sequence gap structure object
   * 
   */
  sequence_gap_structure (void) {
    first_nongap    = LONG_MAX;
    last_nongap     = 0L;
    resolved_start = 0L;
    resolved_end   = 0L;
  }
  
};

/**
 * @brief Initialize the random number generator.
 * 
 * @param s The seed for the random number generator.
 */
void init_genrand(unsigned long s);
/**
 * @brief Generate a random 32-bit integer.
 * 
 * @return unsigned long A random 32-bit integer.
 */
unsigned long genrand_int32(void);
/**
 * @brief Compute the TN93 distance between two sequences.
 * 
 * @param s1 The first sequence.
 * @param s2 The second sequence.
 * @param L The length of the sequences.
 * @param matchMode The mode to use for matching ambiguous characters.
 * @param randomize A pointer to an array of random numbers.
 * @param min_overlap The minimum overlap between the two sequences.
 * @param resolvedCounts A pointer to an array of resolved counts.
 * @param zero_diag A flag to indicate whether to zero the diagonal.
 * @param cnt The number of sequences.
 * @param count1 The number of sequences in the first group.
 * @param count2 The number of sequences in the second group.
 * @param gaps1 A pointer to a sequence gap structure for the first sequence.
 * @param gaps2 A pointer to a sequence gap structure for the second sequence.
 * @return double The TN93 distance between the two sequences.
 */
double		computeTN93 (const char * s1, const char *s2,  const unsigned long L, const char matchMode, const long * randomize, const long min_overlap, unsigned long* = NULL, const double = 0.0, const unsigned long cnt = 0, const long count1 = 1, const long count2 = 1, const sequence_gap_structure * = NULL, const sequence_gap_structure * = NULL);

/**
 * @brief Compute the differences between two sequences.
 * 
 * @param s1 The first sequence.
 * @param s2 The second sequence.
 * @param L The length of the sequences.
 * @param matchMode The mode to use for matching ambiguous characters.
 * @param storage A vector to store the differences in.
 * @param gaps1 A pointer to a sequence gap structure for the first sequence.
 * @param gaps2 A pointer to a sequence gap structure for the second sequence.
 * @return long The number of differences between the two sequences.
 */
long   computeDifferences (const char * s1,
                           const char *s2,
                           const unsigned long L,
                           const char matchMode,
                           Vector& storage,
                           const sequence_gap_structure * = NULL,
                           const sequence_gap_structure * = NULL);


/**
 * @brief Get the length of a string in a vector of strings.
 * 
 * @param lengths A vector of string lengths.
 * @param index The index of the string to get the length of.
 * @return long The length of the string.
 */
long stringLength (Vector& lengths, unsigned long index);
/**
 * @brief Get the text of a string in a vector of strings.
 * 
 * @param strings A string buffer containing the strings.
 * @param lengths A vector of string lengths.
 * @param index The index of the string to get the text of.
 * @return char* A pointer to the text of the string.
 */
char* stringText (const StringBuffer& strings, const Vector& lengths, unsigned long index);
/**
 * @brief Add a sequence to a list of sequences.
 * 
 * @param sequences A string buffer to add the sequence to.
 * @param seqLengths A vector of sequence lengths.
 * @param firstSequenceLength The length of the first sequence.
 * @param names A string buffer of sequence names.
 * @param nameLengths A vector of sequence name lengths.
 */
void addASequenceToList (StringBuffer& sequences, Vector& seqLengths, long &firstSequenceLength, StringBuffer& names, Vector& nameLengths);
/**
 * @brief Read a FASTA file.
 * 
 * @param F The file to read from.
 * @param automatonState The state of the automaton.
 * @param names A string buffer to store the names of the sequences in.
 * @param sequences A string buffer to store the sequences in.
 * @param nameLengths A vector to store the lengths of the names in.
 * @param seqLengths A vector to store the lengths of the sequences in.
 * @param firstSequenceLength The length of the first sequence.
 * @param oneByOne A flag to indicate whether to read the sequences one byone.
 * @param sequenceInstances A vector to store the sequence instances in.
 * @param sep The separator to use between the name and the sequence.
 * @param include_prob The probability of including a sequence.
 * @param show_progress A flag to indicate whether to show progress.
 * @return int The number of sequences read.
 */
int readFASTA (FILE* F, char& automatonState,  StringBuffer &names, StringBuffer& sequences, Vector &nameLengths, Vector &seqLengths, long& firstSequenceLength, bool oneByOne = false,  Vector* sequenceInstances = NULL, char sep = ':', double include_prob = 1.0, bool show_progress = false);
/**
 * @brief Dump a sequence to a FASTA file.
 * 
 * @param index The index of the sequence to dump.
 * @param output The file to dump the sequence to.
 * @param firstSequenceLength The length of the first sequence.
 * @param d A pointer to an array of doubles.
 * @param append A flag to indicate whether to append to the file.
 * @param from The starting position of the sequence to dump.
 * @param to The ending position of the sequence to dump.
 */
void dump_sequence_fasta (unsigned long index, FILE* output, long firstSequenceLength, double * d = NULL, bool = false, unsigned long from = 0L, unsigned long to = 0L);
/**
 * @brief Initialize the alphabets.
 * 
 * @param dna A flag to indicate whether to initialize the DNA alphabet.
 * @param map A pointer to a character map.
 * @param id_map A flag to indicate whether to use an identity map.
 */
void initAlphabets(bool = false, char * = NULL, bool id_map = false);
/**
 * @brief Merge two sequences.
 * 
 * @param source The source sequence.
 * @param target The target sequence.
 * @param sequence_length The length of the sequences.
 */
void merge_two_sequences (const char* source, char* target, const long sequence_length);
/**
 * @brief Check if two sequences are a perfect match.
 * 
 * @param source The source sequence.
 * @param target The target sequence.
 * @param sequence_length The length of the sequences.
 * @return long The number of mismatches.
 */
long perfect_match (const char* source, char* target, const long sequence_length);
/**
 * @brief Dump a sequence to a FASTA file.
 * 
 * @param source The sequence to dump.
 * @param sequence_length The length of the sequence.
 * @param output The file to dump the sequence to.
 * @param newln A flag to indicate whether to add a newline.
 * @param append A flag to indicate whether to append to the file.
 * @param from The starting position of the sequence to dump.
 * @param to The ending position of the sequence to dump.
 */
void dump_fasta (const char* source, const long sequence_length, FILE* output, bool newln = true, bool = false, unsigned long from = 0L, unsigned long to = 0L);

/**
 * @brief Reverse complement a sequence.
 * 
 * @param sequence The sequence to reverse complement.
 * @param from The starting position of the sequence to reverse complement.
 * @param to The ending position of the sequence to reverse complement.
 * @return int 0 on success, 1 on failure.
 */
int    reverseComplement (StringBuffer& sequence, unsigned long from, unsigned long to);
/**
 * @brief Describe the gaps in a sequence.
 * 
 * @param source The sequence to describe.
 * @param sequence_length The length of the sequence.
 * @param char_count The number of characters in the alphabet.
 * @return sequence_gap_structure A struct containing the first and last non-gap characters in the sequence.
 */
struct sequence_gap_structure describe_sequence (const char* source, const unsigned long sequence_length, const unsigned long char_count = 4UL);

/**
 * @brief Resolve an ambiguous character.
 * 
 * @param c The character to resolve.
 * @param dna A flag to indicate whether to resolve the character as DNA.
 * @param protein A flag to indicate whether to resolve the character as a protein.
 * @return const long* A pointer to an array of resolved characters.
 */
const long * resolve_char (unsigned char, bool = false, bool = true);
/**
 * @brief Get the number of resolutions for an ambiguous character.
 * 
 * @param c The character to get the number of resolutions for.
 * @param dna A flag to indicate whether to get the number of resolutions for the character as DNA.
 * @return const double The number of resolutions for the character.
 */
const double resolution_count (unsigned char, bool = false);
/**
 * @brief Unmap a character.
 * 
 * @param c The character to unmap.
 * @param dna A flag to indicate whether to unmap the character as DNA.
 * @return const char The unmapped character.
 */
const char unmap_char (unsigned char, bool = false);
/**
 * @brief Unpack a difference into a location and an alternative.
 * 
 * @param diff The difference to unpack.
 * @param location The location of the difference.
 * @param alt The alternative of the difference.
 */
inline void unpack_difference (long diff, long& location, unsigned& alt) {
    location = diff >> 8;
    alt = diff & 0xff;
}


extern StringBuffer names,
       sequences;

extern unsigned char * resolveTheseAmbigs;

extern double   resolve_fraction;

extern Vector       nameLengths,
       seqLengths,
       workingNodes,
       nodeParents;

extern VectorFP distanceEstimates;
extern const  double  resolutionsCount [];
extern char validFlags[];

#endif
