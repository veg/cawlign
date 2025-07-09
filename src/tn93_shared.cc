
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include "tn93_shared.h"

#define  RESOLVE_A 0x01
#define  RESOLVE_C 0x02
#define  RESOLVE_G 0x04
#define  RESOLVE_T 0x08

#define  TN93_MAX_DIST 1000.0

using namespace std;

const  char ValidChars[]       = "ACGTURYSWKMBDHVN?",
            ValidCharsAA[]     = "ACDEFGHIKLMNPQRSTVWYBZX?";
            
unsigned char   * resolveTheseAmbigs = (unsigned char   *)calloc (256,sizeof (unsigned char));
char   * reverseComplementCharacters = (char   *)calloc (256,sizeof (char));

double          resolve_fraction = 1.;
                  
static char empty_string        [] = "";

const long   resolutions [][4] = { {1,0,0,0},
  {0,1,0,0},
  {0,0,1,0},
  {0,0,0,1},
  {0,0,0,1}, // U - 4
  {1,0,1,0}, //RESOLVE_A | RESOLVE_G, // R - 5
  {0,1,0,1}, //RESOLVE_C | RESOLVE_T, // Y - 6
  {0,1,1,0}, //RESOLVE_C | RESOLVE_G, // S - 7
  {1,0,0,1}, //RESOLVE_A | RESOLVE_T, // W - 8
  {0,0,1,1}, //RESOLVE_G | RESOLVE_T, // K - 9
  {1,1,0,0}, //RESOLVE_A | RESOLVE_C, // M - 10
  {0,1,1,1}, // RESOLVE_C | RESOLVE_G | RESOLVE_T, // B - 11
  {1,0,1,1}, //RESOLVE_A | RESOLVE_G | RESOLVE_T, // D - 12
  {1,1,0,1}, //RESOLVE_A | RESOLVE_C | RESOLVE_T, // H - 13
  {1,1,1,0}, // RESOLVE_A | RESOLVE_C | RESOLVE_G, // V - 14
  {1,1,1,1}, // RESOLVE_A | RESOLVE_C | RESOLVE_G | RESOLVE_T , // N - 15
  {1,1,1,1} //RESOLVE_A | RESOLVE_C | RESOLVE_G | RESOLVE_T , // ? - 16
};






/*A.................Ala.................Alanine
B.................Asx.................Aspartic acid or Asparagine
C.................Cys.................Cysteine
D.................Asp.................Aspartic Acid
E.................Glu.................Glutamic Acid
F.................Phe.................Phenylalanine
G.................Gly.................Glycine
H.................His.................Histidine
I.................Ile.................Isoleucine
K.................Lys.................Lysine
L.................Leu.................Leucine
M.................Met.................Methionine
N.................Asn.................Asparagine
P.................Pro.................Proline
Q.................Gln.................Glutamine
R.................Arg.................Arginine
S.................Ser.................Serine
T.................Thr.................Threonine
V.................Val.................Valine
W.................Trp.................Tryptophan
X.................Xaa.................Any amino acid
Y.................Tyr.................Tyrosine
Z.................Glx.................Glutamine or Glutamic acid*/
 
 
// ACDEFGHIKLMNPQRSTVWYBZ?- 

const long   resolutions_AA [][20] = { {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                        {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                        {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                        {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                        {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                        {0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                        {0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                        {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0},
                                        {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0},
                                        {0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0},
                                        {0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0},
                                        {0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0},
                                        {0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
                                        {0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0},
                                        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0},
                                        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0},
                                        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0},
                                        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0},
                                        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0},
                                        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
//                                       A C D E F G H I L K M N P Q R S T V W Y
                                        {0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0}, // B
                                        {0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0}, // Z
                                        {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}, // X
                                        {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}

};

#define N_CHAR 14
#define GAP    255
#define GAP_AA 255

const  double   resolutionsCount[] = { 1.f,
  1.f,
  1.f,
  1.f,
  1.f,
  1./2.f, // R
  1./2.f, // Y
  1./2.f, // S
  1./2.f, // S
  1./2.f, // W
  1./2.f, // K
  1./2.f, // M
  1./3.f, // B
  1./3.f, // D
  1./3.f, // H
  1./3.f, // V
  1./4.f, // N
  1./4.f // ?
};

const  double   resolutionsCount_AA[] = {
    1.f,
    1.f,
    1.f,
    1.f,
    1.f,
    1.f,
    1.f,
    1.f,
    1.f,
    1.f,
    1.f,
    1.f,
    1.f,
    1.f,
    1.f,
    1.f,
    1.f,
    1.f,
    1.f,
    1.f,
    0.5f,
    0.5f,
    0.05f,
    0.05f
};


char validFlags [256];

  //---------------------------------------------------------------

/**
 * Initializes the valid character flags for nucleotide or amino acid sequences.
 *
 * This function sets up the `validFlags` array, which maps valid characters for 
 * either nucleotide or amino acid sequences. The behavior changes depending on the 
 * `doAminoAcid` flag. If `id_map` is true, the flags map to the character's ASCII value; 
 * otherwise, they map to an index.
 *
 * If a subset of nucleotide resolution is provided, it marks ambiguous bases to be resolved.
 *
 * @param doAminoAcid A boolean flag indicating whether to initialize the flags for amino acids (true) or nucleotides (false).
 * @param resolutionSubset A string representing a subset of nucleotide resolution characters, used to handle ambiguous characters. Can be nullptr for amino acids or no subset.
 * @param id_map A boolean flag indicating whether to map characters to their ASCII values (true) or their indices in the valid character array (false).
 */
void initAlphabets (bool doAminoAcid, char * resolutionSubset, bool id_map) {
    
  for (int i = 0; i < 256; i++) {
    validFlags [i] = -1;
    reverseComplementCharacters [i] = -1;
  }
  
  if (doAminoAcid) {
    if (id_map) {
        for (unsigned int i = 0; i < strlen (ValidCharsAA); i++)
          validFlags [(unsigned char)ValidCharsAA[i]] = (unsigned char)ValidCharsAA[i];
    } else {
        for (unsigned int i = 0; i < strlen (ValidCharsAA); i++)
            validFlags [(unsigned char)ValidCharsAA[i]] = i;

    }
    
    
  } else {  
     // set up reverse complements
      
      
      const long   resolutions [][4] = { {1,0,0,0},
        {0,1,0,0},
        {0,0,1,0},
        {0,0,0,1},
        {0,0,0,1}, // U - 4
        {1,0,1,0}, //RESOLVE_A | RESOLVE_G, // R - 5
        {0,1,0,1}, //RESOLVE_C | RESOLVE_T, // Y - 6
        {0,1,1,0}, //RESOLVE_C | RESOLVE_G, // S - 7
        {1,0,0,1}, //RESOLVE_A | RESOLVE_T, // W - 8
        {0,0,1,1}, //RESOLVE_G | RESOLVE_T, // K - 9
        {1,1,0,0}, //RESOLVE_A | RESOLVE_C, // M - 10
        {0,1,1,1}, // RESOLVE_C | RESOLVE_G | RESOLVE_T, // B - 11
        {1,0,1,1}, //RESOLVE_A | RESOLVE_G | RESOLVE_T, // D - 12
        {1,1,0,1}, //RESOLVE_A | RESOLVE_C | RESOLVE_T, // H - 13
        {1,1,1,0}, // RESOLVE_A | RESOLVE_C | RESOLVE_G, // V - 14
        {1,1,1,1}, // RESOLVE_A | RESOLVE_C | RESOLVE_G | RESOLVE_T , // N - 15
        {1,1,1,1} //RESOLVE_A | RESOLVE_C | RESOLVE_G | RESOLVE_T , // ? - 16
      };
      
     reverseComplementCharacters['A'] = 'T';
     reverseComplementCharacters['C'] = 'G';
     reverseComplementCharacters['G'] = 'C';
     reverseComplementCharacters['T'] = 'A';
     reverseComplementCharacters['N'] = 'N';
     reverseComplementCharacters['?'] = '?';
     reverseComplementCharacters['R'] = 'Y';
     reverseComplementCharacters['Y'] = 'R';
     reverseComplementCharacters['S'] = 'W';
     reverseComplementCharacters['W'] = 'S';
     reverseComplementCharacters['B'] = 'V';
     reverseComplementCharacters['V'] = 'B';
     reverseComplementCharacters['D'] = 'H';
     reverseComplementCharacters['H'] = 'D';

     if (id_map) {
         for (unsigned int i = 0; i < strlen (ValidChars); i++)
          validFlags [(unsigned char)ValidChars[i]] = (unsigned char)ValidChars[i];
     } else {
         for (unsigned int i = 0; i < strlen (ValidChars); i++)
             validFlags [(unsigned char)ValidChars[i]] = i;
     }
      
    if (resolutionSubset) {
      unsigned long subset_length = strlen (resolutionSubset);
      for (unsigned long rc = 0; rc < subset_length; rc++) {
        unsigned char rcc = toupper( (resolutionSubset[rc]));
        if (validFlags[rcc] > 3) {
          resolveTheseAmbigs[(unsigned char)validFlags[rcc]] = 1;
        }
      }
    }
  }
  
}

  //---------------------------------------------------------------
/**
 * Returns the character mapped to the given index in the valid character set.
 *
 */
const char unmap_char (unsigned char c, bool do_aa) {
  return do_aa ? ValidCharsAA[c] : ValidChars[c];
}

  //---------------------------------------------------------------

#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* prototypes */


/**
 * Initializes mt[N] (Mersenne Twister pseudorandom number generator) with a given seed.
 *
 */
void init_genrand(unsigned long s)
{
  mt[0]= s & 0xffffffffUL;
  for (mti=1; mti<N; mti++) {
    mt[mti] =
    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
    /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
    /* In the previous versions, MSBs of the seed affect   */
    /* only MSBs of the array mt[].                        */
    /* 2002/01/09 modified by Makoto Matsumoto             */
    mt[mti] &= 0xffffffffUL;
    /* for >32 bit machines */
  }
}


/**
 * Generates a random 32-bit unsigned integer within the range [0, 0xffffffff].
 */
unsigned long genrand_int32(void)
{
  unsigned long y;
  static unsigned long mag01[2]= {0x0UL, MATRIX_A};
  /* mag01[x] = x * MATRIX_A  for x=0,1 */
  
  if (mti >= N) { /* generate N words at one time */
    int kk;
    
    if (mti == N+1)   /* if init_genrand() has not been called, */
      init_genrand(5489UL); /* a default initial seed is used */
    
    // Generate numbers using bitwise operations
    for (kk=0; kk<N-M; kk++) {
      y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
      mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    for (; kk<N-1; kk++) {
      y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
      mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
    mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];
    
    mti = 0; /* Reset the index after generating the numbers */
  }
  
  y = mt[mti++]; /* Get the next number from the state array */
  
  /* Tempering */
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);
  
  return y;
}


  //---------------------------------------------------------------

/**
 * Returns the count of resolutions for a given character.
 */
const double resolution_count (unsigned char c, bool do_aa) {
  
  if (do_aa) {
    return resolutionsCount_AA[c];
  } 
   return resolutionsCount [c];
  
}

  //---------------------------------------------------------------

/**
 * Calculates the length of a string within a vector of lengths.
 */
long stringLength (Vector& lengths, unsigned long index)
{
  if (index < lengths.length() - 1)
    return lengths.value(index+1) - lengths.value(index) - 1;
  
  return -1;
}

  //---------------------------------------------------------------

/**
 * Retrieves the text of a string from a buffer based on its index. 
 * The buffer contains concatenated strings, and their lengths are stored in a 
 * separate vector.
 */
char* stringText (const StringBuffer& strings, const Vector& lengths, unsigned long index)
{
  if (index < lengths.length() - 1L) 
    return strings.getString() + lengths.value(index);  
  return empty_string;
}




  //---------------------------------------------------------------

/**
 * Adds a sequence to a list of sequences, ensuring consistent length.
 *
 * @param sequences A buffer containing concatenated sequences.
 * @param seqLengths A vector storing the lengths of each sequence.
 * @param firstSequenceLength A reference to store the length of the first sequence.
 * @param names A buffer containing concatenated sequence names.
 * @param nameLengths A vector storing the lengths of each sequence name.
 */
void addASequenceToList (StringBuffer& sequences, Vector& seqLengths, long &firstSequenceLength, StringBuffer& names, Vector& nameLengths)
{
  sequences.appendChar ('\0');
  seqLengths.appendValue (sequences.length());
  if (seqLengths.length() == 2)
      {
    firstSequenceLength = stringLength (seqLengths, 0);
    
    if (firstSequenceLength <= 0)
        {
      cerr << "First sequence length must be positive." << endl;
      exit (1);
        }
      }
  else
      {
    if (stringLength (seqLengths, seqLengths.length()-2) != firstSequenceLength)
        {
      cerr << "All sequences must have the same length (" << firstSequenceLength << "), but sequence '" << stringText (names, nameLengths, nameLengths.length()-2) << "' had length " << stringLength (seqLengths, seqLengths.length()-2);
      exit (1);
        }
      }
}


//---------------------------------------------------------------

/**
 * Reads sequences from a FASTA file and processes them based on the specified options.
 *
 * This function reads sequences from a given FASTA file stream (`F`) and processes them 
 * into buffers for names, sequences, and lengths. It supports both batch and one-by-one 
 * sequence processing, and allows filtering based on an inclusion probability. Sequence 
 * instances (if provided) track how many times each sequence should be included based 
 * on a delimiter character. The function also supports progress reporting for large files.
 *
 * @param F The file stream pointer from which the FASTA sequences are read.
 * @param automatonState The current state of the finite automaton controlling the reading process (0, 1, or 2).
 * @param names A buffer to store sequence names.
 * @param sequences A buffer to store the actual sequences.
 * @param nameLengths A vector to store the lengths of the sequence names.
 * @param seqLengths A vector to store the lengths of the sequences.
 * @param firstSequenceLength A reference to store the length of the first sequence for length consistency checking.
 * @param oneByOne If true, the function processes sequences one at a time, resetting buffers after each sequence.
 * @param sequenceInstances A vector to store the number of instances for each sequence (if applicable). Can be `nullptr`.
 * @param sep The separator character used to detect sequence instances in sequence names.
 * @param include_prob A probability between 0 and 1 to determine whether sequences should be included. If less than 1, sequences are randomly filtered.
 * @param progress If true, progress information is printed during processing.
 * @return An integer indicating the status of the read operation:
 *         - 0: End of file or successful read of all sequences.
 *         - 1: An error occurred during processing.
 *         - 2: Sequence successfully read in one-by-one mode.
 *         - 3: Sequence successfully read and processed at the end of the file.
 */
int readFASTA (FILE* F, char& automatonState,  StringBuffer &names, 
               StringBuffer& sequences, Vector &nameLengths, Vector &seqLengths, 
               long& firstSequenceLength, bool oneByOne, 
               Vector* sequenceInstances, char sep,
               double include_prob, bool progress) {
  
  unsigned long up_to = 0L;
  
  if (oneByOne) {
    sequences.resetString();
    names.resetString();
    include_prob = 1.;
  }
  
  if (include_prob < 1.) {
    up_to = RAND_RANGE * include_prob;
  }
    
  time_t before, after;
  if (progress) {
    time(&before);
  }

  
  bool include_me = true;
  long read_counter = 0;
    
  flockfile(F);
    
  try {
  
      while (1) {
        int currentC = getc_unlocked (F);
        //cout << "State: " << int(automatonState) << "/'" << char(currentC) << "'" << endl;
        if (feof_unlocked (F))
          break;
          
        switch (automatonState) {
          case 0: {
            if (currentC == '>' || currentC == '#') {
              automatonState = 1;
              if (sequenceInstances == NULL && include_prob < 1.) {
                include_me = genrand_int32() < up_to;
              }
            }
            break;
          }
          case 1: {
            if (currentC == '\n' || currentC == '\r') {
              if (include_me) {
                  names.appendChar   ('\0');
     
                long this_name_l;
                if (oneByOne) {
                  this_name_l = names.length ()-1;
                  
                } else {
                  nameLengths.appendValue (names.length());
                  this_name_l = stringLength (nameLengths, nameLengths.length()-2);
                }

                
                if (this_name_l <= 0) {
                  throw std::string("Sequence names must be non-empty.");
                }
                automatonState = 2;

                if (sequenceInstances) {
                    unsigned long count = 1L;
                   if (this_name_l >= 3) {
                      long sep_loc = 0, ll = names.length();
                      for (sep_loc = 2; sep_loc < this_name_l; sep_loc ++) {
                        if (names.getChar(ll-sep_loc-1) == sep) {
                          break;
                        }
                      }
                      if (sep_loc < this_name_l) {
                        count = atoi (names.getString() + (ll-sep_loc));
                      }
                      if (count < 1) {
                        count = 1;
                      }
                      
                      unsigned long resampled_prob = 0UL;
                      
                      if (include_prob < 1.) {
                        for (long k = 0; k < count; k ++) {
                          resampled_prob += genrand_int32() < up_to;
                        }
                      } else {
                        resampled_prob = count;
                      }
                     
                     //cerr << count << " -> " << resampled_prob << endl;
                     
                      if (resampled_prob == 0UL) {
                        if (oneByOne) {
                          names.resetString();
                        } else {
                          names.reset_length (nameLengths.value (nameLengths.length()-2));
                        }
                        nameLengths.remove(nameLengths.length()-1);
                        include_me = false;
                        continue;
                      } else {
                        count = resampled_prob;
                      }
                      
                    }
                    if (oneByOne) {
                      sequenceInstances->resetVector();
                    }
                    sequenceInstances->appendValue(count);
                  }
               }
              
                
            }
            else {
              if (include_me) {
                names.appendChar(currentC);
              }
            }
            break;
          }
          case 2: {
            currentC = toupper (currentC);
            if (validFlags [currentC] >= 0) {
              if (include_me)
                //cout << "Append " << (char)currentC << ":" << currentC << endl;
                sequences.appendChar (currentC);
            }
            else {
              if (currentC == '>' || currentC == '#') {
                automatonState = 1;
                if (include_me) {
                  if (oneByOne) {
                    if (firstSequenceLength == 0) {
                      firstSequenceLength = sequences.length()-1;
                    }
                    //cerr << endl << "Returning a sequence" << firstSequenceLength << endl;
                    automatonState = 0;
                    sequences.appendChar ('\0');
                    ungetc (currentC, F);
                    funlockfile (F);
                    return 2;
                  }
                  addASequenceToList (sequences, seqLengths, firstSequenceLength, names, nameLengths);
                  read_counter++;
                  if (progress && read_counter % 1024 == 0) {
                      time(&after);
                      cerr << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                              "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                              "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bProgress"
                              ":"
                           << setw(8) << read_counter << " sequences read (" << setw(12) << std::setprecision(3)
                           << read_counter / difftime(after, before) << " seqs/sec)";

                      after = before;
                  }
                }
                if (sequenceInstances == NULL && include_prob < 1.) {
                  include_me = genrand_int32() < up_to;
                } else {
                  include_me = true;
                }
                  
              }
            }
            break;
          }
        }
      }
      
      if (automatonState == 2 || (oneByOne && automatonState == 0)) {
          if (include_me) {
              if (oneByOne) {
                  if (firstSequenceLength == 0) {
                      firstSequenceLength = sequences.length()-1;
                  }
                  automatonState = 0;
                  sequences.appendChar ('\0');
                  funlockfile (F);
                  return 3;
              } else {
                  addASequenceToList (sequences, seqLengths, firstSequenceLength, names, nameLengths);
              }
          }
          automatonState = 1;
      } else {
        char err[256] ;
        snprintf (err, 255, "Unexpected end of file: state %d", automatonState);
        throw std::string(err);
      }
  } catch (std::string const err) {
      cerr << err << endl;
      funlockfile (F);
      return 1;
  }
  return 0;
}

//---------------------------------------------------------------

/**
 * Reverse complement a nucleotide sequence (jn place).
 *
 * This function reads a capitalized IUPAC compliant string and performs reverse complementation
 * It will generate return 1 if invalid characters are present in the sequence to RC; the original string data will be corrupted via partial RC
 *
 * @param sequences The string buffer object to use
 * @param from start index of the string in the buffer
 * @param to end index of the string in the buffer (inclusive!)
 * @return 0 if success, otherwise non-zero
 *
 */
int reverseComplement (StringBuffer& sequence, unsigned long from, unsigned long to) {
    
    for (unsigned long i = from, j = to; i <= j; i++, j--) {
        char leading_char  = reverseComplementCharacters[sequence.getChar(i)],
             trailing_char = reverseComplementCharacters[sequence.getChar(j)];
    
        if (leading_char > 0 && trailing_char > 0) {
            sequence.setChar (i, trailing_char);
            sequence.setChar (j, leading_char);
        } else {
            return 1;
        }
    }
    
    
    return 0;
}




