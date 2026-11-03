
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include "tn93_shared.h"

using namespace std;

const  char ValidChars[]       = "ACGTURYSWKMBDHVN?",
            ValidCharsAA[]     = "ACDEFGHIKLMNPQRSTVWYBZX?";
            
char   * reverseComplementCharacters = (char   *)calloc (256,sizeof (char));

int8_t validFlags [256];

  //---------------------------------------------------------------

/**
 * Initializes the valid character flags for nucleotide or amino acid sequences.
 *
 * This function sets up the `validFlags` array, which maps each valid character to its
 * index in the nucleotide or amino acid alphabet (and every other character to -1), and,
 * for nucleotides, the reverse complement table used by `reverseComplement`.
 *
 * @param doAminoAcid A boolean flag indicating whether to initialize the flags for amino acids (true) or nucleotides (false).
 */
void initAlphabets (bool doAminoAcid) {
    
  for (int i = 0; i < 256; i++) {
    validFlags [i] = -1;
    reverseComplementCharacters [i] = 0; // 0 = no complement
  }
  
  if (doAminoAcid) {
    for (unsigned int i = 0; i < strlen (ValidCharsAA); i++)
        validFlags [(unsigned char)ValidCharsAA[i]] = i;
  } else {  
     // set up reverse complements
     reverseComplementCharacters['A'] = 'T';
     reverseComplementCharacters['C'] = 'G';
     reverseComplementCharacters['G'] = 'C';
     reverseComplementCharacters['T'] = 'A';
     reverseComplementCharacters['U'] = 'A';
     reverseComplementCharacters['N'] = 'N';
     reverseComplementCharacters['?'] = '?';
     reverseComplementCharacters['R'] = 'Y';
     reverseComplementCharacters['Y'] = 'R';
     reverseComplementCharacters['S'] = 'S';
     reverseComplementCharacters['W'] = 'W';
     reverseComplementCharacters['K'] = 'M';
     reverseComplementCharacters['M'] = 'K';
     reverseComplementCharacters['B'] = 'V';
     reverseComplementCharacters['V'] = 'B';
     reverseComplementCharacters['D'] = 'H';
     reverseComplementCharacters['H'] = 'D';

     for (unsigned int i = 0; i < strlen (ValidChars); i++)
         validFlags [(unsigned char)ValidChars[i]] = i;
  }
  
}

//---------------------------------------------------------------

/**
 * Reads the next sequence from a FASTA file.
 *
 * Sequences are read one at a time: `names` and `sequences` are reset on entry and
 * receive the next record (both NUL-terminated). Characters that are not in the
 * current alphabet (see `initAlphabets`) are skipped.
 *
 * @param F The file stream pointer from which the FASTA sequences are read.
 * @param automatonState The current state of the finite automaton controlling the reading process (0, 1, or 2).
 * @param names A buffer to store the sequence name.
 * @param sequences A buffer to store the sequence.
 * @param firstSequenceLength Set to (sequence length - 1) if it is 0 on entry.
 * @return An integer indicating the status of the read operation:
 *         - 1: An error occurred during processing.
 *         - 2: Sequence successfully read; more may follow.
 *         - 3: Sequence successfully read at the end of the file.
 */
int readFASTA (FILE* F, char& automatonState,  StringBuffer &names, 
               StringBuffer& sequences, long& firstSequenceLength) {
  
  sequences.resetString();
  names.resetString();
    
  flockfile(F);
    
  try {
  
      while (1) {
        int currentC = getc_unlocked (F);
        if (feof_unlocked (F))
          break;
          
        switch (automatonState) {
          case 0: {
            if (currentC == '>' || currentC == '#') {
              automatonState = 1;
            }
            break;
          }
          case 1: {
            if (currentC == '\n' || currentC == '\r') {
              names.appendChar   ('\0');
              if (names.length () - 1 <= 0) {
                throw std::string("Sequence names must be non-empty.");
              }
              automatonState = 2;
            }
            else {
              names.appendChar(currentC);
            }
            break;
          }
          case 2: {
            currentC = toupper (currentC);
            if (validFlags [currentC] >= 0) {
              sequences.appendChar (currentC);
            }
            else {
              if (currentC == '>' || currentC == '#') {
                if (firstSequenceLength == 0) {
                  firstSequenceLength = sequences.length()-1;
                }
                automatonState = 0;
                sequences.appendChar ('\0');
                ungetc (currentC, F);
                funlockfile (F);
                return 2;
              }
            }
            break;
          }
        }
      }
      
      if (automatonState == 2 || automatonState == 0) {
          if (firstSequenceLength == 0) {
              firstSequenceLength = sequences.length()-1;
          }
          automatonState = 0;
          sequences.appendChar ('\0');
          funlockfile (F);
          return 3;
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
    
    if (sequence.length() == 0 || from > to || to >= (unsigned long)sequence.length()) return 0;
    unsigned long i = from, j = to;
    while (i < j) {
        char leading_char  = reverseComplementCharacters[(unsigned char)sequence.getChar(i)],
             trailing_char = reverseComplementCharacters[(unsigned char)sequence.getChar(j)];

        if (leading_char > 0 && trailing_char > 0) {
            sequence.setChar (i, trailing_char);
            sequence.setChar (j, leading_char);
        } else {
            return 1;
        }
        i++;
        j--;
    }
    if (i == j) {
        char c = reverseComplementCharacters[(unsigned char)sequence.getChar(i)];
        if (c > 0) {
            sequence.setChar (i, c);
        } else {
            return 1;
        }
    }
    
    
    return 0;
}
