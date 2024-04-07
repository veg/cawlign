#ifndef SCORING_H
#define SCORING_H

#include <iostream>
#include "alignment.h"
#include "argparse.hpp"
#include "tn93_shared.h"

using namespace std;
using namespace argparse;


class CawalignSimpleScores {
    public:
        CawalignSimpleScores  (
                                const char * _alphabet,
                                const cawlign_fp * _scoring_matrix,
                                const cawlign_fp _open_gap_reference,
                                const cawlign_fp _open_gap_query,
                                const cawlign_fp _extend_gap_reference,
                                const cawlign_fp _extend_gap_query
                              );
    
        CawalignSimpleScores  (ConfigParser*);
        CawalignSimpleScores  (void) {D=0;};
        virtual ~CawalignSimpleScores (void) {};
    
        StringBuffer           alphabet;
        /*
            ordered characters that are included in the scoring matrix
        */
        unsigned int          D;
        // the number of characters in the string
    
        long                  char_map [255];
        /*
            for each ASCII character, this will map the character to the corresping entry the scoring matrix
            all characters NOT in `alphabet` get mapped to index D (the 'not defined' character)
        */
    
        VectorFP         scoring_matrix;
        /*
                A (D+1)x(D+1) scoring matrix where element (i,j) gives the score of matching (or mis-matching)
                the D-th row/column is for matchign a character NOT in the alphabet
                While generally symmetric, an asymmetric matrix can also be meaningful if there is some reason to have
                substitutions in reference/query weigted differently
        */
    
        cawlign_fp              open_gap_reference,
                            open_gap_query,
                            extend_gap_query,
                            extend_gap_reference;
    
        char                gap_char;
        /* gap open and extend character*/
    
        void                _init_alphabet (long not_found = -1);
        
};

class CawalignCodonScores : public CawalignSimpleScores {
    public:
    
        CawalignCodonScores  (ConfigParser*);
        virtual ~CawalignCodonScores (void) {};
    
        // compute how many nucleotides are different between the two codons encoded as 0-63 integers
        static int nucleotide_diff (long, long);
        
        Vector                translation_table;
        // codon (0-63 index) to single letter amino-acid code translation table
    
        // partial score tables
        VectorFP          s3x1,
                              s3x2,
                              s3x4,
                              s3x5;
        
        // the cost of introducing frameshits
        cawlign_fp                frameshift_cost,
        // the penalty for synonymous substitutions, per nucleotide change
                              synonymous_penalty;
        
        // ordered amino-acid scoring tables
        StringBuffer          amino_acids;
    
        int                   stop_codon_index;
        int                   mismatch_index;
    
    
    
};


extern const char   kNucleotideAlphabet[];
extern const cawlign_fp kNucScoring[];

#endif
