
#include "scoring.hpp"

#ifdef _OPENMP
    #include <omp.h>
#endif

using namespace std;
using namespace argparse;

const char   kNucleotideAlphabet[] = "ACGT";

const cawlign_fp kNucScoring[] = {5.,-4.,-4.,-4.,-5.,
                             -4., 5.,-1.,-4.,-5.,
                             -4.,-4., 5.,-4.,-5.,
                             -4.,-4.,-4., 5.,-5.,
                             -5.,-5.,-5.,-5.,1.};




void CawalignSimpleScores::_init_alphabet (long not_found) {
    for (int i = 0; i < 255; i++) {
        char_map [i] = not_found;
    }
    for (int i = 0; i < D; i++) {
        char_map [(unsigned char)alphabet.getChar(i)] = i;
    }
}

CawalignSimpleScores::CawalignSimpleScores (
                                      const char * _alphabet,
                                      const cawlign_fp * _scoring_matrix,
                                      const cawlign_fp _open_gap_reference,
                                      const cawlign_fp _open_gap_query,
                                      const cawlign_fp _extend_gap_reference,
                                      const cawlign_fp _extend_gap_query
                                      ) :
            alphabet(),
            scoring_matrix(),
            open_gap_reference       (_open_gap_reference),
            open_gap_query           (_open_gap_query),
            extend_gap_reference     (_extend_gap_reference),
            extend_gap_query         (_extend_gap_query),
            gap_char ('-')
{
    alphabet.appendBuffer(_alphabet);
    D = (unsigned int)alphabet.length();
    if (D == 0) {
        ERROR_NO_USAGE ("Empty alphabet");
    }
    _init_alphabet (D);
    scoring_matrix.appendValues(_scoring_matrix, (D+1)*(D+1));
        
}

CawalignSimpleScores::CawalignSimpleScores (
                                      ConfigParser * settings
                                      ) :
            alphabet(),
            scoring_matrix(),
            gap_char ('-')
{
    string _alph = settings->aConfig<string>("ALPHABET", "alphabet");
    alphabet.appendBuffer(_alph.c_str());
    D = (unsigned int)alphabet.length();
    if (D == 0) {
        ERROR_NO_USAGE ("Empty/missing alphabet");
    }
    _init_alphabet ();
    vector<cawlign_fp> _scores = settings->aConfigVec<cawlign_fp>("MATRIX", "cost");
    if (_scores.size() != (D+1) * (D+1)) {
        ERROR_NO_USAGE ("The dimension of the cost matrix is incorrect");
    }
    scoring_matrix.appendValues((cawlign_fp*)_scores.data(), (D+1)*(D+1));
    open_gap_reference  = settings->aConfig<cawlign_fp>("PARAMETERS", "open_deletion");
    open_gap_query      = settings->aConfig<cawlign_fp>("PARAMETERS", "open_insertion");
    extend_gap_reference  = settings->aConfig<cawlign_fp>("PARAMETERS", "extend_deletion");
    extend_gap_query      = settings->aConfig<cawlign_fp>("PARAMETERS", "extend_insertion");
}

int CawalignCodonScores::nucleotide_diff (long c1, long c2) {
    long diff = c1 ^ c2; // exclusive OR to set differences
    return ((diff & 0x03) ? 1 : 0) + ((diff & 0xC) ? 1 : 0) + ((diff & 0x30) ? 1 : 0);
}

StringBuffer codon_string (long c1) {
    StringBuffer codon;
    codon.appendChar(kNucleotideAlphabet[(c1 & 0x30) >> 4]);
    codon.appendChar(kNucleotideAlphabet[(c1 & 0x0C) >> 2]);
    codon.appendChar(kNucleotideAlphabet[(c1 & 0x03)]);
    return codon;
}

CawalignCodonScores::CawalignCodonScores (ConfigParser * settings) {
    alphabet.appendBuffer(kNucleotideAlphabet);
    
    D = 4;
    _init_alphabet ();
    D = 64;
    
    gap_char = '-';
    string _alph = settings->aConfig<string>("CODE", "aminoacids");
    
    const  long aaD = _alph.size();
    if (aaD < 21) {
        ERROR_NO_USAGE ("Incomplete amino-acid alphabet");
    }
    
    char allowed_aa [255] {-1};
    for (int i = 0; i < aaD; i++) {
        allowed_aa[toupper(_alph.at(i))] = i;
    }
    
    stop_codon_index = allowed_aa['X'];
    mismatch_index = allowed_aa['*'];
    
    if (stop_codon_index < 0) {
        ERROR_NO_USAGE ("Could not find the stop codon character 'X' (CODE:aminoacids)");
    }
    
    if (mismatch_index < 0) {
        ERROR_NO_USAGE ("Could not find the unresolved character '*' (CODE:aminoacids)");
    }
    
    vector<cawlign_fp> _scores = settings->aConfigVec<cawlign_fp>("MATRIX", "cost");
    if (_scores.size() != (aaD) * (aaD)) {
        ERROR_NO_USAGE ("The dimension of the cost matrix is incorrect (MATRIX:cost)");
    }
    
    cawlign_fp max_score = 0.;
    cawlign_fp min_score = 1e100;
    for (int i = 0; i < _scores.size(); i++) {
        cawlign_fp score = _scores[i];
        if (score > max_score) {
            max_score = score;
        }
        if (score < min_score) {
            min_score = score;
        }
    }
    
    vector<string> _translations = settings->aConfigVec<string>("CODE", "translations");
    
    if (_translations.size() != 64) {
        ERROR_NO_USAGE ("Expected a vector with 64 translations (CODE:translations)");
    }
    
    for (int i = 0; i < 64; i++) {
        const string token = _translations[i];
        if (token.size () != 1) {
            ERROR_NO_USAGE ("All entries in CODE:translations must have length 1");
        }
        const char translation_token = allowed_aa[token[0]];
        if (translation_token < 0) {
            ERROR_NO_USAGE ("All entries in CODE:translations must be present in CODE:aminoacids");
        }
        translation_table.appendValue(translation_token);
    }
    
    
    synonymous_penalty = 1.;
    /* first, define a 65x65 scoring matrix for all pairs of codons + sink (unresolved) state
     
     for resolved codons, the score is defined as the score of the corresponding amino-acid table
     for synonymous codons, the score is defined as the match score for the corresponding amino-acid match, with
     "synonymous_penalty" subtracted for each nucleotide substitution
     
     e.g. is Serine to Serine score is 7, then the scores for some of the serine codon pairs will be like this
     
     AGC : AGC = 7
     AGC : AGT = 7 - synonymous_penalty
     AGC : TCT = 7 - 3 x synonymous_penalty
     */
    
    const cawlign_fp _large_penalty = -1.e4;
    for (long codon1 = 0; codon1 < 64; codon1++) {
        long codon1_translation = translation_table.value(codon1);
        for (long codon2 = 0; codon2 < 64; codon2++) {
            long codon2_translation = translation_table.value(codon2);
            if ((codon1_translation == stop_codon_index || codon2_translation == stop_codon_index) && codon2_translation != codon1_translation) {
                scoring_matrix.appendValue (_large_penalty);
            } else {
                cawlign_fp score = _scores[ codon1_translation*aaD + codon2_translation ] ;
                //if (codon2_translation == codon1_translation) {
                    if (codon1 != codon2) {
                        score -= 0.5;
                    }
                //}
                scoring_matrix.appendValue (score);
            }
            //printf ("%g ", scoring_matrix.value (scoring_matrix.length()-1));
        }
        // score codon1 vs "unresolved"
        scoring_matrix.appendValue (0.);
        //printf ("%g\n", scoring_matrix.value (scoring_matrix.length()-1));
    }
    
    // add unresolved vs codon
    for (long codon1 = 0; codon1 < 64; codon1++) {
        //long codon1_translation = translation_table.value(codon1);
        scoring_matrix.appendValue (0.);
    }
    
    // add unresolved vs unresolved
    scoring_matrix.appendValue (_scores[mismatch_index + aaD*mismatch_index]);
    //printf ("\n%ld %ld %g", 64, 64, scoring_matrix.value (scoring_matrix.length()-1));
    
    /* populate various partial scoring matrices*/
    
    /* c3x1 is the (65x12) matrix  which scores matching three characters in the reference with a single character in the query, e.g.
     
     AAA  CGC
     -A-  C--
     
     The scoring for this goes as follows. For each of the codons (+ a sink state for all unresolved codons, e.g. NNA)
     compute 12 scores, each possible nucleotide in each possible position, indexed as such
     
     0:      A--
     1:      -A-
     2:      --A
     3:      C--
     ...
     11 :    --T
     
     The score is obtained by taking all 16 possible filler resolutions, e.g.
     
     score (codon X, A--) = Max (score (codon X, Aaa), score (codon X, Aac), score (codon X, Aag), .. score (codon X, Att))
     the score for a pair of codons
     
     c3x2 is the (65x48) matrix which scores matching three characters in the reference with a two character in the query, e.g.
     
     AAA  CGC
     -AA  C-C
     
     Scoring follows the same pattern as in c3x1, i.e.
     
     score (codon X, AA-) = Max (score (codon X, AAa), score (codon X, AAc), ... score (codon X, AAt))
     
     the ordering of the scores in each row is like this
     
     0: AA-
     1: A-A
     2: -AA
     3: AC-
     4: A-C
     5 :-AC
     ...
     
     
     */
    
    for ( long thisCodon = 0; thisCodon < 64; thisCodon ++ ) {
        for ( long d1 = 0; d1 < 4; d1 ++ ) {
            cawlign_fp max100 = -1e100;
            cawlign_fp max010 = -1e100;
            cawlign_fp max001 = -1e100;
            
            for ( long d2 = 0; d2 < 4; d2 += 1 ) {
                long partialCodon = 4 * d1 + d2;
                cawlign_fp max110 = -1e100;
                cawlign_fp max101 = -1e100;
                cawlign_fp max011 = -1e100;
                
                for ( long d3 = 0; d3 < 4; d3 += 1 ) {
                    long thisCodon2 = 4 * partialCodon + d3;
                    cawlign_fp thisScore = scoring_matrix.value(thisCodon * 65 +  thisCodon2);
                    
                    // this is the trivial and stupid way of doing it, but it should work
                    
                    for ( long i = 0; i < 10; i++) {
                        s3x5.appendValue (thisScore);
                    }
                    
                    for ( long i = 0; i < 4; i++) {
                        s3x4.appendValue (thisScore);
                    }
                    
                    // d1 is 1
                    max100 = MAX( max100, scoring_matrix.value( thisCodon * 65 + 16 * d1 + 4 * d2 + d3  ));
                    max010 = MAX( max010, scoring_matrix.value( thisCodon * 65 + 16 * d2 + 4 * d1 + d3  ));
                    max001 = MAX( max001, scoring_matrix.value( thisCodon * 65 + 16 * d2 + 4 * d3 + d1  ));
                    
                    // d1 and d2 are 1
                    max110 = MAX( max110, scoring_matrix.value( thisCodon * 65 + 16 * d1 + 4 * d2 + d3  ));
                    max101 = MAX( max101, scoring_matrix.value( thisCodon * 65 + 16 * d1 + 4 * d3 + d2  ));
                    max011 = MAX( max011, scoring_matrix.value( thisCodon * 65 + 16 * d3 + 4 * d1 + d2  ));
                }
                
                s3x2.appendValue (max110);
                s3x2.appendValue (max101);
                s3x2.appendValue (max011);
                
            }
            
            s3x1.appendValue (max100);
            s3x1.appendValue (max010);
            s3x1.appendValue (max001);
        }
    }
    
    // pad all thhe partial score matrices with "0" for the unresolved codons
    auto pad_vector = [] (VectorFP & v, long how_many) {
        for (long i = 0; i < how_many; i++) {
            v.appendValue(0.);
        }
    };
    
    pad_vector (s3x1, 10);
    pad_vector (s3x2, 48);
    pad_vector (s3x4, 256);
    pad_vector (s3x5, 640);
    
    
    open_gap_reference    = settings->aConfig<cawlign_fp>("PARAMETERS", "open_deletion");
    open_gap_query        = settings->aConfig<cawlign_fp>("PARAMETERS", "open_insertion");
    extend_gap_reference  = settings->aConfig<cawlign_fp>("PARAMETERS", "extend_deletion");
    extend_gap_query      = settings->aConfig<cawlign_fp>("PARAMETERS", "extend_insertion");
    frameshift_cost       = settings->aConfig<cawlign_fp>("PARAMETERS", "frameshift_cost");
    
    cawlign_fp indel_cost = MAX(max_score, -min_score),
               ext_cost = 3.*(max_score-min_score) / 40.;
    
    if (frameshift_cost < 0.) {
        frameshift_cost = 3.*indel_cost;
    }
    if (open_gap_reference < 0.) {
        open_gap_reference = 2.*indel_cost;
    }
    if (open_gap_query < 0.) {
        open_gap_query = 2.*indel_cost;
    }
    
    if (extend_gap_query < 0.) {
        extend_gap_query = ext_cost;
    }
    
    if (extend_gap_reference< 0.) {
        extend_gap_reference = ext_cost;
    }
}
