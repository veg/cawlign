/*
 
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
 Art FY Poon    (apoon42@uwo.ca)
 Steven Weaver (sweaver@temple.edu)
 
 Module Developers:
 Lance Hepler (nlhepler@gmail.com)
 Martin Smith (martin.audacis@gmail.com)
 
 Significant contributions from:
 Spencer V Muse (muse@stat.ncsu.edu)
 Simon DW Frost (sdf22@cam.ac.uk)
 
 Permission is hereby granted, free of charge, to any person obtaining a
 copy of this software and associated documentation files (the
 "Software"), to deal in the Software without restriction, including
 without limitation the rights to use, copy, modify, merge, publish,
 distribute, sublicense, and/or sell copies of the Software, and to
 permit persons to whom the Software is furnished to do so, subject to
 the following conditions:
 
 The above copyright notice and this permission notice shall be included
 in all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 
 */

#include <cctype>
#include <cstring>


#include <math.h>
#include <stdio.h>

#include "alignment.h"

#define MAX_OP(a, b) ( (( a ) > ( b )) ? ( a ) : ( b ) )

//____________________________________________________________________________________

/**
    The following numerical codes define "moves" available to the algorithm working in codon space for pairwise alignment.
    The first (reference( sequence (represented by ROWS in the dynamic programming matrix) is the "blessed" reference sequence, i.e. it's in frame (no stop codons, divisble by 3)
    The second (query) sequence (represented by COLUMNS in the dynamic programming matrix) is the one being mapped to the reference, and it admits out-of-frame indels
 
    The mnemonic for each macro is as follows
        HY_XXX_YYY
                Len (X) = Len (Y) = the number of characters output in the final alignment
                Each '1' in the pattern corresponds to an actual alignment character being consumed from the input, and '0' - to a gap (indel)
 
        For example, i.f the pattern is
 
            HY_1101_1111, and the current context of the strings being aligned is
 
        REF: ... AGTACA ...
        QRY: ...AAGTACA...
 
            then the output will receive
        A-GT
        AaGT
            and the counter will advance 3 characters in the reference and 4 characters in the query
 
*/


#define HY_ALIGNMENT_TYPES_COUNT 24

// match or skip whole codons
#define HY_111_111 0
#define HY_111_000 1
#define HY_000_111 2

// match 3 in the ref to 1 in the query
#define HY_111_100 3
#define HY_111_010 4
#define HY_111_001 5

#define HY_3X1_START 3
#define HY_3X1_COUNT 3

// match 3 in the ref to 2 in the query
#define HY_111_110 6
#define HY_111_101 7
#define HY_111_011 8

/** 
 all operations consuming 3 chars in the reference and 2 in the query start with index 6 and there are 3 of them
*/
#define HY_3X2_START 6
#define HY_3X2_COUNT 3

// match 3 in the ref to 4 in the query
#define HY_1110_1111 9
#define HY_1101_1111 10
#define HY_1011_1111 11
#define HY_0111_1111 12

/**
 all operations consuming 3 chars in the reference and 4 in the query start with index 9 and there are 4 of them
*/
#define HY_3X4_START 9
#define HY_3X4_COUNT 4

// match 3 in the ref to 5 in the query

#define HY_11100_11111 13
#define HY_11010_11111 14
#define HY_11001_11111 15
#define HY_10110_11111 16
#define HY_10101_11111 17
#define HY_10011_11111 18
#define HY_01110_11111 19
#define HY_01101_11111 20
#define HY_01011_11111 21
#define HY_00111_11111 22

/**
 all operations consuming 3 chars in the reference and 5 in the query start with index 13 and there are 10 of them
*/
#define HY_3X5_START   13
#define HY_3X5_COUNT   10


// the local align move (i.e. take a direct shortcut to any internal position in the alignment matrix)

#define HY_LOCAL_ALIGN_SHORTCUT 23

//____________________________________________________________________________________


/**
 * @name CodonAlignStringsStep
 * Perform a single step in the codon alignment algorithm; returns the best operation and updates DP matrices
 *
 * @param score_matrix an NxN score matrix (N = 65, to accommodate all 64 codons plus an explicit "unresolved" character)
 * @param reference the reference string remapped with characters remapped to indices in score_matrix
 * @param query the query string remapped with characters remapped to indices in score_matrix
 * @param r index in reference (counted out in codons)
 * @param q index in query (counted out in nucleotides)
 * @param score_cols the number of columns in score_matrix
 * @param char_count the number of valid characters (nucleotides)
 * @param miscall_cost the cost of introducing an out-of-frame indel
 * @param open_deletion the cost of opening an indel in the reference 
 * @param extend_deletion the cost of extending an indel in the reference
 * @param open_inserion the cost of opening an indel in the query
 * @param extend_inserion the cost of extending an indel in the query
 * @param cost_matrix the main DP matrix
 * @param cost_stride the number of columns in the DP matrix
 * @param insertion_matrix the DP matrix for affine insertions
 * @param deletion_matrix the DP matrix for affine deletions
 * @param codon3x5 the 3-ref 5-qry scoring matrix
 * @param codon3x4 the 3-ref 4-qry scoring matrix
 * @param codon3x2 the 3-ref 2-qry scoring matrix
 * @param codon3x1 the 3-ref 1-qry scoring matrix
 * @param do_local if TRUE, perform a local alignment (no prefix/suffix indel cost)
 *
 * @return the best scoring alignment operation

 */

long CodonAlignStringsStep( cawlign_fp * const score_matrix
                          , long * const reference
                          , long * const query
                          , const long r
                          , const long q
                          , const long score_cols
                          , const long char_count
                          , const cawlign_fp miscall_cost
                          , const cawlign_fp open_insertion
                          , const cawlign_fp open_deletion
                          , const cawlign_fp extend_insertion
                          , const cawlign_fp extend_deletion
                          , const cawlign_fp * cost_matrix
                          , const long cost_stride
                          , cawlign_fp * const insertion_matrix
                          , cawlign_fp * const deletion_matrix
                          , cawlign_fp const *codon3x5
                          , cawlign_fp const *  codon3x4
                          , cawlign_fp const *  codon3x2
                          , cawlign_fp const*  codon3x1
                          , const    bool  do_local
                          )
{
    /**
     * r is CODON position in the reference,
     * q is NUC position in the query,
     * curr is a pointer to the current position in the scoring matrix,
     * prev is a pointer to the previous CODON in the scoring matrix
     * NOTE: moving by score_cols in the scoring matrix changes the CODON
     *       position in the scoring matrix, as we're only interested in CODON
     *       alignments to the reference
     * rpos is the position of r in the reference
     * do_local is true if we wish to perform a local alignment 
     */
    const long curr = ( r ) * score_cols + q, // where we currently are in the DP matrix
               prev = curr - score_cols, // upstream a codon in the reference for the DP matrix
               // offsets are strides in the corresponding scoring matrices
               offset3x5 = HY_3X5_COUNT * char_count * char_count * char_count, // both 3x5 and 3x4 are
               offset3x4 = HY_3X4_COUNT * char_count * char_count * char_count, // full codons
               offset3x2 = HY_3X2_COUNT * char_count * char_count,
               offset3x1 = HY_3X1_COUNT * char_count,
               rpos = r * 3; // the real position in R
    
    // mutable vars
    long r_codon = -1,
         q_codon = -1,
         best_choice = 0,
         //i,
         choice,
         partial_codons[ 10 ];
    // we need to multiply by 3 to get the NUC position
    
    // 3x5 codon specifications (inverted indices)
    static long const codon_spec_3x5[ 10 ][ 3 ] = {
        { 5, 4, 3 }, // 11100
        { 5, 4, 2 }, // 11010
        { 5, 4, 1 }, // 11001
        { 5, 3, 2 }, // 10110
        { 5, 3, 1 }, // 10101
        { 5, 2, 1 }, // 10011
        { 4, 3, 2 }, // 01110
        { 4, 3, 1 }, // 01101
        { 4, 2, 1 }, // 01011
        { 3, 2, 1 }  // 00111
    };
    // 3x4 codon specifications (inverted indices)
    static long const long codon_spec_3x4[ 4 ][ 3 ] = {
        { 4, 3, 2 }, // 1110
        { 4, 3, 1 }, // 1101
        { 4, 2, 1 }, // 1011
        { 3, 2, 1 }  // 0111
    };

    long    local_shortcut_came_from_this_move = -1;

    cawlign_fp choices[ HY_ALIGNMENT_TYPES_COUNT ],
           max_score = -INFINITY,
           penalty;

    // store the scores of our choices in choices,
    // pre-initialize to -infinity
    for (long i = 0; i < HY_ALIGNMENT_TYPES_COUNT; i++ ) {
        choices[ i ] = -INFINITY;
    }
    
    // if we're at least a CODON away from the edge...
    // (psst, r is CODONs remember?)
    if ( r >= 1 ) {
        // if we're doing affine gaps (deletions)
        if ( deletion_matrix ) {
            choices[ HY_111_000 ] = MAX_OP(
                score_matrix[ prev ] - open_deletion,
                deletion_matrix[ prev ] - ( r > 1 ? extend_deletion : open_deletion )
            );
            deletion_matrix[ curr ] = choices[ HY_111_000 ];
        } else {
            choices[ HY_111_000 ] = score_matrix[ prev ] - open_deletion;
        }

        r_codon = ( reference[ rpos - 3 ]   * char_count
                  + reference[ rpos - 2 ] ) * char_count
                  + reference[ rpos - 1 ] ;

        if ( r_codon < 0 ) {
            // anything other than a fully resolved codon gets mapped to the "generic unresolved" character
            r_codon = cost_stride - 1;
            /* 20240219: SLKP optimization note this needs to be adjusted; let's not print stuff to stderr */
            //fprintf (stderr, "*** NIL CODON *** %ld %ld %ld %ld\n", reference[ rpos - 3 ], reference[ rpos - 2], reference[ rpos - 1], r_codon);
        }
    }

    // if we're at least 1 codon away from the edge
    if ( q >= 3 ) {
        // if we're doing affine gaps (insertions)
        if ( insertion_matrix ) {
            choices[ HY_000_111 ] = MAX_OP(
                score_matrix[ curr - 3 ] - open_insertion,
                insertion_matrix[ curr - 3 ] - ( q > 3 ? extend_insertion : open_insertion )
            );
            insertion_matrix[ curr ] = choices[ HY_000_111 ];
        } else {
            choices[ HY_000_111 ] = score_matrix[ curr - 3 ] - open_insertion;
        }

        q_codon = ( query[ q - 3 ]   * char_count
                  + query[ q - 2 ] ) * char_count
                  + query[ q - 1 ] ;

        if ( q_codon < 0 ) {
           // anything other than a fully resolved codon gets mapped to the "generic unresolved" character
           q_codon = cost_stride - 1;
        }
    }

    // if q_codon and r_codon both exist, set the score equal to match
    if ( q_codon >= 0 && r_codon >= 0 ) {
        const cawlign_fp move_cost = cost_matrix[ r_codon * cost_stride + q_codon ];
        choices[ HY_111_111 ] = score_matrix[ prev - 3 ] + move_cost;
        //if (do_local && choices [HY_LOCAL_ALIGN_SHORTCUT] < move_cost) {
        //    local_shortcut_came_from_this_move = HY_111_111;
        //    choices [HY_LOCAL_ALIGN_SHORTCUT]  = move_cost;
        //}
    }

    // we disallow partial moves in the reference, so those used to be here but are now gone

    // HERE BE DRAGONS!!!!

    // miscall matches, starting with 3x5, then 3x4, then 3x2, finally 3x1
    if ( r_codon >= 0 ) {
        // 3x5 partial codons
        if ( q >= 5 ) {
            // fill in the partial codons table
            // use a 10x1 array to allow for load hoisting,
            // we don't want to be bouncing cachelines in this critical inner loop
            for (long i = 0; i < HY_3X5_COUNT; ++i ) {
                partial_codons[ i ] = ( query[ q - codon_spec_3x5[ i ][ 0 ] ]   * char_count
                                      + query[ q - codon_spec_3x5[ i ][ 1 ] ] ) * char_count
                                      + query[ q - codon_spec_3x5[ i ][ 2 ] ] ;
            }
            // go over each choice, fill it in
            for (long i = 0; i < HY_3X5_COUNT; ++i ) {
                if ( partial_codons[ i ] >= 0 ) {
                    // this partial codon is resolved
                    choice = HY_3X5_START + i;
                    // if we have a cawlign_fp ragged edge, don't penalize
                    if ( ( q == 5              && choice == HY_00111_11111 )
                      || ( q == score_cols - 1 && choice == HY_11100_11111 ) )
                        penalty = 0.;
                    // if we have a single ragged edge, penalize by a single miscall
                    // we don't have to worry about specifying each case here,
                    // as the 00111_11111 case takes preference above,
                    // so we don't have to explicitly avoid it
                    else if ( q == 5 && choice >= HY_01110_11111 )
                        penalty = miscall_cost;
                    // if we have a single ragged edge, penalize by a single miscall
                    // unfortunately these cases are spread out,
                    // so we have to enumerate them explicitly here
                    else if ( q == score_cols - 1
                           && ( choice == HY_11010_11111
                             || choice == HY_10110_11111
                             || choice == HY_01110_11111 ) )
                        penalty = miscall_cost;
                    // if we don't fit into any of these special cases,
                    // the miscall penalty is cawlign_fp (as we're matching 3 to 5)
                    else
                        penalty = 2. * miscall_cost;
                        
                    const cawlign_fp move_cost = codon3x5[ r_codon * offset3x5 + HY_3X5_COUNT * partial_codons[ i ] + i ];
                    
                    choices[ choice ] = score_matrix[ prev - 5 ] - penalty + move_cost;
                    //if (do_local && choices [HY_LOCAL_ALIGN_SHORTCUT] < move_cost) {
                    //    local_shortcut_came_from_this_move = choice;
                    //    choices [HY_LOCAL_ALIGN_SHORTCUT]  = move_cost;
                    //}
                }
            }
        }

        // 3x4 partial codons
        if ( q >= 4 ) {
            // fill in partial codons table
            for (long i = 0; i < HY_3X4_COUNT; ++i ) {
                partial_codons[ i ] = ( query[ q - codon_spec_3x4[ i ][ 0 ] ]   * char_count
                                      + query[ q - codon_spec_3x4[ i ][ 1 ] ] ) * char_count
                                      + query[ q - codon_spec_3x4[ i ][ 2 ] ] ;
            }
            // fill in choices
            for (long i = 0; i < HY_3X4_COUNT; ++i ) {
                if ( partial_codons[ i ] >= 0 ) {
                    choice = HY_3X4_START + i;
                    // if we have a ragged edge,
                    // penalize it not at all
                    if ( ( q == 4              && choice == HY_0111_1111 )
                      || ( q == score_cols - 1 && choice == HY_1110_1111 ) )
                        penalty = 0.;
                    // otherwise it's just a single miscall penalty
                    else
                        penalty = miscall_cost;

                    const cawlign_fp move_cost = codon3x4[ r_codon * offset3x4 + HY_3X4_COUNT * partial_codons[ i ] + i ];

                    choices[ choice ] = score_matrix[ prev - 4 ] - penalty + move_cost;
                    //if (do_local && choices [HY_LOCAL_ALIGN_SHORTCUT] < move_cost) {
                    //    local_shortcut_came_from_this_move = choice;
                    //    choices [HY_LOCAL_ALIGN_SHORTCUT]  = move_cost;
                    //}
                }
            }
        }

        // 3x2
        if ( q >= 2 ) {
            // only a single partial codon
            partial_codons[ 0 ] = query[ q - 2 ] * char_count
                                + query[ q - 1 ] ;
            // fill in choices
            if ( partial_codons[ 0 ] >= 0 ) {
                for (long i = 0; i < HY_3X2_COUNT; ++i ) {
                    choice = HY_3X2_START + i;
                    // if we have a ragged edge at the beginning or end,
                    // respectively, don't penalize it
                    if ( ( q == 2              && choice == HY_111_011 )
                      || ( q == score_cols - 1 && choice == HY_111_110 ) )
                        penalty = 0.;
                    // otherwise it's just a single miscall penalty
                    else
                        penalty = miscall_cost;
                        
                                                                
                    const cawlign_fp move_cost = codon3x2[ r_codon * offset3x2 + HY_3X2_COUNT * partial_codons[ 0 ] + i ];

                    choices[ choice ] = score_matrix[ prev - 2 ] - penalty + move_cost;
                    //if (do_local && choices [HY_LOCAL_ALIGN_SHORTCUT] < move_cost) {
                    //    local_shortcut_came_from_this_move = choice;
                    //    choices [HY_LOCAL_ALIGN_SHORTCUT]  = move_cost;
                    //}
                }
            }
        }

        // 3x1
        if ( q >= 1 ) {
            // only a single partial codon
            partial_codons[ 0 ] = query[ q - 1 ];
            // fill in choices
            if ( partial_codons[ 0 ] >= 0 ) {
                for (long i = 0; i < HY_3X1_COUNT; ++i ) {
                    choice = HY_3X1_START + i;
                    // if we have a cawlign_fp ragged edge,
                    // don't enforce a miscall penalty
                    if ( ( q == 1              && choice == HY_111_001 )
                      || ( q == score_cols - 1 && choice == HY_111_100 ) )
                        penalty = 0.;
                    // if we have a single ragged edge,
                    // enforce only a single miscall penalty
                    else if ( ( q == 1              && choice == HY_111_010 )
                           || ( q == score_cols - 1 && choice == HY_111_010 ) )
                        penalty = miscall_cost;
                    // otherwise we need a cawlign_fp miscall penalty,
                    // for the two positions we're inserting
                    else
                        penalty = 2. * miscall_cost;

                    const cawlign_fp move_cost = codon3x1[ r_codon * offset3x1 + HY_3X1_COUNT * partial_codons[ 0 ] + i ];

                    choices[ choice ] = score_matrix[ prev - 1 ] - penalty + move_cost;
                    //if (do_local && choices [HY_LOCAL_ALIGN_SHORTCUT] < move_cost) {
                    //    local_shortcut_came_from_this_move = choice;
                    //    choices [HY_LOCAL_ALIGN_SHORTCUT]  = move_cost;
                    //}
                }
            }
        }
    }

    // find the best possible choice
    //if (do_local) {
        for (long i = 0; i < HY_ALIGNMENT_TYPES_COUNT ; ++i ) {
            if ( choices[ i ] > max_score ) {
                best_choice = i;
                max_score = choices[ i ];
            }
        }
    /*} else {
        for (long i = 0; i < HY_ALIGNMENT_TYPES_COUNT - 1 ; ++i ) {
            if ( choices[ i ] > max_score ) {
                best_choice = i;
                max_score = choices[ i ];
            }
        }
    }*/
    
    
    //fprintf( stderr, "\nscore: %.3g best: %ld\n", max_score, best_choice );

    // assign the score to the current position
    score_matrix[ curr ] = max_score;

    //if (do_local && best_choice == HY_LOCAL_ALIGN_SHORTCUT) {
    //    return -local_shortcut_came_from_this_move - 1;
    //}
    return best_choice;
}

//____________________________________________________________________________________

/**
 * @name BacktrackAlign
 * Perform a single backtracking step in the dynamic programming matrix for single character alignment (nucleotides or amino-acids)
 * @param edit_ops an array of chars (only 3 states used) to record which of the three possible steps is taken
 *             0  : match    (ref: X, qry : Y)
 *             -1 : deletion (ref: X, qry : -)
 *             1  : insterion (ref: -, qry: Y)
 * @param edit ptr: current index into the edit_ops array (gets ++ with each step)
 * @param r       current character index in the reference
 * @param q       current character index in the query
 * @param deletion cost of the deletion (-1) op
 * @param insertion cost of the insertion (-1) op
 * @param match cost of the match  (0) op
 */

inline void BacktrackAlign( signed char * const edit_ops
                          , long & edit_ptr
                          , long & r
                          , long & q
                          , cawlign_fp deletion
                          , cawlign_fp insertion
                          , cawlign_fp match
                          ) {
    if ( match >= deletion && match >= insertion ) {
        --r;
        --q;
        edit_ops[ edit_ptr++ ] = 0;
    } else if ( deletion >= insertion ) {
        --r;
        edit_ops[ edit_ptr++ ] = -1;
    } else {
        --q;
        edit_ops[ edit_ptr++ ] = 1;
    }
}

//____________________________________________________________________________________

/**
 * @name BacktrackAlignCodon
 * Perform a single backtracking step in the dynamic programming matrix for codon alignment
 * @param edit_ops an array of chars (only 3 states used) to record which of the three possible steps is taken
 *             0  : match    (ref: X, qry : Y)
 *             -1 : deletion (ref: X, qry : -)
 *             -2 : deletion with frameshift (ref Xx, qry Y-)
 *             1  : insterion (ref: -, qry: Y)
 *             2 : deletion with frameshift (ref X-, qry Yy)
 * @param edit ptr: current index into the edit_ops array (gets ++ with each step)
 * @param r       current character index in the reference
 * @param q       current character index in the query
 * @param code  operation code; the operation yielding the maximum score for this step
*/

inline void BacktrackAlignCodon( signed char * const edit_ops
                               , long & edit_ptr
                               , long & r
                               , long & q
                               , const long code
                               )
{
    /** SLKP 20240319:
          Optimization notes (if warranted)
                1. this switch statement needs to be profiled to see if the ordering of the cases
                2. there's also an opportinity to use bitmasks
            
    */
    
    switch ( code ) {
            // match
        case HY_111_111:
            r -= 3; q -= 3;
            edit_ops [edit_ptr]   = 0;
            edit_ops [edit_ptr+1] = 0;
            edit_ops [edit_ptr+2] = 0;
            edit_ptr += 3;
            return;
            
            // deletion
        case HY_111_000:
            r -= 3;
            edit_ops [edit_ptr]   = -1;
            edit_ops [edit_ptr+1] = -1;
            edit_ops [edit_ptr+2] = -1;
            edit_ptr += 3;
           return;
            
            // insertion
        case HY_000_111:
            q -= 3;
            edit_ops [edit_ptr]   = 1;
            edit_ops [edit_ptr+1] = 1;
            edit_ops [edit_ptr+2] = 1;
            edit_ptr += 3;
            return;
    }
    
    unsigned char r_str[ 5 ] = { 1, 1, 1, 1, 1 },
         // which characters are we taking from the reference (up to 5); 1 take, 0 leave
            q_str[ 5 ] = { 1, 1, 1, 1, 1 };
         // which characters are we taking from the query (up to 5)
    int     idx         = 2;
        // how many characters are we taking (max over both strings); index so the actual number is +1
        // can be 2,3, or 4


    switch (code) {
        // 3x2
        case HY_111_110:
            q_str[2] = 0;
            break;
        case HY_111_101:
            q_str[1] = 0;
            break;
        case HY_111_011:
            q_str[0] = 0;
            break;

        // 3x1
        case HY_111_100:
            q_str[1] = 0;
            q_str[2] = 0;
            break;
        case HY_111_010:
            q_str[0] = 0;
            q_str[2] = 0;
            break;
        case HY_111_001:
            q_str[0] = 0;
            q_str[1] = 0;
            break;

        // 3x4
        case HY_1110_1111:
            r_str[3] = 0;
            idx = 3;
            break;
        case HY_1101_1111:
            r_str[2] = 0;
            idx = 3;
            break;
        case HY_1011_1111:
            r_str[1] = 0;
            idx = 3;
            break;
        case HY_0111_1111:
            r_str[0] = 0;
            idx = 3;
            break;

        // 3x5
        case HY_11100_11111:
            r_str[3] = 0;
            r_str[4] = 0;
            idx = 4;
            break;
        case HY_11010_11111:
            r_str[2] = 0;
            r_str[4] = 0;
            idx = 4;
            break;
        case HY_11001_11111:
            r_str[2] = 0;
            r_str[3] = 0;
            idx = 4;
            break;
        case HY_10110_11111:
            r_str[1] = 0;
            r_str[4] = 0;
            idx = 4;
            break;
        case HY_10101_11111:
            r_str[1] = 0;
            r_str[3] = 0;
            idx = 4;
            break;
        case HY_10011_11111:
            r_str[1] = 0;
            r_str[2] = 0;
            idx = 4;
            break;
        case HY_01110_11111:
            r_str[0] = 0;
            r_str[4] = 0;
            idx = 4;
            break;
        case HY_01101_11111:
            r_str[0] = 0;
            r_str[3] = 0;
            idx = 4;
            break;
        case HY_01011_11111:
            r_str[0] = 0;
            r_str[2] = 0;
            idx = 4;
            break;
        case HY_00111_11111:
            r_str[0] = 0;
            r_str[1] = 0;
            idx = 4;
            break;
    }

    for ( ; idx >= 0 ; --idx ) {
        if ( r_str[ idx ] ) {
            if ( q_str[ idx ] ) {
                --r;
                --q;
                edit_ops[ edit_ptr++ ] = 0;
            } else {
                --r;
                edit_ops[ edit_ptr++ ] = -2;
            }
        } else {
            --q;
            edit_ops[ edit_ptr++ ] = 2;
        }
    }
}

//____________________________________________________________________________________

/**
 * @name MatchScore
 * Looks up the score for matching two characters at the current positions in REF and QRY sequences, assuming they are valid
 * Increments the score parameter by that value (which of course can also be negative)
 * If either character is invalid, the implied score is 0
 *
 * @param r_str the reference string
 * @param q_str the query string
 * @param r index in the reference (1-based)
 * @param q index in the query (1-based)
 * @param char_map a map from valid character to the encoding used by the score matrix
 * @param cost_matrix the pairwise character cost matrix
 * @param cost_stride the column dimension of the cost matrix
 * @param score the cumulative alignment score (will be modified)

*/

//____________________________________________________________________________________

inline void MatchScore( char const * r_str
                      , char const * q_str
                      , const long r
                      , const long q
                      , long * char_map
                      , const cawlign_fp * cost_matrix
                      , const long cost_stride
                      , cawlign_fp & score
                      )
{
    const long r_char = char_map[ (unsigned char)r_str[ r - 1 ] ];
    if ( r_char >= 0 ) {
        const long q_char = char_map[ (unsigned char) q_str[ q - 1 ] ];
        if ( q_char >= 0 ) {
            score += cost_matrix[ r_char * cost_stride + q_char ];
        }
    }
}

//____________________________________________________________________________________


/**
 * @name AlignStrings
 * Performs a pairwise string alignment; the nature of the alignment algoritm invoked will depend on option flags
 * returns the score of the alignment
 * r_res and q_res parameters will the receive aligned strings they will be allocated and populated by this function
 * 
 * @param r_str the reference string
 * @param q_str the query string
 * @param _r_len the length of the reference string (-1 to compute here)
 * @param _q_len the length of the query string (-1 to compute here)
 * @param r_res will receive the resulting REFERENCE alignment
 * @param q_res will receive the resulting QUERY alignment
 * @param char_map a char->long map which takes a character from the input string and maps into in [-1,N-1] where N is the dimension of the scoring matrix; an "indvalid" character is mapped to -1
 * @param cost_matrix an NxN (row-major) matrix of pairwise character scores
 * @param cost_stride N (column dimension of cost_matrix)
 * @param gap the character to use for an indel (e.g. '-')
 * @param open_insertion the cost of opening a gap in the REFERENCE
 * @param extend_insertion the cost of extending a gap in the REFERENCE
 * @param open_deletion the cost of opening a gap in the QUERY
 * @param extend_deletion the cost of extending a gap in the QUERY
 * @param miscall_cost the cost of introducing any out-of-frame indels
 * @param do_local if TRUE, local alignment will be performed (no prefix or suffix gap cost)
 * @param do_affine if TRUE, use the affine gap penalty
 * @param do_codon if TRUE, do a codon-aware alignment
 * @param char_count number of alphabet (N-1, where N == cost_stride)
 * @param codon3x5 (do_codon == TRUE), the cost of aligning codon X with any of the 10 options where 3 refs align to 5 qry
 * @param codon3x4 (do_codon == TRUE), the cost of aligning codon X with any of the 4 options where 3 refs align to 4 qry
 * @param codon3x2 (do_codon == TRUE), the cost of aligning codon X with any of the 3 options where 3 refs align to 2 qry
 * @param codon3x1 (do_codon == TRUE), the cost of aligning codon X with any of the 3 options where 3 refs align to 1 qry
 * @param do_true_local if TRUE, perform a true local (best scoring substings) alignment
 * @param report_ref_insertions if TRUE (default), will report deletions relative to the
 * @param score_matrix_cache if provided, use this to store the scoring matrix (assumed to have sufficient size)
 * @param insertion_matrix_cache if provided, use this to store the insertion (affine gaps) score matrix (assumed to have sufficient size)
 * @param deletion_matrix_cache if provided, use this to store the  deletion (affine gaps) score scoring matrix (assumed to have sufficient size)
 * @return the alignment score

*/
cawlign_fp AlignStrings( char const * r_str
                   , char const * q_str
                   , const long _r_len
                   , const long _q_len
                   , char * & r_res
                   , char * & q_res
                   , long * char_map
                   , cawlign_fp const * cost_matrix
                   , const long cost_stride
                   , const char gap
                   , cawlign_fp open_insertion
                   , cawlign_fp extend_insertion
                   , cawlign_fp open_deletion
                   , cawlign_fp extend_deletion
                   , cawlign_fp miscall_cost
                   , const bool do_local
                   , const bool do_affine
                   , const bool do_codon
                   , const long char_count
                   , const cawlign_fp * codon3x5
                   , const cawlign_fp * codon3x4
                   , const cawlign_fp * codon3x2
                   , const cawlign_fp * codon3x1
                   , const bool do_true_local
                   , const bool report_ref_insertions
                   , cawlign_fp* score_matrix_cache
                   , cawlign_fp* insertion_matrix_cache
                   , cawlign_fp* deletion_matrix_cache
                   )
{
    const unsigned long r_len = _r_len >= 0 ? _r_len : strlen( r_str ),
                        q_len = _q_len >= 0 ? _q_len : strlen( q_str ),
                        ref_stride = ( do_codon ? 3 : 1 ),
                        // when doing codon alignments, the algorithm will advance in the reference by full codons
                        score_rows = r_len / ref_stride + 1,
                        // R = number of rows in the dynamic programming matrix
                        score_cols = q_len + 1;
                        // C = number of columns in the dynamic programming matrix
                        // the extra row/column is to accommodate prefix gaps
                        

    
    if ( do_codon && ( r_len % 3 != 0 ) ) {
        return -INFINITY;
    }

    cawlign_fp score = 0.;

    // handle some edge cases,
    // return early if possible
    if ( score_rows <= 1 ) {
        // the reference string is EMPTY
        if ( score_cols > 1 ) {
            // the query string is NON-EMPTY
            r_res = new char[ q_len + 1 ];
            q_res = new char[ q_len + 1 ];
            // no ref, just query, which remains untouched
            memcpy( q_res, q_str, q_len + 1 );
            // ref full of gaps
            memset( r_res, gap, sizeof( char ) * q_len );
            // null terminate
            r_res[ q_len ] = '\0';
            q_res[ q_len ] = '\0';
            // compute the score for this "all-gaps" alignment
            if ( ! do_local ) {
                if ( do_affine )
                    score = -open_insertion - ( q_len - 1 ) * extend_insertion;
                else
                    score = -open_insertion * q_len;
            }
        }
    } else {
        // the reference string is NON-EMPTY
        if ( score_cols <= 1 ) {
            // the query string is EMPTY
            r_res = new char[ r_len + 1 ];
            q_res = new char[ r_len + 1 ];
            // no query, just ref, which remains untouched
            memcpy( r_res, r_str, r_len + 1 );
            // ref full of gaps
            memset( q_res, gap, sizeof( char ) * r_len );
            // null terminate
            r_res[ r_len ] = '\0';
            q_res[ r_len ] = '\0';
            // if do local, score is 0
            if ( ! do_local ) {
                if ( do_affine )
                    score = -open_deletion - ( r_len - 1 ) * extend_deletion;
                else
                    score = -open_deletion * r_len;
            }
        } else {
            long edit_ptr = 0;
            // don't forget the optional termination character
            
            signed char * const edit_ops = (signed char * const)alloca (sizeof (signed char) * ( r_len + q_len ));
            // allocate enough storage to store the best scoring path through the DP matrix

            cawlign_fp * const score_matrix = score_matrix_cache ?  score_matrix_cache : new cawlign_fp[ score_rows * score_cols ],
                   * const insertion_matrix = do_affine ? (insertion_matrix_cache ? insertion_matrix_cache : new cawlign_fp[ score_rows * score_cols ]) : NULL,
                   * const deletion_matrix  = do_affine ? (deletion_matrix_cache ? deletion_matrix_cache : new cawlign_fp[ score_rows * score_cols ]) : NULL;
            
            // allocate the DP matrix and, if using affine gaps, also allocate insertion and deletion DP matrices

            // encode each string using the character map (char_map)
            long * r_enc = NULL,
                 * q_enc = NULL;

 
            if ( do_codon ) {
                // encode both strings in [0, charCount],
                // i.e. something that can be directly looked up in the scoring matrix
                r_enc = (long*)alloca (sizeof (long) *  r_len );
                q_enc = (long*)alloca (sizeof (long) *  q_len );
                
                for (long i = 0; i < r_len; ++i ) {
                    r_enc[ i ] = char_map[ (unsigned char) r_str[ i ] ];
                }
                for (long i = 0; i < q_len; ++i ) {
                    q_enc[ i ] = char_map[ (unsigned char) q_str[ i ] ];
                }
            }
            
            //memset (score_matrix, 0, sizeof (cawlign_fp) * score_rows * score_cols);
            

            //if ( do_affine ) {
                //memset (insertion_matrix, 0, sizeof (cawlign_fp) * score_rows * score_cols);
                //memset (deletion_matrix , 0, sizeof (cawlign_fp) * score_rows * score_cols);
            //}

            score_matrix [ 0 ] = 0.;
            // pre-initialize the values in the various matrices
            if ( ! do_local ) {
                // full global alignment, i.e. indels at the beginning and end ARE penalized
                cawlign_fp cost;
                 // initialize gap costs in first column and first row
                // they are 0 for local alignments, so ignore

                if ( do_affine ) {
 
                    // first handle insertions
                    cost = -open_insertion;
                    insertion_matrix[ 0 ] = cost;

                    for (long i = 1; i < score_cols; ++i, cost -= extend_insertion ) {
                        score_matrix[ i ] = cost;
                        insertion_matrix[ i ] = cost;
                        deletion_matrix[ i ] = cost;
                    }

                    // then deletions
                    cost = -open_deletion;
                    deletion_matrix[ 0 ] = cost;

                    /** 20240219 : SLKP optimization note; may be faster to do 3 loops because of memory locality */
                    for (long i = score_cols; i < score_rows * score_cols; i += score_cols, cost -= extend_deletion ) {
                        score_matrix[ i ] = cost;
                        insertion_matrix[ i ] = cost;
                        deletion_matrix[ i ] = cost;
                    }
                } else {
                    // no affine gaps
                    if ( ! do_codon ) {
                        cost = -open_insertion;
                        for (long i = 1; i < score_cols; ++i, cost -= open_insertion )
                            score_matrix[ i ] = cost;

                        cost = -open_deletion;
                        for (long i = score_cols; i < score_rows * score_cols; i += score_cols, cost -= open_deletion )
                            score_matrix[ i ] = cost;

                        // handle the do_local, do_codon case
                    } else {
                        cost = -open_insertion;
                        /** 20240219 : SLKP optimization note; surely the next two loops don't need to do integer remainer at each iteration */
                        for (long i = 1; i < score_cols; ++i, cost -= open_insertion )
                            score_matrix[ i ] = cost - ( i % 3 != 1 ? miscall_cost : 0 );

                        cost = -open_deletion;
                        for (long i = score_cols, j = 0; i < score_rows * score_cols; i += score_cols, cost -= open_insertion, ++j )
                            score_matrix[ i ] = cost - ( j % 3 != 0 ? miscall_cost : 0 );
                    }
                }
             } else {
                 // here we're doing a local alignment,
                 // the costs of opening a deletion or an insertion
                 // remain the same no matter how far down the ref or query
                 // we've traveled, respectively

                if ( do_affine ) {
                    deletion_matrix [ 0 ] = 0.;
                    insertion_matrix [ 0 ] = 0.;
                    if ( do_codon ) {
                        /** 20240219 : SLKP optimization note; surely the next two loops don't need to do integer remainer at each iteration */

                        // XXX: should we be including the frameshift penalty here? I think not
                        // fill in the first row of the affine deletion matrix
                        // with the deletion cost plus the miscall penalty
                        
                        
                        for (long i = 1; i < score_cols; ++i ) {
                            deletion_matrix[ i ] = -open_deletion - ( i % 3 != 1 ? miscall_cost : 0 );
                            insertion_matrix [ i ] = 0.;
                            score_matrix [ i ] = 0.;
                        }

                        // fill in the first column of the affine insertion matrix
                        // with the insertion cost plus the miscall penalty
                        for (long i = score_cols, j = 0; i < score_rows * score_cols; i += score_cols, ++j ) {
                            
                            cawlign_fp is = -open_insertion - ( j % 3 != 0 ? miscall_cost : 0 );
                            
                            insertion_matrix[ i ] = is;
                            insertion_matrix[ i + 1 ] = 0.;//is;
                            insertion_matrix[ i + 2] = 0.;// is;

                            deletion_matrix [ i ] = 0.;
                            score_matrix [ i ] = 0.;
                        }
                        
                        for (long i = score_cols + 1; i < score_rows * score_cols; i += score_cols ) {
                        }
                        

                    } else {
                        // fill in the first row of the affine deletion matrix
                        // with the deletion cost
                        for (long i = 1; i < score_cols; ++i ) {
                            deletion_matrix[ i ] = -open_deletion;
                            insertion_matrix [ i ] = 0.;
                            score_matrix [ i ] = 0.;
                        }

                        // fill in the first column of the affine insertion matrix
                        // with the insertion cost
                        for (long i = score_cols; i < score_rows * score_cols; i += score_cols ) {
                            insertion_matrix[ i ] = -open_insertion;
                            deletion_matrix [ i ] = 0.;
                            score_matrix [ i ] = 0.;
                        }
                    }
                } else {
                    // init zeros in the score matrix for the non-affine case
                    for (long i = 1; i < score_cols; ++i ) {
                        score_matrix [ i ] = 0.;
                    }

                    // fill in the first column of the affine insertion matrix
                    // with the insertion cost
                    for (long i = score_cols; i < score_rows * score_cols; i += score_cols ) {
                        score_matrix [ i ] = 0.;
                    }
                }
            }

            if ( do_codon ) {
                /** populate the dynamic programming matrix here */
                for (long i = 1; i < score_rows; ++i )
                    for (long j = 1; j < score_cols; ++j )
                        CodonAlignStringsStep( score_matrix
                                             , r_enc
                                             , q_enc
                                             , i
                                             , j
                                             , score_cols
                                             , char_count
                                             , miscall_cost
                                             , open_insertion
                                             , open_deletion
                                             , extend_insertion
                                             , extend_deletion
                                             , cost_matrix
                                             , cost_stride
                                             , insertion_matrix
                                             , deletion_matrix
                                             , codon3x5
                                             , codon3x4
                                             , codon3x2
                                             , codon3x1
                                             , do_true_local
                                             );
                // not doing codon alignment
            } else {
                /** populate the dynamic programming matrix here */
                for (long i = 1; i < score_rows; ++i ) {
                    const long r_char = char_map[ (unsigned char) r_str[ i - 1 ] ];
                    for (long j = 1; j < score_cols; ++j ) {
                        
                        const long curr = ( i ) * score_cols + j,
                                   prev = ( i - 1 ) * score_cols + j;

                        // ref but not query is deletion
                        // query but not ref is insertion
                        cawlign_fp deletion  = score_matrix[ prev ] - open_deletion,
                               insertion = score_matrix[ curr - 1 ] - open_insertion,
                               match     = score_matrix[ prev - 1 ];

                        // if there is a match bonus or penalty, add it in
                        if ( r_char >= 0 ) {
                            const long q_char = char_map[ (unsigned char) q_str[ j - 1 ] ];
                            if ( q_char >= 0 ) {
                                match += cost_matrix[ r_char * cost_stride + q_char ];
                            }
                        }

                        // if we're doing affine gaps,
                        // look up potential moves in the affine gap matrices
                        // 20240219: make sure this works as intended
                        
                        if ( do_affine ) {
                            deletion  = MAX_OP( deletion,
                                             deletion_matrix[ prev ] - ( i > 1 ? extend_deletion : open_deletion ) ),
                            insertion = MAX_OP( insertion,
                                             insertion_matrix[ curr - 1 ] - ( j > 1 ? extend_insertion : open_insertion ) ),
                            // store the values back in the gap matrices
                            deletion_matrix[ curr ] = deletion;
                            insertion_matrix[ curr ] = insertion;
                        }

                        score_matrix[ curr ] = MAX_OP( match, MAX_OP( deletion, insertion ) );
                    }
                }
            }

            // set these indices to point at the ends
            // of the ref and query, respectively
            long index_R = r_len;
            long index_Q = q_len;
            
            bool took_local_shortcut = false;

            // grab maximum score from the last entry in the table
            score = score_matrix[ score_rows * score_cols - 1 ];

            // if we're doing a local alignment,
            if ( do_true_local) {
                // find the best score in the matrix
                // except for the first row/first column 
                // and start backtracking from there
                const cawlign_fp * score_row = score_matrix + score_cols;
                for (long m = 1; m < score_rows; m ++)  {
                    for (long k = 1; k < score_cols; k ++) {
                        if ( score_row[ k ] > score ) {
                            score = score_row[ k ];
                            index_R = ref_stride * m;
                            index_Q = k;
                        }
                    }
                    score_row += score_cols;
                }
                               
            } else 
                // find the best score in the last row and column of the scoring matrix
                // and start backtracking from there ( if it's better than the score
                // we've already found, that is )
                if ( do_local ) {
                    // grab the best score from the last column of the score matrix,
                    // skipping the very last entry ( we already checked it )
                    
                    for (long k = score_cols - 1; k < score_rows * score_cols - 1; k += score_cols )
                        if ( score_matrix[ k ] > score ) {
                            score = score_matrix[ k ];
                            // if do_codon, k / score_cols indexes into the codon space
                            // of the reference, which is resolved by multiplication
                            // by ref_stride ( which is 3 ), otherwise this
                            // directly indexes into the reference
                            index_R = ref_stride * ( k / score_cols );
                        }

                    // grab the best score from the last row of the score matrix,
                    // skipping the very last entry ( we already checked it )
                    for (long k = ( score_rows - 1 ) * score_cols; k < score_rows * score_cols - 1; ++k )
                        if ( score_matrix[ k ] > score ) {
                            score = score_matrix[ k ];
                            // if we've found a better score here,
                            // don't forget to reset the ref index
                            index_R = r_len;
                            // remove the initial value!
                            index_Q = k - ( score_rows - 1 ) * score_cols;
                        }

                    // fill in the edit_ops with the difference
                    // between r_len and i
                    for (long k = index_R; k < r_len; ++k ) {
                        edit_ops[ edit_ptr++ ] = -1;
                    }

                    // fill in the edit_ops with the difference
                    // between q_len and j
                    for (long k = index_Q; k < q_len; ++k ) {
                        edit_ops[ edit_ptr++ ] = 1;
                    }
                }

            // backtrack now

            /*
            // prints the score matrix
            for ( long m = 0; m < score_rows; ++m ) {
               for ( long n = 0; n < score_cols; ++n ) {
                   if ( n > 0 )
                       fprintf( stderr, "," );
                   fprintf( stderr, "% 3.3g", score_matrix[ m * score_cols + n ] );
               }
               fprintf( stderr, "\n" );
            }
            fprintf( stderr, "\n" );
            */

            if ( do_codon ) {
                // if either index hits 0, we're done
                // or if both indices fall below 3, we're done
                while ( index_R && index_Q && ( index_R >= 3 || index_Q >= 3 ) && !took_local_shortcut ) {
                    // perform a step
                    long code = CodonAlignStringsStep( score_matrix
                                                           , r_enc
                                                           , q_enc
                                                           // divide by 3 to index into codon space
                                                           , ( index_R / 3 )
                                                           , index_Q
                                                           , score_cols
                                                           , char_count
                                                           , miscall_cost
                                                           , open_insertion
                                                           , open_deletion
                                                           , extend_insertion
                                                           , extend_deletion
                                                           , cost_matrix
                                                           , cost_stride
                                                           , insertion_matrix
                                                           , deletion_matrix
                                                           , codon3x5
                                                           , codon3x4
                                                           , codon3x2
                                                           , codon3x1
                                                           , do_true_local
                                                           );

                    // alter edit_ops and decrement i and j
                    // according to the step k we took
                    if (do_true_local && code < 0) {
                         code = -code - 1;
                         took_local_shortcut = true;
                    }

                    
                    BacktrackAlignCodon( edit_ops, edit_ptr, index_R, index_Q, code );
                    

                    // if anything drops below 0, something bad happened
                    if ( index_R < 0 || index_Q < 0 ) {
                        //delete [] edit_ops;
                        return -INFINITY;
                    }

                    // handle the affine cases
                    if ( do_affine ) {
                        // divide by 3 to index into codon space
                        long k = ( index_R / 3 ) * score_cols + index_Q;
                        // reference matched but not query, a deletion
                        if ( code == HY_111_000 ) {
                            // while deletion is preferential to match
                            while ( index_R >= 3
                                 && score_matrix[ k ] - open_deletion
                                 <= deletion_matrix[ k ] - extend_deletion ) {
                                // take a codon out of the reference
                                index_R -= 3;
                                edit_ops[ edit_ptr++ ] = -1;
                                edit_ops[ edit_ptr++ ] = -1;
                                edit_ops[ edit_ptr++ ] = -1;
                                // move up a row in the score_matrix
                                // which is a codon in the reference
                                k -= score_cols;
                            }
                            // query matched but not reference, insertion
                        } else if ( code == HY_000_111 ) {
                            // while insertion is preferential to match
                            while ( index_Q >= 3
                                 && score_matrix[ k ] - open_insertion
                                 <= insertion_matrix[ k ] - extend_insertion ) {
                                // take a codon out of the query
                                index_Q -= 3;
                                edit_ops[ edit_ptr++ ] = 1;
                                edit_ops[ edit_ptr++ ] = 1;
                                edit_ops[ edit_ptr++ ] = 1;
                                // move up 3 in the score_matrix
                                // which is a codon in the query
                                k -= 3;
                            }
                        }
                    }
                }
            } else {
                if ( do_affine ) {
                    while ( index_R && index_Q ) {
                        long curr = ( index_R ) * score_cols + index_Q,
                             prev = ( index_R - 1 ) * score_cols + index_Q,
                             best_choice = 0;

                        // check the current affine scores and the match score
                        cawlign_fp scores[ 3 ] = {
                            deletion_matrix[ curr ],
                            insertion_matrix[ curr ],
                            score_matrix[ prev - 1 ]
                        }, max_score = scores[ best_choice ];

                        MatchScore( r_str, q_str, index_R, index_Q, char_map, cost_matrix, cost_stride, scores[2] );

                        // look at choice other than 0
                        if (scores[1] > max_score) {
                            max_score = scores[1];
                            best_choice = 1;
                        }
                        if (scores[2] > max_score) {
                            best_choice = 2;
                        }
                        
                        switch ( best_choice ) {
                        case 0:
                            // we have at least 1 deletion
                            --index_R;
                            edit_ops[ edit_ptr++ ] = -1;
                            // deletion is travel in the reference but not query,
                            // look at scores back in the reference,
                            // and while they are better for the deletion case,
                            // move backwards in the reference
                            while ( index_R
                                 && score_matrix[ curr - score_cols ] - open_deletion
                                 <= deletion_matrix[ curr - score_cols ] - extend_deletion
                                  ) {
                                --index_R;
                                edit_ops[ edit_ptr++ ] = -1;
                                curr -= score_cols;
                            }
                            break;

                        case 1:
                            // we have at least 1 insertion
                            --index_Q;
                            edit_ops[ edit_ptr++ ] = 1;
                            // insertion is travel in the query but not the reference,
                            // look at scores back in the query,
                            // and while they are better than for the insertion case,
                            // move backwards in the query
                            while ( index_Q
                                 && score_matrix[ curr - 1 ] - open_insertion
                                 <= insertion_matrix[ curr - 1 ] - extend_insertion
                                  ) {
                                --index_Q;
                                edit_ops[ edit_ptr++ ] = 1;
                                --curr;
                            }
                            break;

                        case 2:
                            // it's a match! move back in both
                            --index_R;
                            --index_Q;
                            edit_ops[ edit_ptr++ ] = 0;
                            break;
                        }
                    }
                    // no affine gaps, no codons
                } else {
                    while ( index_R && index_Q ) {
                        const long curr = ( index_R ) * score_cols + index_Q,
                                   prev = ( index_R - 1 ) * score_cols + index_Q;

                        cawlign_fp deletion  = score_matrix[ prev ] - open_deletion,
                               insertion = score_matrix[ curr - 1 ] - open_insertion,
                               match     = score_matrix[ prev - 1 ];

                        MatchScore( r_str, q_str, index_R, index_Q, char_map, cost_matrix, cost_stride, match );
                        BacktrackAlign( edit_ops, edit_ptr, index_R, index_Q, deletion, insertion, match );
                    }
                }
            }


            
            if (!took_local_shortcut) {
                // for anything that remains,
                // don't forget it!!!
                // reference
                
                
                while ( --index_R >= 0 )
                    edit_ops[ edit_ptr++ ] = -1;

                // then query
                while ( --index_Q >= 0 )
                    edit_ops[ edit_ptr++ ] = 1;
            }

            if ( edit_ptr > 0 ) {
                // reset indices to 0
                if (!took_local_shortcut){ 
                    index_Q = index_R = 0;
                }

                // rebuild the strings from the edit_ops
                // with room for the null terminator
                r_res = new char[ edit_ptr + 1 ];
                q_res = new char[ edit_ptr + 1 ];
                --edit_ptr;
                long k = 0;
                for (; edit_ptr >= 0; --edit_ptr, ++k ) {
                    switch ( edit_ops[ edit_ptr ] ) {
                            // match! include characters from both strings
                        case 0:
                            r_res[ k ] = r_str[ index_R++ ];
                            q_res[ k ] = q_str[ index_Q++ ];
                            break;
                            // insertion!
                        case 1:
                            if (report_ref_insertions) {
                                r_res[ k ] = gap;
                                q_res[ k ] = q_str[ index_Q++ ];
                            } else {
                                k --; index_Q++;
                            }
                            break;
                        case 2:
                            if (report_ref_insertions) {
                                r_res[ k ] = gap;
                                q_res[ k ] = tolower( q_str[ index_Q++ ] );
                            } else {
                                k --; index_Q++;
                            }
                            break;
                        case -1:
                            r_res[ k ] = r_str[ index_R++ ];
                            q_res[ k ] = gap;
                            break;
                        case -2:
                            r_res[ k ] = tolower( r_str[ index_R++ ] );
                            q_res[ k ] = gap;
                            break;
                    }
                }
                // make sure to null-terminate
                r_res[ k ] = '\0';
                q_res[ k ] = '\0';
            }

#ifdef ALIGN_DEBUG
            _String alignDebug( "alignScoreMatrix" );
            _Variable * ad = CheckReceptacle( &alignDebug, empty, false );
            ad->SetValue( score_matrix, true );

            // grab the affine matrices too,
            // if that's what we're doing
            if ( do_affine ) {
                _String alignDebug( "alignScoreMatrixG1" );
                _Variable * ad = CheckReceptacle( &alignDebug, empty, false );
                ad->SetValue( insertion_matrix, true );
                alignDebug = ( "alignScoreMatrixG2" );
                ad = CheckReceptacle( &alignDebug, empty, false );
                ad->SetValue( deletion_matrix, true );
            }
#endif
            //delete [] edit_ops;
            if (score_matrix != score_matrix_cache) {
                delete [] score_matrix;
            }

            if ( do_affine ) {
                if (insertion_matrix != insertion_matrix_cache) {
                    delete [] insertion_matrix;
                }
                if (deletion_matrix != deletion_matrix_cache) {
                    delete [] deletion_matrix;
                }
            }

            /*
            if (do_codon) {
                delete [] r_enc;
                delete [] q_enc;
            }
            */
        }
    }

    return score;
}

//____________________________________________________________________________________

#define _ALIGNMENT_NOLOCAL      0x00
#define _ALIGNMENT_LOCAL_START  0x01
#define _ALIGNMENT_LOCAL_END    0x02

//____________________________________________________________________________________

/**
 * Computes the alignment score between two sequences using affine or linear gap penalties.
 *
 * Calculates the alignment score for two sequences, `s1` and `s2`, based on the 
 * specified scoring matrix and gap penalties. Optionally supports local alignments (ignoring 
 * prefix and suffix gaps) and affine gap penalties (distinct gap opening and extension costs).
 * Allows alignment of reversed sequences and returns the computed score along 
 * with how the alignment was achieved.
 *
 * @param s1 The first sequence (as a character array).
 * @param s2 The second sequence (as a character array).
 * @param s1L Length of the first sequence (`s1`).
 * @param s2L Length of the second sequence (`s2`).
 * @param from1 The starting index in the first sequence for the alignment.
 * @param from2 The starting index in the second sequence for the alignment.
 * @param to1 The ending index (exclusive) in the first sequence for the alignment.
 * @param to2 The ending index (exclusive) in the second sequence for the alignment.
 * @param rev1 If `true`, align the first sequence in reverse.
 * @param rev2 If `true`, align the second sequence in reverse.
 * @param cmap Character-to-position mapping in the scoring matrix.
 * @param ccost Scoring matrix for character-to-character match/mismatch scores.
 * @param mapL Dimension of the scoring matrix `ccost`.
 * @param gopen Gap opening penalty for the first sequence.
 * @param gextend Gap extension penalty for the first sequence.
 * @param gopen2 Gap opening penalty for the second sequence.
 * @param gextend2 Gap extension penalty for the second sequence.
 * @param doLocal If `true`, perform local alignment (ignoring prefix and suffix gaps).
 * @param doAffine If `true`, use affine gap penalties (distinct opening and extension penalties).
 * @param scoreMatrix Matrix to store the last row of the alignment scores.
 * @param gapScore1 Matrix to store gap penalties for the first sequence (used if `doAffine` is `true`).
 * @param gapScore2 Matrix to store gap penalties for the second sequence (used if `doAffine` is `true`).
 * @param secondGap If `true`, align with a gap in the second sequence.
 * @param howAchieved Array that records how each position in the alignment was achieved (match, insert, or delete).
 * @return The alignment score for the two sequences.
 */
cawlign_fp   CostOnly   (  const char * s1,
                        const char * s2,
                        const long s1L,
                        const long s2L,
                        long from1,
                        long from2,
                        long to1,
                        long to2,
                        bool rev1,
                        bool rev2,
                        long*  cmap,
                        cawlign_fp const * ccost,
                        const  long    mapL,
                        cawlign_fp gopen,
                        cawlign_fp gextend,
                        cawlign_fp gopen2,
                        cawlign_fp gextend2,
                        bool doLocal,
                        bool doAffine,
                        cawlign_fp * scoreMatrix,
                        cawlign_fp * gapScore1,
                        cawlign_fp * gapScore2,
                        char secondGap,
                        char * howAchieved
                     )
{
    cawlign_fp   score    = 0.;

    long         s1Length = to1-from1,
                 s2Length = to2-from2;


    bool         doLocal1S = false,
                 doLocal1E = false,
                 doLocal2S = false,
                 doLocal2E = false;

    if (doLocal) {
        if (rev1) {
            doLocal1S = (to1==s1L);
            doLocal1E = from1 == 0L;
            //doLocal1 = (to1==s1->sLength)*_ALIGNMENT_LOCAL_START + (from1 == 0)*_ALIGNMENT_LOCAL_END;
        } else {
            doLocal1E = (to1==s1L);
            doLocal1S = from1 == 0L;
            //doLocal1 = (from1==0)*_ALIGNMENT_LOCAL_START + (to1==s1->sLength)*_ALIGNMENT_LOCAL_END;
        }
        if (rev2) {
            doLocal2E = from2 == 0L;
            doLocal2S = (to2==s2L);
            //doLocal2 = (to2==s2->sLength)*_ALIGNMENT_LOCAL_START + (from2 == 0)*_ALIGNMENT_LOCAL_END;
        } else {
            doLocal2S = from2 == 0L;
            doLocal2E = (to2==s2L);
            //doLocal2 = (from2==0)*_ALIGNMENT_LOCAL_START + (to2==s2->sLength)*_ALIGNMENT_LOCAL_END;
        }
    }

    if (s1Length) {
        // first string not empty
        if (s2Length) {
            // second string not empty
            cawlign_fp          aux2;
            long            colCount = s2Length+1;

            scoreMatrix[0] = 0.;
            if (doAffine) {
                gapScore1[0] = gapScore2[0] = 0.;
            }


            if (doLocal1S == 0) {
                cawlign_fp cost = -gopen;
                if (doAffine) {
                    for (long k=1; k < colCount; k++, cost-=gextend) {
                        scoreMatrix [k]  = cost;
                        gapScore1   [k]  = cost;
                        gapScore2   [k]  = cost;
                    }
                } else
                    for (long m=1; m < colCount; m++, cost-=gopen) {
                        scoreMatrix [m] = cost;
                    }
            } else {
                for (long k=1; k < colCount; k++) {
                    scoreMatrix [k] = 0.;
                }

                if (doAffine) {
                    for (long k=1; k < colCount; k++) {
                        gapScore1[k] = 0;
                        gapScore2[k] = -(secondGap==1?gextend2:gopen2);
                        // prefix gaps in the second sequence
                    }
                    gapScore1[0] = -gopen;
                }
            }

            // long mapL = ccost->GetVDim();
            // how many valid characters

            if (doAffine) {
                aux2 = 0.;

                if (doLocal2S == 0) {
                    gapScore1[0] = gapScore2[0] = -(secondGap==1?gextend2:gopen2);
                }

                from2 --;
                from1 --;
                for (long r=1; r<=s1Length; r++) { // iterate by rows
                    long      c1 = cmap[s1[rev1?(to1-r):(from1+r)]];

                    if (doLocal2S) {
                        aux2        = 0.;
                    } else {
                        if (r>1) {
                            aux2         = -((r-2)*gextend2 + (secondGap==1?gextend2:gopen2));
                        }
                        scoreMatrix[0] = -((secondGap==1?gextend2:gopen2) + (r-1)*gextend2);
                    }

                    for (long c=1; c<=s2Length; c++) { // iterate by columns
                        cawlign_fp gscore1  ,           // gap in 2nd
                        gscore2  ,           // gap in 1st
                        gscore3  = aux2,     // no gap
                        t;

                        // if secondGap == 2, then we MUST _start_ with a gap in the 2nd sequence


                        if (doLocal1E && r == s1Length) {
                            //gscore2 = MAX(scoreMatrix.theData[c-1],gapScore1->theData[c-1]);
                            gscore2 = scoreMatrix[c-1];
                            if (gapScore1[c-1] > gscore2) {
                                gscore2 = gapScore1[c-1];
                            }
                        } else {
                            gscore2 = scoreMatrix[c-1]-gopen;
                            t       = gapScore1[c-1]-((c>1)?gextend:gopen);
                            if (t > gscore2) {
                                gscore2 = t;
                            }
                        }

                        if (doLocal2E && c == s2Length) {
                            //gscore1 = MAX(scoreMatrix.theData[c],gapScore2->theData[c]);
                            gscore1 = scoreMatrix[c];
                            if (gscore1 < gapScore2[c]) {
                                gscore1 = gapScore2[c];
                            }
                        } else {
                            //gscore1 = MAX(scoreMatrix.theData[c]-gopen2,gapScore2->theData[c]-((r>1)?gextend2:gopen2));
                            gscore1 = scoreMatrix [c]-gopen2;
                            t       = gapScore2[c]-((r>1)?gextend2:gopen2);
                            if (t > gscore1) {
                                gscore1 = t;
                            }
                        }
                        // either open a new gap from a character; or continue an existing one
                        // if this is the second row, then we start a gap in the second sequence -|

                        if (c1>=0) {
                            long       c2 = cmap[s2[rev2?(to2-c):(from2+c)]];

                            if (c2>=0) {
                                gscore3 += ccost[c1*mapL+c2];
                            }
                        }

                        aux2                    = scoreMatrix[c];
                        char                      option = 0;
                        t                       = gscore2;


                        if (r > 1 || secondGap == 0) {
                            if (gscore1 > gscore2) {
                                t = gscore1;
                                option                 = 1;
                            }
                            if (gscore3 > t) {
                                t                      = gscore3;
                                option                 = 2;
                            }
                        }
                        scoreMatrix[c] = t;
                        if (howAchieved) {
                            howAchieved[c] = option;
                        }

                        //if (rev2 && secondGap==2 && c == s2Length)
                        //  gscore1 = MAX(scoreMatrix.theData[c]-gextend2,gapScore2->theData[c]-((r>1)?gextend2:gopen2));

                        gapScore2 [c]  = gscore1;
                        gapScore1 [c]  = gscore2;

                    }

                    if (doLocal2S && r < s1Length) {
                        gapScore1[0] -= gextend2;
                        gapScore2[0] -= gextend2;
                    }
                }
            } else {
                // populate the cost matrix row by row
                aux2 = 0.;
                for (long r=1; r<=s1Length; r++) {
                    if (doLocal2S) {
                        aux2        = 0.;
                    } else {
                        scoreMatrix[0] = -(gopen2 * r);
                        if (r>1) {
                            aux2         = -((r-1)*gopen2);
                        }
                    }

                    //printf ("%d: %g\t", r, scoreMatrix.theData[0]);
                    long      c1 = cmap[s1[rev1?(to1-r):(from1+r-1)]];

                    for (long c=1; c<=s2Length; c++) {
                        cawlign_fp score1 = scoreMatrix[c], // gap in 2nd
                        score2 = scoreMatrix[c-1],  // gap in 1st
                        score3 = aux2;

                        if (c < s2Length || doLocal2E == 0) {
                            score1 -= gopen2;
                        }
                        if (r < s1Length || doLocal1E == 0) {
                            score2 -= gopen;
                        }

                        if (c1>=0) {
                            long       c2 = cmap[s2[rev2?(to2-c):(from2+c-1)]];

                            if (c2>=0) {
                                score3 += ccost [c1*mapL+c2];
                            }
                        }

                        aux2                    = scoreMatrix[c];
                        char                    option = 0;
                        scoreMatrix[c]  = score1;
                        if (score2 > score1) {
                            scoreMatrix[c] = score2;
                            option                 = 1;
                        }
                        if (score3 > scoreMatrix[c]) {
                            scoreMatrix[c] = score3;
                            option                 = 2;
                        }
                        if (howAchieved) {
                            howAchieved[c] = option;
                        }
                    }
                    //printf ("\n");
                }

            }
            score = scoreMatrix[s2Length];
        } else { // 2nd string empty
            if ((doLocal2S || doLocal2E) == false) {
                if (doAffine) {
                    score = gopen2+gextend2*s1Length;
                } else {
                    score = gopen2 * s1Length;
                }
            }
        }
    } else // first string empty
        if (s2Length) { // second string not empty
            if ((doLocal1S || doLocal1E) == false) {
                score = -gopen;

                scoreMatrix[0] = 0.0;
                if (doAffine) {
                    gapScore1[0] = gapScore2[0] = 0.0;
                    for (long k = 1; k <= s2Length; k++, score-=gextend) {
                        scoreMatrix[k] = gapScore1[k] = gapScore2[k] = score;
                    }

                    score += gextend;
                } else {
                    for (long k = 1; k <= s2Length; k++, score-=gopen) {
                        scoreMatrix[k] = score;
                    }
                    score += gopen;
                }
            } else {
                for (long k = 0; k <= s2Length; k++) {
                    scoreMatrix[k] = 0.;
                }
                if (doAffine) {
                    for (long k = 0; k <= s2Length; k++) {
                        gapScore1[k] = 0.;
                        gapScore2[k] = 0.;
                    }
                }
            }

        }

    return score;
}

//____________________________________________________________________________________

/**
 * Performs sequence alignment in linear space, reducing memory usage for long sequences.
 *
 * This function implements the Needleman-Wunsch algorithm with affine or linear gap penalties
 * while keeping memory usage low by only storing a small portion of the alignment matrix. It uses
 * a divide-and-conquer approach to calculate the optimal alignment between `s1` and `s2`.
 * It supports both local alignments (ignoring prefix and suffix gaps) and affine gap penalties.
 *
 * @param s1 The first sequence (as a character array).
 * @param s2 The second sequence (as a character array).
 * @param s1L Length of the first sequence (`s1`).
 * @param s2L Length of the second sequence (`s2`).
 * @param cmap Character-to-position mapping in the scoring matrix.
 * @param ccost Scoring matrix for character-to-character match/mismatch scores.
 * @param costD Dimension of the scoring matrix `ccost`.
 * @param gopen Gap opening penalty for the first sequence.
 * @param gextend Gap extension penalty for the first sequence.
 * @param gopen2 Gap opening penalty for the second sequence.
 * @param gextend2 Gap extension penalty for the second sequence.
 * @param doLocal If `true`, perform local alignment (ignoring prefix and suffix gaps).
 * @param doAffine If `true`, use affine gap penalties (distinct opening and extension penalties).
 * @param ops Array to store the optimal alignment operations (edit steps).
 * @param scoreCheck Variable to check the final score of the alignment.
 * @param from1 Starting index in the first sequence.
 * @param to1 Ending index (exclusive) in the first sequence.
 * @param from2 Starting index in the second sequence.
 * @param to2 Ending index (exclusive) in the second sequence.
 * @param buffer Matrix storage for the alignment computation.
 * @param parentGapLink Keeps track of the previous gap state (to decide whether to continue or close a gap).
 * @param ha Array that stores how the alignment was achieved (match, insert, or delete).
 * @return The alignment score for the two sequences.
 */
cawlign_fp      LinearSpaceAlign (   const char    *s1,
                                  const char    *s2,
                                  const long    s1L,
                                  const long    s2L,

                                  long        * cmap,
                                  cawlign_fp const* ccost,
                                  const long    costD,
                                  cawlign_fp        gopen,
                                  cawlign_fp        gextend,
                                  cawlign_fp        gopen2,
                                  cawlign_fp        gextend2,
                                  bool          doLocal,
                                  bool          doAffine,
                                  long*         ops,
                                  cawlign_fp        scoreCheck,
                                  long          from1,
                                  long          to1,
                                  long          from2,
                                  long          to2,
                                  cawlign_fp        **buffer,
                                  char          parentGapLink,
                                  char          *ha
                                  )
{
    if (to2 == from2 || to1 == from1) {
        return 0;
    }

    long                    midpoint = (from1 + to1)/2,
    span     = to2-from2,
    span1     = to1-from1;

    if                      (span1 > 1) {
        CostOnly                (s1,
                                 s2,
                                 s1L,
                                 s2L,
                                 from1,
                                 from2,
                                 midpoint,
                                 to2,
                                 false,
                                 false,
                                 cmap,
                                 ccost,
                                 costD,
                                 gopen,
                                 gextend,
                                 gopen2,
                                 gextend2,
                                 doLocal,
                                 doAffine,
                                 buffer[0],
                                 buffer[1], 
                                 buffer[2],
                                 parentGapLink>=2, 
                                 ha);
        CostOnly                (s1,s2,s1L,s2L, midpoint,from2,to1,to2,true,true,  cmap,ccost,costD, gopen,gextend,gopen2,gextend2,doLocal,doAffine,buffer[3], buffer[4], buffer[5],   2*(parentGapLink%2), ha+s2L+1UL);
    } else {
        CostOnly                (s1,s2,s1L,s2L,from1,from2,to1,to2,false,false,cmap,ccost,costD, gopen,gextend,gopen2,gextend2,doLocal,doAffine,buffer[0], buffer[1], buffer[2], (parentGapLink>=2), ha);
    }

    cawlign_fp     maxScore = -1e100;
    long       maxIndex = 0;
    bool       gapLink  = false;
    char       alignmentKind    = 0;

    cawlign_fp    gapOffsetScore   = gopen2-gextend2;
    
    if (!doAffine) {
        if (span1 > 1) {
            for (long k = 0; k <= span; k++) {
                cawlign_fp currentScore = buffer[0][k] + buffer[3][span-k];
                if (currentScore > maxScore) {
                    maxScore = currentScore;
                    maxIndex = k;
                }
            }
        } else { // handle the case of a single row span correctly
            for (long k = 0; k <= span; k++) {
                cawlign_fp currentScore     = buffer[0][k];

                if (! doLocal || to1 != s1L ) {
                    currentScore -= gopen*(span-k);
                }

                if (currentScore > maxScore) {
                    maxScore        = currentScore;
                    alignmentKind   = ha[k];
                    maxIndex = k;
                }
            }
        }
    } else {
        if (span1 > 1) {
            // two cases here: no-gap link
            // or gap-to-gap link

            for (long k = 0; k <= span; k++) {
                cawlign_fp currentScoreNoGap    = buffer[0][k] + buffer[3][span-k],
                currentScoreWithGap2  = buffer[2][k] + buffer[5][span-k] + gapOffsetScore;


                if (doAffine && (((from1 == 0 || from2==0) && k == 0) || ((to1 == s1L || to2 == s2L) && k == span))) {
                    currentScoreWithGap2 -= gapOffsetScore;
                }

                if (currentScoreNoGap > maxScore) {
                    maxScore = currentScoreNoGap;
                    maxIndex = k;
                    gapLink  = false;
                }
                if (currentScoreWithGap2 > maxScore) {
                    maxScore = currentScoreWithGap2;
                    maxIndex = k;
                    gapLink  = true;
                }
                /*printf ("[%d %d %d %d] (%d) %d %g %g: %g %g / %g %d\n", from1, to1, from2, to2, parentGapLink,  k,
                 buffer[0]->theData[k],  buffer[3]->theData[span-k], buffer[2]->theData[k],  buffer[5]->theData[span-k],
                 maxScore, maxIndex);*/

            }

        } else { // handle the case of a single row span correctly
            if (parentGapLink == 1) {
                maxIndex      = span;
                maxScore      = buffer[2][span];
                alignmentKind = 1;
            } else {
                for (long k = 0; k <= span; k++) {
                    cawlign_fp currentScoreNoGap    = buffer[0][k],
                    currentScoreWithGap2  = buffer[2][k];

                    if (!doLocal || to1 != s1L) // indel in sequence 1
                        if (span-k) {
                            currentScoreNoGap       -= gopen;
                            currentScoreWithGap2    -= gopen;
                            if (span-k>1) {
                                currentScoreNoGap    -= gextend*(span-k-1);
                                currentScoreWithGap2 -= gextend*(span-k-1);
                            }
                        }

                    /*printf ("[%d %d %d %d] %d %g %g: %g %g / %g %d\n", from1, to1, from2, to2, k,
                     buffer[0]->theData[k],  buffer[2]->theData[k], currentScoreNoGap, currentScoreWithGap2,
                     maxScore, maxIndex);*/

                    if (currentScoreNoGap > maxScore) {
                        maxScore = currentScoreNoGap;
                        maxIndex = k;
                        alignmentKind   = ha[k];
                    }
                    
                    if (currentScoreWithGap2 > maxScore) {
                        maxScore = currentScoreWithGap2;
                        maxIndex = k;
                        alignmentKind   = 0;
                    }
                }
            }
        }
    }

    if (span1 == 1) {
        if (alignmentKind == 2) {
            ops[from1+1] = from2+maxIndex-1;
        } else if (alignmentKind == 0 && maxIndex == 0/*&& to2 == s2->sLength && to1 == s1->sLength*/) {
            ops[from1+1]  = -3;
        }
    } else {

        cawlign_fp check1 = buffer[0][maxIndex],
               check2 = buffer[3][span-maxIndex];

        if (span1>1) {
            if (maxIndex > 0) {
                char gapCode = gapLink;
                if (parentGapLink >= 2) {
                    gapCode += 2;
                }
                LinearSpaceAlign (s1,s2,s1L,s2L, cmap,ccost,costD, gopen,gextend,gopen2,gextend2,doLocal,doAffine,ops,check1, from1, midpoint, from2, from2 + maxIndex, buffer, gapCode, ha);
            } else if (from2 == 0)
                for (long k = from1; k < midpoint; k++) {
                    ops[k+1] = -3;
                }

            if (maxIndex < span) {
                char gapCode = 2*gapLink;
                if (parentGapLink % 2 == 1) {
                    gapCode ++;
                }
                LinearSpaceAlign (s1,s2,s1L, s2L, cmap,ccost,costD, gopen,gextend,gopen2,gextend2,doLocal,doAffine,ops,check2, midpoint, to1, from2 + maxIndex, to2, buffer, gapCode, ha);
            }
        }
    }
    return maxScore;
}
