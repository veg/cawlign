/*
 
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (spond@ucsd.edu)
 Art FY Poon    (apoon42@uwo.ca)
 Steven Weaver (sweaver@ucsd.edu)
 
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

#ifndef __ALIGNMENT_HEADER_FILE__

#define __ALIGNMENT_HEADER_FILE__

typedef   float     cawlign_fp;

cawlign_fp AlignStrings( char const * r_str
                   , char const * q_str
                   , const long _r_len
                   , const long _q_len
                   , char * & r_res
                   , char * & q_res
                   , long * char_map
                   , const cawlign_fp * cost_matrix
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
                   , const bool do_true_local = false
                   , const bool report_ref_insertions = true
                   , cawlign_fp* score_matrix_cache = nullptr
                   , cawlign_fp* insertion_matrix_cache = nullptr
                   , cawlign_fp* deletion_matrix_cache = nullptr
                   , const long* resolution_map = nullptr

                   );

cawlign_fp LinearSpaceAlign(    const char * s1           // first string
                           , const char * s2           // second string
                           , const long s1L
                           , const long s2L
                           , long*  cmap     // char -> position in scoring matrix mapper
                           , const cawlign_fp * ccost        // NxN matrix of edit distances on characters
                           , const long costD
                           , cawlign_fp gopen       // the cost of opening a gap in sequence 1
                           , cawlign_fp gextend     // the cost of extending a gap in sequence 1 (ignored unless doAffine == true)
                           , cawlign_fp gopen2      // the cost of opening a gap in sequence 2
                           , cawlign_fp gextend2    // the cost of opening a gap in sequence 2   (ignored unless doAffine == true)
                           , bool doLocal           // ignore prefix and suffix gaps
                           , bool doAffine          // use affine gap penalties
                           , long * ops      // edit operations for the optimal alignment
                           , cawlign_fp scoreCheck  // check the score of the alignment
                           , long from1
                           , long to1
                           , long from2
                           , long to2
                           , cawlign_fp ** buffer      // matrix storage,
                           , char parentGapLink
                           , char * ha
                           );

#endif
