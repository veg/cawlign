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
#include <memory>
#include <vector>
#include <algorithm>
#include <mutex>

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
#endif

#include <math.h>
#include <stdio.h>

#include "alignment.h"

#define MAX_OP(a, b) ( (( a ) > ( b )) ? ( a ) : ( b ) )

struct RefKmerEntry {
    int32_t pos1;
    int32_t pos2;
};

static const int8_t kCharTo2Bit[256] = {
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1, 3, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1, 3, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1
};

static std::mutex g_ref_kmer_mutex;
static std::shared_ptr<const std::vector<RefKmerEntry>> g_cached_ref_table;
static unsigned long g_cached_ref_len = 0;
static uint64_t g_cached_ref_hash = 0;

static uint64_t ComputeStringHash64(const char* str, unsigned long len) {
    uint64_t hash = 14695981039346656037ULL;
    for (unsigned long i = 0; i < len; ++i) {
        hash ^= (uint64_t)(unsigned char)str[i];
        hash *= 1099511628211ULL;
    }
    return hash;
}

static std::shared_ptr<const std::vector<RefKmerEntry>> GetRefKmerTable(const char* r_str, unsigned long r_len) {
    uint64_t hash = ComputeStringHash64(r_str, r_len);
    std::lock_guard<std::mutex> lock(g_ref_kmer_mutex);
    if (g_cached_ref_table && g_cached_ref_len == r_len && g_cached_ref_hash == hash) {
        return g_cached_ref_table;
    }
    auto table = std::make_shared<std::vector<RefKmerEntry>>(1 << 20, RefKmerEntry{-1, -1});
    if (r_len >= 10) {
        for (long r = 0; r <= (long)r_len - 10; ++r) {
            uint32_t key = 0;
            bool ok = true;
            for (int k = 0; k < 10; ++k) {
                int8_t b = kCharTo2Bit[(unsigned char)r_str[r + k]];
                if (b < 0) { ok = false; break; }
                key = (key << 2) | b;
            }
            if (ok) {
                if ((*table)[key].pos1 < 0) {
                    (*table)[key].pos1 = (int32_t)r;
                } else if ((*table)[key].pos2 < 0) {
                    (*table)[key].pos2 = (int32_t)r;
                }
            }
        }
    }
    g_cached_ref_table = table;
    g_cached_ref_len = r_len;
    g_cached_ref_hash = hash;
    return table;
}

static bool is_band_disabled() {
    static const bool disabled = (getenv("CAWLIGN_NO_BAND") != nullptr);
    return disabled;
}

struct KmerHit {
    int32_t r;
    int32_t q;
    int32_t d;
};

static bool BuildCorridor(
    const char* q_str, unsigned long q_len,
    unsigned long r_len, unsigned long score_rows,
    long W,
    const std::vector<RefKmerEntry>& ref_kmer_table,
    std::vector<long>& j_min,
    std::vector<long>& j_max
) {
    if (q_len < 150 || r_len < 300) return false;

    std::vector<KmerHit> raw_hits;
    raw_hits.reserve(q_len / 3);

    const long diag_offset = (long)r_len;
    const long diag_range = (long)(r_len + q_len + 1);
    std::vector<int> diag_votes(diag_range, 0);

    for (long q = 0; q <= (long)q_len - 10; q += 3) {
        uint32_t key = 0;
        bool ok = true;
        for (int k = 0; k < 10; ++k) {
            int8_t b = kCharTo2Bit[(unsigned char)q_str[q + k]];
            if (b < 0) { ok = false; break; }
            key = (key << 2) | b;
        }
        if (!ok) continue;
        RefKmerEntry e = ref_kmer_table[key];
        if (e.pos1 >= 0) {
            int d1 = (int)q - e.pos1;
            raw_hits.push_back({e.pos1, (int32_t)q, d1});
            diag_votes[d1 + diag_offset]++;
        }
        if (e.pos2 >= 0) {
            int d2 = (int)q - e.pos2;
            raw_hits.push_back({e.pos2, (int32_t)q, d2});
            diag_votes[d2 + diag_offset]++;
        }
    }

    if (raw_hits.empty()) return false;

    int max_votes = 0;
    int best_d = 0;
    for (size_t d = 0; d < diag_votes.size(); ++d) {
        if (diag_votes[d] > max_votes) {
            max_votes = diag_votes[d];
            best_d = (int)d - (int)diag_offset;
        }
    }
    if (max_votes < 5) return false;

    const long W_nt = W;
    std::vector<KmerHit> filtered_hits;
    filtered_hits.reserve(raw_hits.size());
    for (const auto& h : raw_hits) {
        if (diag_votes[h.d + diag_offset] >= 5 && std::labs(h.d - best_d) <= W_nt) {
            filtered_hits.push_back(h);
        }
    }

    if (filtered_hits.size() < 10) return false;

    std::sort(filtered_hits.begin(), filtered_hits.end(), [](const KmerHit& a, const KmerHit& b) {
        return a.r < b.r || (a.r == b.r && a.q < b.q);
    });

    std::vector<std::pair<long, long>> chain;
    chain.reserve(filtered_hits.size());
    for (const auto& h : filtered_hits) {
        if (chain.empty() || (h.r > chain.back().first && h.q >= chain.back().second)) {
            chain.push_back({h.r, h.q});
        }
    }

    if (chain.size() < 10) return false;

    std::vector<std::pair<long, long>> anchors;
    anchors.reserve(chain.size());
    for (const auto& p : chain) {
        long row = p.first / 3;
        long col = p.second;
        if (anchors.empty() || row > anchors.back().first) {
            anchors.push_back({row, col});
        }
    }

    if (anchors.size() < 10) return false;

    const long r_span = (anchors.back().first - anchors.front().first) * 3;
    const long q_span = anchors.back().second - anchors.front().second;
    const long min_len = (long)std::min(r_len, q_len);

    if (r_span < 0.5 * min_len && q_span < 0.5 * min_len) return false;

    const long lead_unanchored_r = anchors.front().first * 3;
    const long trail_unanchored_r = (long)r_len - anchors.back().first * 3;
    const long lead_unanchored_q = anchors.front().second;
    const long trail_unanchored_q = (long)q_len - anchors.back().second;

    if (lead_unanchored_q - lead_unanchored_r > W) return false;
    if (trail_unanchored_q - trail_unanchored_r > W) return false;

    if (lead_unanchored_r > 150 && lead_unanchored_q > 150) return false;
    if (trail_unanchored_r > 150 && trail_unanchored_q > 150) return false;

    for (size_t a = 0; a + 1 < anchors.size(); ++a) {
        long gap_r = (anchors[a + 1].first - anchors[a].first) * 3;
        long gap_q = anchors[a + 1].second - anchors[a].second;
        long gap_diag = std::labs(gap_q - gap_r);
        if (gap_diag > W || (gap_r > 300 && gap_q > 300)) return false;
    }

    std::vector<long> j_star(score_rows);
    long i0 = anchors[0].first;
    long q0 = anchors[0].second;
    for (long i = 0; i < i0; ++i) {
        j_star[i] = std::max(0L, q0 - 3 * (i0 - i));
    }
    for (size_t a = 0; a + 1 < anchors.size(); ++a) {
        long ia = anchors[a].first, qa = anchors[a].second;
        long ib = anchors[a + 1].first, qb = anchors[a + 1].second;
        for (long i = ia; i < ib; ++i) {
            j_star[i] = qa + (i - ia) * (qb - qa) / (ib - ia);
        }
    }
    long i_last = anchors.back().first;
    long q_last = anchors.back().second;
    for (long i = i_last; i < (long)score_rows; ++i) {
        j_star[i] = std::min((long)q_len, q_last + 3 * (i - i_last));
    }

    j_min.resize(score_rows);
    j_max.resize(score_rows);
    for (long i = 0; i < (long)score_rows; ++i) {
        j_min[i] = std::max(0L, j_star[i] - W);
        j_max[i] = std::min((long)q_len, j_star[i] + W);
    }

    return true;
}


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

#define HY_MAXIMUM_RESOLUTIONS  8

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

static long const codon_spec_3x4[ 4 ][ 3 ] = {
    { 4, 3, 2 }, // 1110
    { 4, 3, 1 }, // 1101
    { 4, 2, 1 }, // 1011
    { 3, 2, 1 }  // 0111
};

static inline long construct_resolutions_helper (long c1, long c2, long c3, long * store, const long * resolutions, const long char_count) {
    if (resolutions) {
        long res1[4] = {0}, res2[4] = {0}, res3[4] = {0};
        long res_count [3] = {0,0,0};
        if (c1 >= 0) {
            res1[c1] = 1;
            res_count [0] = 1;
        } else {
            if (c1 == -1) return 0;
            c1 = -c1 - 2;
            res1[0] = resolutions[c1*4];
            res1[1] = resolutions[c1*4+1];
            res1[2] = resolutions[c1*4+2];
            res1[3] = resolutions[c1*4+3];
            res_count [0] = res1[0] + res1[1] + res1[2] + res1[3];
        }
        if (c2 >= 0) {
            res2[c2] = 1;
            res_count [1] = 1;
        } else {
            if (c2 == -1) return 0;
            c2 = -c2 - 2;
            res2[0] = resolutions[c2*4];
            res2[1] = resolutions[c2*4+1];
            res2[2] = resolutions[c2*4+2];
            res2[3] = resolutions[c2*4+3];
            res_count [1] = res2[0] + res2[1] + res2[2] + res2[3];
        }
        if (c3 >= 0) {
            res3[c3] = 1;
            res_count [2] = 1;
        } else {
            if (c3 == -1) return 0;
            c3 = -c3 - 2;
            res3[0] = resolutions[c3*4];
            res3[1] = resolutions[c3*4+1];
            res3[2] = resolutions[c3*4+2];
            res3[3] = resolutions[c3*4+3];
            res_count [2] = res3[0] + res3[1] + res3[2] + res3[3];
        }
        long total_res = res_count[0] * res_count[1] * res_count[2];
        if (total_res > HY_MAXIMUM_RESOLUTIONS || total_res == 0) return 0;
        long res_index = 0;
        for (long i = 0; i < 4; i++) {
            if (res1[i]) {
                for (long j = 0; j < 4; j++) {
                    if (res2[j]) {
                        for (long k = 0; k < 4; k++) {
                            if (res3[k]) {
                                store [res_index++] = (i*char_count + j)*char_count + k;
                            }
                        }
                    }
                }
            }
        }
        return res_index;
    }
    return 0;
}

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
                          , cawlign_fp& score
                          , long const * resolutions
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
    
    auto construct_codon_code = [&] (long c1, long c2, long c3) -> long {
        if (c1 >= 0 && c2 >= 0 && c3 >= 0) {
            return ((c1 * char_count + c2) * char_count + c3);
        }
        return -1;
    };
    
    auto construct_resolutions = [&] (long c1, long c2, long c3, long * store) -> long {
        return construct_resolutions_helper (c1, c2, c3, store, resolutions, char_count);
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
            
            //printf ("%ld/%ld %g %g \n", r, q, deletion_matrix[ prev ] , deletion_matrix[ curr ]);
        } else {
            choices[ HY_111_000 ] = score_matrix[ prev ] - open_deletion;
        }

        r_codon = construct_codon_code ( reference[ rpos - 3 ] ,reference[ rpos - 2 ], reference[ rpos - 1 ]) ;

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

        q_codon = construct_codon_code ( query[ q - 3 ], query[ q - 2 ], query[ q - 1 ]);
    }

    // if q_codon and r_codon both exist, set the score equal to match
    if ( r_codon >= 0 && q >= 3) {
        if (q_codon < 0) {
            // check the resolutions if they are available
            long resolutions [HY_MAXIMUM_RESOLUTIONS],
            res_count = construct_resolutions (query[ q - 3 ], query[ q - 2 ], query[ q - 1 ], resolutions);
            if (res_count) {
                cawlign_fp score = -1.e100;
                for (long i = 0; i < res_count; i++) {
                   cawlign_fp score_res = cost_matrix[ r_codon * cost_stride + resolutions[i]];
                    if (score_res > score) {
                        score = score_res;
                    }
                }
                //score /= res_count;
                choices[ HY_111_111 ] = score_matrix[ prev - 3 ] + score;
            } else {
                choices[ HY_111_111 ] = score_matrix[ prev - 3 ] + cost_matrix[ r_codon * cost_stride + cost_stride - 1 ];
            }
            //q_codon = cost_stride - 1;
        } else {
            
            if (q_codon >= 0) {
                const cawlign_fp move_cost = cost_matrix[ r_codon * cost_stride + q_codon ];
                choices[ HY_111_111 ] = score_matrix[ prev - 3 ] + move_cost;
            }
        }
        //if (do_local && choices [HY_LOCAL_ALIGN_SHORTCUT] < move_cost) {
        //    local_shortcut_came_from_this_move = HY_111_111;
        //    choices [HY_LOCAL_ALIGN_SHORTCUT]  = move_cost;
        //}
    }

    // we disallow partial moves in the reference, so those used to be here but are now gone

    // HERE BE DRAGONS!!!!

    // miscall matches, starting with 3x5, then 3x4, then 3x2, finally 3x1
    if ( r_codon >= 0 ) {
#if defined(__ARM_NEON2__)
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
            for (long i = 0; i < HY_3X5_COUNT; i+=4 ) {
                 cawlign_fp p[4];
                 long       c[4];
                 
                 for (long k=0; k<4; k++) {
                    if (i+k < HY_3X5_COUNT) {
                        c[k] = HY_3X5_START + i + k;
                        if ( ( q == 5              && c[k] == HY_00111_11111 )
                          || ( q == score_cols - 1 && c[k] == HY_11100_11111 ) )
                            p[k] = 0.;
                        // if we have a single ragged edge, penalize by a single miscall
                        // we don't have to worry about specifying each case here,
                        // as the 00111_11111 case takes preference above,
                        // so we don't have to explicitly avoid it
                        else if ( q == 5 && c[k] >= HY_01110_11111 )
                            p[k] = miscall_cost;
                        // if we have a single ragged edge, penalize by a single miscall
                        // unfortunately these cases are spread out,
                        // so we have to enumerate them explicitly here
                        else if ( q == score_cols - 1
                               && ( c[k] == HY_11010_11111
                                 || c[k] == HY_10110_11111
                                 || c[k] == HY_01110_11111 ) )
                            p[k] = miscall_cost;
                        // if we don't fit into any of these special cases,
                        // the miscall penalty is cawlign_fp (as we're matching 3 to 5)
                        else
                            p[k] = 2. * miscall_cost;
                    }
                 }
                
                float32x4_t penalties = vld1q_f32(p);
                float32x4_t move_costs = vld1q_f32(codon3x5 + r_codon * offset3x5 + HY_3X5_COUNT * partial_codons[ i ] + i);
                float32x4_t score_matrix_vector = vdupq_n_f32 (score_matrix[ prev - 5 ]);
                float32x4_t choice_scores = vsubq_f32(vaddq_f32(score_matrix_vector, move_costs), penalties);
                
                for (long k=0; k<4; k++) {
                    if (i+k < HY_3X5_COUNT) {
                        choices [c[k]] = choice_scores[k];
                    }
                }

            }
        }
#else
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
#endif

#if defined(__ARM_NEON2__)
        // 3x4 partial codons
        if ( q >= 4 ) {
            // fill in partial codons table
            for (long i = 0; i < HY_3X4_COUNT; ++i ) {
                partial_codons[ i ] = ( query[ q - codon_spec_3x4[ i ][ 0 ] ]   * char_count
                                      + query[ q - codon_spec_3x4[ i ][ 1 ] ] ) * char_count
                                      + query[ q - codon_spec_3x4[ i ][ 2 ] ] ;
            }
            // fill in choices
            cawlign_fp p[4];
            long       c[4];

            for (long i = 0; i < HY_3X4_COUNT; ++i ) {
                if ( partial_codons[ i ] >= 0 ) {
                    c[i] = HY_3X4_START + i;
                    // if we have a ragged edge,
                    // penalize it not at all
                    if ( ( q == 4              && c[i] == HY_0111_1111 )
                      || ( q == score_cols - 1 && c[i] == HY_1110_1111 ) )
                        p[i] = 0.;
                    // otherwise it's just a single miscall penalty
                    else
                        p[i] = miscall_cost;
                }
            }
            
            float32x4_t penalties = vld1q_f32(p);
            float32x4_t move_costs = vld1q_f32(codon3x4 + r_codon * offset3x4 + HY_3X4_COUNT * partial_codons[ 0 ] );
            float32x4_t score_matrix_vector = vdupq_n_f32 (score_matrix[ prev - 4 ]);
            float32x4_t choice_scores = vsubq_f32(vaddq_f32(score_matrix_vector, move_costs), penalties);
            
            for (long k=0; k<4; k++) {
                if (partial_codons[k] >= 0) {
                    choices [c[k]] = choice_scores[k];
                }
            }

        }
#else
        // 3x4 partial codons
        if ( q >= 4 ) {
            // fill in partial codons table
            for (long i = 0; i < HY_3X4_COUNT; ++i ) {
                partial_codons[ i ] = ( query[ q - codon_spec_3x4[ i ][ 0 ] ]   * char_count
                                      + query[ q - codon_spec_3x4[ i ][ 1 ] ] ) * char_count
                                      + query[ q - codon_spec_3x4[ i ][ 2 ] ] ;
            }
            // fill in choices
            const cawlign_fp previous_score = score_matrix[ prev - 4 ];
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

                    choices[ choice ] = previous_score - penalty + codon3x4[ r_codon * offset3x4 + HY_3X4_COUNT * partial_codons[ i ] + i ];
                    //if (do_local && choices [HY_LOCAL_ALIGN_SHORTCUT] < move_cost) {
                    //    local_shortcut_came_from_this_move = choice;
                    //    choices [HY_LOCAL_ALIGN_SHORTCUT]  = move_cost;
                    //}
                }
            }
        }
#endif
#if defined(__ARM_NEON2__)
        // 3x2
        if ( q >= 2 ) {
            // only a single partial codon
            partial_codons[ 0 ] = query[ q - 2 ] * char_count
                                + query[ q - 1 ] ;
            // fill in choices
            if ( partial_codons[ 0 ] >= 0 ) {
                cawlign_fp p[4];
                long       c[4];

                for (long i = 0; i < HY_3X2_COUNT; ++i ) {
                    c[i] = HY_3X2_START + i;
                    // if we have a ragged edge at the beginning or end,
                    // respectively, don't penalize it
                    if ( ( q == 2              && c[i] == HY_111_011 )
                      || ( q == score_cols - 1 && c[i] == HY_111_110 ) )
                        p[i] = 0.;
                    // otherwise it's just a single miscall penalty
                    else
                        p[i] = miscall_cost;
                }
                
                float32x4_t penalties = vld1q_f32(p);
                float32x4_t move_costs = vld1q_f32(codon3x2 + r_codon * offset3x2 + HY_3X2_COUNT * partial_codons[ 0 ]);
                float32x4_t score_matrix_vector = vdupq_n_f32 (score_matrix[ prev - 2 ]);
                float32x4_t choice_scores = vsubq_f32(vaddq_f32(score_matrix_vector, move_costs), penalties);

                for (long k=0; k<HY_3X2_COUNT; k++) {
                    choices [c[k]] = choice_scores[k];
                }
            }
        }
#else
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
#endif
#if defined(__ARM_NEON2__)
        // 3x1
        if ( q >= 1 ) {
            // only a single partial codon
            partial_codons[ 0 ] = query[ q - 1 ];
            // fill in choices
            if ( partial_codons[ 0 ] >= 0 ) {
                cawlign_fp p[4];
                long       c[4];

                for (long i = 0; i < HY_3X1_COUNT; ++i ) {
                    c[i] = HY_3X1_START + i;
                    // if we have a cawlign_fp ragged edge,
                    // don't enforce a miscall penalty
                    if ( ( q == 1              && c[i] == HY_111_001 )
                      || ( q == score_cols - 1 && c[i] == HY_111_100 ) )
                        p[i] = 0.;
                    // if we have a single ragged edge,
                    // enforce only a single miscall penalty
                    else if ( ( q == 1              && c[i] == HY_111_010 )
                           || ( q == score_cols - 1 && c[i] == HY_111_010 ) )
                        p[i] = miscall_cost;
                    // otherwise we need a cawlign_fp miscall penalty,
                    // for the two positions we're inserting
                    else
                        p[i] = 2. * miscall_cost;
                }
                
                float32x4_t penalties = vld1q_f32(p);
                float32x4_t move_costs = vld1q_f32(codon3x1 + r_codon * offset3x1 + HY_3X1_COUNT * partial_codons[ 0 ]);
                float32x4_t score_matrix_vector = vdupq_n_f32 (score_matrix[ prev - 1 ]);
                float32x4_t choice_scores = vsubq_f32(vaddq_f32(score_matrix_vector, move_costs), penalties);

                for (long k=0; k<HY_3X1_COUNT; k++) {
                    choices [c[k]] = choice_scores[k];
                }
            }
        }
#else
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
#endif
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
    score = max_score;
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
    //printf ("%ld (%ld/%ld) %ld\n", edit_ptr, r, q, code);
    
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
 * @param resolution_map if provided, lists ambiguity resolutions for codon-aware alignments
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
                   , const long* resolution_map
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

            uint8_t* move_matrix = nullptr;
            cawlign_fp * score_matrix = nullptr,
                       * insertion_matrix = nullptr,
                       * deletion_matrix = nullptr;

            if (!do_codon) {
                score_matrix = score_matrix_cache ? score_matrix_cache : new cawlign_fp[ score_rows * score_cols ];
                insertion_matrix = do_affine ? (insertion_matrix_cache ? insertion_matrix_cache : new cawlign_fp[ score_rows * score_cols ]) : NULL;
                deletion_matrix  = do_affine ? (deletion_matrix_cache ? deletion_matrix_cache : new cawlign_fp[ score_rows * score_cols ]) : NULL;
            }
            
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
            //    memset (insertion_matrix, 0, sizeof (cawlign_fp) * score_rows * score_cols);
            //    memset (deletion_matrix , 0, sizeof (cawlign_fp) * score_rows * score_cols);
            //}

            if ( ! do_codon ) {
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
                        cost = -open_insertion;
                        for (long i = 1; i < score_cols; ++i, cost -= open_insertion )
                            score_matrix[ i ] = cost;

                        cost = -open_deletion;
                        for (long i = score_cols; i < score_rows * score_cols; i += score_cols, cost -= open_deletion )
                            score_matrix[ i ] = cost;
                    }
                 } else {
                     // here we're doing a local alignment,
                     // the costs of opening a deletion or an insertion
                     // remain the same no matter how far down the ref or query
                     // we've traveled, respectively

                    if ( do_affine ) {
                        deletion_matrix [ 0 ] = 0.;
                        insertion_matrix [ 0 ] = 0.;
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
            }

            long index_R = r_len;
            long index_Q = q_len;
            bool took_local_shortcut = false;

            if ( do_codon ) {
                std::vector<cawlign_fp> row_storage(5 * score_cols + score_rows);
                cawlign_fp* sc_p = row_storage.data();
                cawlign_fp* sc_c = sc_p + score_cols;
                cawlign_fp* del_p = sc_c + score_cols;
                cawlign_fp* del_c = del_p + score_cols;
                cawlign_fp* ins_c = del_c + score_cols;
                cawlign_fp* last_col = ins_c + score_cols;

                // Precomputing 3x bounds
                const long offset3x5 = HY_3X5_COUNT * char_count * char_count * char_count;
                const long offset3x4 = HY_3X4_COUNT * char_count * char_count * char_count;
                const long offset3x2 = HY_3X2_COUNT * char_count * char_count;
                const long offset3x1 = HY_3X1_COUNT * char_count;

                cawlign_fp max_val_3x5[65], max_val_3x4[65], max_val_3x2[65], max_val_3x1[65];
                for (long c = 0; c < cost_stride; ++c) {
                    max_val_3x5[c] = -INFINITY;
                    max_val_3x4[c] = -INFINITY;
                    max_val_3x2[c] = -INFINITY;
                    max_val_3x1[c] = -INFINITY;
                    for (long k = 0; k < offset3x5; ++k) max_val_3x5[c] = MAX_OP(max_val_3x5[c], codon3x5[c * offset3x5 + k]);
                    for (long k = 0; k < offset3x4; ++k) max_val_3x4[c] = MAX_OP(max_val_3x4[c], codon3x4[c * offset3x4 + k]);
                    for (long k = 0; k < offset3x2; ++k) max_val_3x2[c] = MAX_OP(max_val_3x2[c], codon3x2[c * offset3x2 + k]);
                    for (long k = 0; k < offset3x1; ++k) max_val_3x1[c] = MAX_OP(max_val_3x1[c], codon3x1[c * offset3x1 + k]);
                }

                // Precomputing r_codons
                int16_t* const r_codons = (int16_t*)alloca(sizeof(int16_t) * score_rows);
                for (long r = 1; r < score_rows; ++r) {
                    long rpos = r * 3;
                    long c1 = r_enc[rpos - 3], c2 = r_enc[rpos - 2], c3 = r_enc[rpos - 1];
                    if (c1 >= 0 && c2 >= 0 && c3 >= 0) r_codons[r] = (c1 * char_count + c2) * char_count + c3;
                    else r_codons[r] = cost_stride - 1;
                }

                // Precomputing q_codons and ambiguity resolutions
                struct QCodonRes {
                    uint8_t count;
                    uint8_t codons[HY_MAXIMUM_RESOLUTIONS];
                };
                QCodonRes* const q_res = (QCodonRes*)alloca(sizeof(QCodonRes) * score_cols);
                int16_t* const q_codons = (int16_t*)alloca(sizeof(int16_t) * score_cols);
                for (long q = 0; q < score_cols; ++q) {
                    q_res[q].count = 0;
                    if (q >= 3) {
                        long c1 = q_enc[q - 3], c2 = q_enc[q - 2], c3 = q_enc[q - 1];
                        if (c1 >= 0 && c2 >= 0 && c3 >= 0) {
                            q_codons[q] = (c1 * char_count + c2) * char_count + c3;
                        } else {
                            long res[HY_MAXIMUM_RESOLUTIONS];
                            long res_count = construct_resolutions_helper(c1, c2, c3, res, resolution_map, char_count);
                            if (res_count == 1) {
                                q_codons[q] = (int16_t)res[0];
                            } else if (res_count > 1) {
                                q_codons[q] = -1;
                                q_res[q].count = (uint8_t)res_count;
                                for (long k = 0; k < res_count; ++k) {
                                    q_res[q].codons[k] = (uint8_t)res[k];
                                }
                            } else {
                                q_codons[q] = cost_stride - 1;
                            }
                        }
                    } else {
                        q_codons[q] = -1;
                    }
                }

                auto ref_kmer_table_ptr = GetRefKmerTable(r_str, r_len);

                constexpr long kBandWidth = 120;
                constexpr int kBandShift = 8;
                constexpr long kBandStride = 1 << kBandShift;
                static_assert(kBandStride >= 2 * kBandWidth + 1, "kBandStride must fit 2 * kBandWidth + 1");

                static thread_local std::vector<uint8_t> tl_band_move_matrix;
                static thread_local std::vector<uint8_t> tl_full_move_matrix;

                std::vector<long> j_min, j_max;
                bool do_band = (do_local && !do_true_local && do_affine);
                if (is_band_disabled()) do_band = false;
                if (do_band) {
                    do_band = BuildCorridor(q_str, q_len, r_len, score_rows, kBandWidth, *ref_kmer_table_ptr, j_min, j_max);
                }

                bool banded_done = false;
                if (do_band) {
                    const size_t band_bytes = (size_t)score_rows << kBandShift;
                    if (tl_band_move_matrix.size() < band_bytes) {
                        tl_band_move_matrix.resize(band_bytes);
                    }
                    move_matrix = tl_band_move_matrix.data();

                    // BANDED EXECUTION
                    std::fill(row_storage.begin(), row_storage.end(), -INFINITY);

                    // Row 0
                    for (long j = j_min[0]; j <= j_max[0]; ++j) {
                        sc_p[j] = 0.0f;
                        del_p[j] = -open_deletion - (j % 3 != 1 ? miscall_cost : 0.0f);
                    }
                    last_col[0] = (j_max[0] == (long)score_cols - 1) ? sc_p[score_cols - 1] : -INFINITY;

                    for (long i = 1; i < (long)score_rows; ++i) {
                        const long r_codon = r_codons[i];
                        const cawlign_fp* const r_cost_row = cost_matrix + r_codon * cost_stride;
                        const cawlign_fp* const r_3x5_row = codon3x5 + r_codon * offset3x5;
                        const cawlign_fp* const r_3x4_row = codon3x4 + r_codon * offset3x4;
                        const cawlign_fp* const r_3x2_row = codon3x2 + r_codon * offset3x2;
                        const cawlign_fp* const r_3x1_row = codon3x1 + r_codon * offset3x1;
                        const cawlign_fp m5 = max_val_3x5[r_codon];
                        const cawlign_fp m4 = max_val_3x4[r_codon];
                        const cawlign_fp m2 = max_val_3x2[r_codon];
                        const cawlign_fp m1 = max_val_3x1[r_codon];

                        const long row_curr = i << kBandShift;
                        const long j_min_i = j_min[i];

                        // Clear previously written intervals in sc_c, del_c (from row i-2) and ins_c (from row i-1)
                        if (i >= 2) {
                            long old_start = j_min[i - 2], old_end = j_max[i - 2];
                            for (long j = old_start; j <= old_end; ++j) {
                                sc_c[j] = -INFINITY;
                                del_c[j] = -INFINITY;
                            }
                        }
                        long ins_old_start = j_min[i - 1], ins_old_end = j_max[i - 1];
                        for (long j = ins_old_start; j <= ins_old_end; ++j) {
                            ins_c[j] = -INFINITY;
                        }

                        long j_start = j_min[i];
                        long j_end = j_max[i];

                        if (j_start == 0) {
                            sc_c[0] = 0.0f;
                            del_c[0] = 0.0f;
                            ins_c[0] = -open_insertion - ((i - 1) % 3 != 0 ? miscall_cost : 0.0f);
                            ins_c[1] = 0.0f;
                            ins_c[2] = 0.0f;
                        }

                        // del_c
                        const cawlign_fp ext_del = (i > 1 ? extend_deletion : open_deletion);
                        for (long j = std::max(1L, j_start); j <= j_end; ++j) {
                            del_c[j] = MAX_OP(sc_p[j] - open_deletion, del_p[j] - ext_del);
                        }

                        auto compute_cell_general = [&](long j) {
                            const cawlign_fp del_score = del_c[j];

                            cawlign_fp ins_score = -INFINITY;
                            if (j >= 3) {
                                ins_score = MAX_OP(
                                    sc_c[j - 3] - open_insertion,
                                    ins_c[j - 3] - (j > 3 ? extend_insertion : open_insertion)
                                );
                                ins_c[j] = ins_score;
                            }

                            cawlign_fp match_score = -INFINITY;
                            if (r_codon >= 0 && j >= 3) {
                                const long qc = q_codons[j];
                                if (qc >= 0) {
                                    match_score = sc_p[j - 3] + r_cost_row[qc];
                                } else {
                                    const auto& qr = q_res[j];
                                    cawlign_fp best_res = -1e100;
                                    for (uint8_t k = 0; k < qr.count; ++k) {
                                        const cawlign_fp sc = r_cost_row[qr.codons[k]];
                                        if (sc > best_res) best_res = sc;
                                    }
                                    match_score = sc_p[j - 3] + best_res;
                                }
                            }

                            cawlign_fp max_score = match_score;
                            uint8_t best_choice = HY_111_111;
                            if (del_score > max_score) {
                                max_score = del_score;
                                best_choice = HY_111_000;
                            }
                            if (ins_score > max_score) {
                                max_score = ins_score;
                                best_choice = HY_000_111;
                            }

                            if (r_codon >= 0) {
                                if (j >= 1) {
                                    const cawlign_fp min_p = (j == 1 || j == (long)score_cols - 1) ? 0.0f : miscall_cost;
                                    if (sc_p[j - 1] - min_p + m1 > max_score) {
                                        const long p0 = q_enc[j - 1];
                                        if (p0 >= 0) {
                                            const cawlign_fp prev_sc = sc_p[j - 1];
                                            for (long k = 0; k < HY_3X1_COUNT; ++k) {
                                                cawlign_fp penalty;
                                                const long choice = HY_3X1_START + k;
                                                if ((j == 1 && choice == HY_111_001) || (j == (long)score_cols - 1 && choice == HY_111_100))
                                                    penalty = 0.0f;
                                                else if ((j == 1 && choice == HY_111_010) || (j == (long)score_cols - 1 && choice == HY_111_010))
                                                    penalty = miscall_cost;
                                                else
                                                    penalty = 2.0f * miscall_cost;
                                                const cawlign_fp sc = prev_sc - penalty + r_3x1_row[HY_3X1_COUNT * p0 + k];
                                                if (sc > max_score) {
                                                    max_score = sc;
                                                    best_choice = (uint8_t)choice;
                                                }
                                            }
                                        }
                                    }
                                }

                                if (j >= 2) {
                                    const cawlign_fp min_p = (j == 2 || j == (long)score_cols - 1) ? 0.0f : miscall_cost;
                                    if (sc_p[j - 2] - min_p + m2 > max_score) {
                                        const long p0 = q_enc[j - 2] * char_count + q_enc[j - 1];
                                        if (p0 >= 0) {
                                            const cawlign_fp prev_sc = sc_p[j - 2];
                                            for (long k = 0; k < HY_3X2_COUNT; ++k) {
                                                cawlign_fp penalty;
                                                const long choice = HY_3X2_START + k;
                                                if ((j == 2 && choice == HY_111_011) || (j == (long)score_cols - 1 && choice == HY_111_110))
                                                    penalty = 0.0f;
                                                else
                                                    penalty = miscall_cost;
                                                const cawlign_fp sc = prev_sc - penalty + r_3x2_row[HY_3X2_COUNT * p0 + k];
                                                if (sc > max_score) {
                                                    max_score = sc;
                                                    best_choice = (uint8_t)choice;
                                                }
                                            }
                                        }
                                    }
                                }

                                if (j >= 4) {
                                    const cawlign_fp min_p = (j == 4 || j == (long)score_cols - 1) ? 0.0f : miscall_cost;
                                    if ((sc_p[j - 4] - min_p) + m4 > max_score) {
                                        const cawlign_fp prev_sc = sc_p[j - 4];
                                        for (long k = 0; k < HY_3X4_COUNT; ++k) {
                                            const long c1 = q_enc[j - codon_spec_3x4[k][0]];
                                            const long c2 = q_enc[j - codon_spec_3x4[k][1]];
                                            const long c3 = q_enc[j - codon_spec_3x4[k][2]];
                                            const long p = (c1 * char_count + c2) * char_count + c3;
                                            if (p >= 0) {
                                                cawlign_fp penalty = (j == 4 && (HY_3X4_START + k == HY_0111_1111)) ||
                                                                     (j == (long)score_cols - 1 && (HY_3X4_START + k == HY_1110_1111))
                                                                     ? 0.0f
                                                                     : miscall_cost;
                                                const cawlign_fp sc = prev_sc - penalty + r_3x4_row[HY_3X4_COUNT * p + k];
                                                if (sc > max_score) {
                                                    max_score = sc;
                                                    best_choice = (uint8_t)(HY_3X4_START + k);
                                                }
                                            }
                                        }
                                    }
                                }

                                if (j >= 5) {
                                    const cawlign_fp min_p = (j == 5 || j == (long)score_cols - 1) ? 0.0f : miscall_cost;
                                    if ((sc_p[j - 5] - min_p) + m5 > max_score) {
                                        const cawlign_fp prev_sc = sc_p[j - 5];
                                        for (long k = 0; k < HY_3X5_COUNT; ++k) {
                                            const long c1 = q_enc[j - codon_spec_3x5[k][0]];
                                            const long c2 = q_enc[j - codon_spec_3x5[k][1]];
                                            const long c3 = q_enc[j - codon_spec_3x5[k][2]];
                                            const long p = (c1 * char_count + c2) * char_count + c3;
                                            if (p >= 0) {
                                                cawlign_fp penalty;
                                                const long code = HY_3X5_START + k;
                                                if ((j == 5 && code == HY_00111_11111) || (j == (long)score_cols - 1 && code == HY_11100_11111)) {
                                                    penalty = 0.0f;
                                                } else if (j == 5 && code >= HY_01110_11111) {
                                                    penalty = miscall_cost;
                                                } else if (j == (long)score_cols - 1 && (code == HY_11010_11111 || code == HY_11001_11111 || code == HY_10110_11111 || code == HY_10101_11111 || code == HY_10011_11111)) {
                                                    penalty = miscall_cost;
                                                } else {
                                                    penalty = 2.0f * miscall_cost;
                                                }
                                                const cawlign_fp sc = prev_sc - penalty + r_3x5_row[HY_3X5_COUNT * p + k];
                                                if (sc > max_score) {
                                                    max_score = sc;
                                                    best_choice = (uint8_t)code;
                                                }
                                            }
                                        }
                                    }
                                }
                            }

                            sc_c[j] = max_score;

                            uint8_t move_byte = best_choice;
                            const cawlign_fp del_thresh = extend_deletion - open_deletion;
                            const cawlign_fp ins_thresh = extend_insertion - open_insertion;
                            if (max_score > -1e30f) {
                                move_byte |= ((del_score - max_score >= del_thresh) ? 0x80 : 0)
                                           | ((ins_score - max_score >= ins_thresh) ? 0x40 : 0);
                            }
                            move_matrix[row_curr + (j - j_min_i)] = move_byte;
                        };

                        for (long j = std::max(1L, j_start); j <= std::min(5L, j_end); ++j) {
                            compute_cell_general(j);
                        }

                        const cawlign_fp pen_2x = 2.0f * miscall_cost;
                        const cawlign_fp b1 = m1 - pen_2x;
                        const cawlign_fp b2 = m2 - miscall_cost;
                        const cawlign_fp b4 = m4 - miscall_cost;
                        const cawlign_fp b5 = m5 - pen_2x;
                        const cawlign_fp max_b = MAX_OP(MAX_OP(b1, b2), MAX_OP(b4, b5));
                        const cawlign_fp del_thresh = extend_deletion - open_deletion;
                        const cawlign_fp ins_thresh = extend_insertion - open_insertion;

                        long fast_start = std::max(6L, j_start);
                        long fast_end = std::min((long)score_cols - 2, j_end);
                        if (fast_start <= fast_end) {
                            for (long j = fast_start; j <= fast_end; ++j) {
                                const cawlign_fp del_score = del_c[j];

                                cawlign_fp ins_score = MAX_OP(sc_c[j - 3] - open_insertion, ins_c[j - 3] - extend_insertion);
                                ins_c[j] = ins_score;

                                cawlign_fp match_score;
                                const long qc = q_codons[j];
                                if (__builtin_expect(qc >= 0, 1)) {
                                    match_score = sc_p[j - 3] + r_cost_row[qc];
                                } else {
                                    const auto& qr = q_res[j];
                                    cawlign_fp best_res = -1e100;
                                    for (uint8_t k = 0; k < qr.count; ++k) {
                                        const cawlign_fp sc = r_cost_row[qr.codons[k]];
                                        if (sc > best_res) best_res = sc;
                                    }
                                    match_score = sc_p[j - 3] + best_res;
                                }

                                const bool del_wins = (del_score > match_score);
                                cawlign_fp max_score = del_wins ? del_score : match_score;
                                uint8_t best_choice = del_wins ? HY_111_000 : HY_111_111;
                                if (__builtin_expect(ins_score > max_score, 0)) {
                                    max_score = ins_score;
                                    best_choice = HY_000_111;
                                }

                                const cawlign_fp bound1 = (sc_p[j - 1] - pen_2x) + m1;
                                const cawlign_fp bound2 = (sc_p[j - 2] - miscall_cost) + m2;
                                const cawlign_fp bound4 = (sc_p[j - 4] - miscall_cost) + m4;
                                const cawlign_fp bound5 = (sc_p[j - 5] - pen_2x) + m5;
                                const cawlign_fp max_bound = MAX_OP(MAX_OP(bound1, bound2), MAX_OP(bound4, bound5));
                                if (__builtin_expect(max_bound > max_score, 0)) {
                                    if (bound1 > max_score) {
                                        const long p0 = q_enc[j - 1];
                                        if (p0 >= 0) {
                                            const cawlign_fp base_sc = sc_p[j - 1] - pen_2x;
                                            for (long k = 0; k < HY_3X1_COUNT; ++k) {
                                                const cawlign_fp sc = base_sc + r_3x1_row[HY_3X1_COUNT * p0 + k];
                                                if (sc > max_score) {
                                                    max_score = sc;
                                                    best_choice = (uint8_t)(HY_3X1_START + k);
                                                }
                                            }
                                        }
                                    }

                                    if (bound2 > max_score) {
                                        const long p0 = q_enc[j - 2] * char_count + q_enc[j - 1];
                                        if (p0 >= 0) {
                                            const cawlign_fp base_sc = sc_p[j - 2] - miscall_cost;
                                            for (long k = 0; k < HY_3X2_COUNT; ++k) {
                                                const cawlign_fp sc = base_sc + r_3x2_row[HY_3X2_COUNT * p0 + k];
                                                if (sc > max_score) {
                                                    max_score = sc;
                                                    best_choice = (uint8_t)(HY_3X2_START + k);
                                                }
                                            }
                                        }
                                    }

                                    if (bound4 > max_score) {
                                        const cawlign_fp base_sc = sc_p[j - 4] - miscall_cost;
                                        for (long k = 0; k < HY_3X4_COUNT; ++k) {
                                            const long c1 = q_enc[j - codon_spec_3x4[k][0]];
                                            const long c2 = q_enc[j - codon_spec_3x4[k][1]];
                                            const long c3 = q_enc[j - codon_spec_3x4[k][2]];
                                            const long p = (c1 * char_count + c2) * char_count + c3;
                                            if (p >= 0) {
                                                const cawlign_fp sc = base_sc + r_3x4_row[HY_3X4_COUNT * p + k];
                                                if (sc > max_score) {
                                                    max_score = sc;
                                                    best_choice = (uint8_t)(HY_3X4_START + k);
                                                }
                                            }
                                        }
                                    }

                                    if (bound5 > max_score) {
                                        const cawlign_fp base_sc = sc_p[j - 5] - pen_2x;
                                        for (long k = 0; k < HY_3X5_COUNT; ++k) {
                                            const long c1 = q_enc[j - codon_spec_3x5[k][0]];
                                            const long c2 = q_enc[j - codon_spec_3x5[k][1]];
                                            const long c3 = q_enc[j - codon_spec_3x5[k][2]];
                                            const long p = (c1 * char_count + c2) * char_count + c3;
                                            if (p >= 0) {
                                                const cawlign_fp sc = base_sc + r_3x5_row[HY_3X5_COUNT * p + k];
                                                if (sc > max_score) {
                                                    max_score = sc;
                                                    best_choice = (uint8_t)(HY_3X5_START + k);
                                                }
                                            }
                                        }
                                    }
                                }

                                sc_c[j] = max_score;

                                uint8_t move_byte = best_choice;
                                if (max_score > -1e30f) {
                                    move_byte |= ((del_score - max_score >= del_thresh) ? 0x80 : 0)
                                               | ((ins_score - max_score >= ins_thresh) ? 0x40 : 0);
                                }
                                move_matrix[row_curr + (j - j_min_i)] = move_byte;
                            }
                        }

                        if (j_end == (long)score_cols - 1) {
                            compute_cell_general(score_cols - 1);
                            last_col[i] = sc_c[score_cols - 1];
                        } else {
                            last_col[i] = -INFINITY;
                        }

                        if (i + 1 < (long)score_rows) {
                            std::swap(sc_p, sc_c);
                            std::swap(del_p, del_c);
                        }
                    }

                    // Terminal score & backtrack
                    score = (j_max[score_rows - 1] == (long)score_cols - 1) ? sc_c[score_cols - 1] : -INFINITY;
                    index_R = r_len;
                    index_Q = q_len;

                    for (long m = 0; m < (long)score_rows - 1; ++m) {
                        if (last_col[m] > score) {
                            score = last_col[m];
                            index_R = ref_stride * m;
                            index_Q = q_len;
                        }
                    }

                    for (long j = j_min[score_rows - 1]; j <= j_max[score_rows - 1] && j < (long)score_cols - 1; ++j) {
                        if (sc_c[j] > score) {
                            score = sc_c[j];
                            index_R = r_len;
                            index_Q = j;
                        }
                    }

                    for (long k = index_R; k < (long)r_len; ++k) edit_ops[edit_ptr++] = -1;
                    for (long k = index_Q; k < (long)q_len; ++k) edit_ops[edit_ptr++] = 1;

                    bool corridor_escaped = false;
                    long term_row = index_R / 3;
                    if (term_row < 0 || term_row >= (long)score_rows || index_Q < j_min[term_row] || index_Q > j_max[term_row]) {
                        corridor_escaped = true;
                    } else {
                        if (index_Q <= j_min[term_row] + 2 && j_min[term_row] > 0 && term_row > 2) corridor_escaped = true;
                        if (index_Q >= j_max[term_row] - 2 && j_max[term_row] < (long)q_len && term_row + 3 < (long)score_rows) corridor_escaped = true;
                    }

                    while (index_R && index_Q && (index_R >= 3 || index_Q >= 3) && !took_local_shortcut && !corridor_escaped) {
                        long cur_row = index_R / 3;
                        if (cur_row < 0 || cur_row >= (long)score_rows || index_Q < j_min[cur_row] || index_Q > j_max[cur_row]) {
                            corridor_escaped = true;
                            break;
                        }

                        if (index_Q <= j_min[cur_row] + 2 && j_min[cur_row] > 0 && cur_row > 2) {
                            corridor_escaped = true;
                            break;
                        }
                        if (index_Q >= j_max[cur_row] - 2 && j_max[cur_row] < (long)q_len && cur_row + 3 < (long)score_rows) {
                            corridor_escaped = true;
                            break;
                        }

                        const uint8_t move_byte = move_matrix[(cur_row << kBandShift) + (index_Q - j_min[cur_row])];
                        const long code = move_byte & 0x1F;

                        BacktrackAlignCodon(edit_ops, edit_ptr, index_R, index_Q, code);

                        if (index_R < 0 || index_Q < 0) {
                            corridor_escaped = true;
                            break;
                        }

                        if (code == HY_111_000) {
                            while (index_R >= 3) {
                                long r_row = index_R / 3;
                                if (r_row < 0 || r_row >= (long)score_rows || index_Q < j_min[r_row] || index_Q > j_max[r_row]) {
                                    corridor_escaped = true;
                                    break;
                                }
                                uint8_t pred_move = move_matrix[(r_row << kBandShift) + (index_Q - j_min[r_row])];
                                if (!(pred_move & 0x80)) break;
                                index_R -= 3;
                                edit_ops[edit_ptr++] = -1;
                                edit_ops[edit_ptr++] = -1;
                                edit_ops[edit_ptr++] = -1;
                            }
                        } else if (code == HY_000_111) {
                            long r_row = index_R / 3;
                            while (index_Q >= 3) {
                                if (r_row < 0 || r_row >= (long)score_rows || index_Q < j_min[r_row] || index_Q > j_max[r_row]) {
                                    corridor_escaped = true;
                                    break;
                                }
                                uint8_t pred_move = move_matrix[(r_row << kBandShift) + (index_Q - j_min[r_row])];
                                if (!(pred_move & 0x40)) break;
                                index_Q -= 3;
                                edit_ops[edit_ptr++] = 1;
                                edit_ops[edit_ptr++] = 1;
                                edit_ops[edit_ptr++] = 1;
                            }
                        }

                        if (corridor_escaped) break;
                    }

                    if (!corridor_escaped) {
                        long aligned_matches = 0;
                        for (long k = 0; k < edit_ptr; ++k) {
                            if (edit_ops[k] == 0) aligned_matches++;
                        }
                        long min_len = (long)std::min(r_len, q_len);
                        if (aligned_matches < 0.5 * min_len || score <= 0.0f) {
                            corridor_escaped = true;
                        }
                    }

                    if (!corridor_escaped) {
                        banded_done = true;
                    } else {
                        // Reset state for unbanded fallback
                        edit_ptr = 0;
                        index_R = r_len;
                        index_Q = q_len;
                        took_local_shortcut = false;
                    }
                }

                if (!banded_done) {
                    const size_t full_bytes = (size_t)score_rows * score_cols;
                    if (tl_full_move_matrix.size() < full_bytes) {
                        tl_full_move_matrix.resize(full_bytes);
                    }
                    move_matrix = tl_full_move_matrix.data();

                    // Initialize Row 0
                sc_p[0] = 0.0f;
                del_p[0] = (!do_local && do_affine) ? -open_deletion : 0.0f;
                if (do_affine) {
                    if (!do_local) {
                        cawlign_fp cost = -open_insertion;
                        for (long j = 1; j < score_cols; ++j, cost -= extend_insertion) {
                            sc_p[j] = cost;
                            del_p[j] = cost;
                        }
                    } else {
                        for (long j = 1; j < score_cols; ++j) {
                            sc_p[j] = 0.0f;
                            del_p[j] = -open_deletion - (j % 3 != 1 ? miscall_cost : 0.0f);
                        }
                    }
                } else {
                    if (!do_local) {
                        cawlign_fp cost = -open_insertion;
                        for (long j = 1; j < score_cols; ++j, cost -= open_insertion) {
                            sc_p[j] = cost - (j % 3 != 1 ? miscall_cost : 0.0f);
                            del_p[j] = 0.0f;
                        }
                    } else {
                        for (long j = 1; j < score_cols; ++j) {
                            sc_p[j] = 0.0f;
                            del_p[j] = 0.0f;
                        }
                    }
                }
                last_col[0] = sc_p[score_cols - 1];

                cawlign_fp best_true_local_score = -INFINITY;
                long best_true_local_R = r_len;
                long best_true_local_Q = q_len;



                for (long i = 1; i < score_rows; ++i) {
                    const long r_codon = r_codons[i];
                    const cawlign_fp* const r_cost_row = cost_matrix + r_codon * cost_stride;
                    const cawlign_fp* const r_3x5_row = codon3x5 + r_codon * offset3x5;
                    const cawlign_fp* const r_3x4_row = codon3x4 + r_codon * offset3x4;
                    const cawlign_fp* const r_3x2_row = codon3x2 + r_codon * offset3x2;
                    const cawlign_fp* const r_3x1_row = codon3x1 + r_codon * offset3x1;
                    const cawlign_fp m5 = max_val_3x5[r_codon];
                    const cawlign_fp m4 = max_val_3x4[r_codon];
                    const cawlign_fp m2 = max_val_3x2[r_codon];
                    const cawlign_fp m1 = max_val_3x1[r_codon];

                    const long row_curr = i * score_cols;

                    // Initialize Col 0
                    if (do_affine) {
                        if (!do_local) {
                            cawlign_fp cost = -open_deletion - (i - 1) * extend_deletion;
                            sc_c[0] = cost;
                            del_c[0] = cost;
                            ins_c[0] = cost;
                            ins_c[1] = -INFINITY;
                            ins_c[2] = -INFINITY;
                        } else {
                            sc_c[0] = 0.0f;
                            del_c[0] = 0.0f;
                            ins_c[0] = -open_insertion - ((i - 1) % 3 != 0 ? miscall_cost : 0.0f);
                            ins_c[1] = 0.0f;
                            ins_c[2] = 0.0f;
                        }
                    } else {
                        if (!do_local) {
                            cawlign_fp cost = -open_deletion - (i - 1) * open_insertion;
                            sc_c[0] = cost - ((i - 1) % 3 != 0 ? miscall_cost : 0.0f);
                            del_c[0] = 0.0f;
                            ins_c[0] = 0.0f;
                        } else {
                            sc_c[0] = 0.0f;
                            del_c[0] = 0.0f;
                            ins_c[0] = 0.0f;
                        }
                    }

                    if (do_affine) {
                        const cawlign_fp ext_del = (i > 1 ? extend_deletion : open_deletion);
                        long j = 1;
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
                        const float32x4_t v_open_del = vdupq_n_f32(open_deletion);
                        const float32x4_t v_ext_del = vdupq_n_f32(ext_del);
                        for (; j + 16 <= score_cols; j += 16) {
                            float32x4_t vp0 = vld1q_f32(sc_p + j);
                            float32x4_t vp1 = vld1q_f32(sc_p + j + 4);
                            float32x4_t vp2 = vld1q_f32(sc_p + j + 8);
                            float32x4_t vp3 = vld1q_f32(sc_p + j + 12);

                            float32x4_t vdp0 = vld1q_f32(del_p + j);
                            float32x4_t vdp1 = vld1q_f32(del_p + j + 4);
                            float32x4_t vdp2 = vld1q_f32(del_p + j + 8);
                            float32x4_t vdp3 = vld1q_f32(del_p + j + 12);

                            vst1q_f32(del_c + j,      vmaxq_f32(vsubq_f32(vp0, v_open_del), vsubq_f32(vdp0, v_ext_del)));
                            vst1q_f32(del_c + j + 4,  vmaxq_f32(vsubq_f32(vp1, v_open_del), vsubq_f32(vdp1, v_ext_del)));
                            vst1q_f32(del_c + j + 8,  vmaxq_f32(vsubq_f32(vp2, v_open_del), vsubq_f32(vdp2, v_ext_del)));
                            vst1q_f32(del_c + j + 12, vmaxq_f32(vsubq_f32(vp3, v_open_del), vsubq_f32(vdp3, v_ext_del)));
                        }
#endif
                        for (; j < score_cols; ++j) {
                            del_c[j] = MAX_OP(
                                sc_p[j] - open_deletion,
                                del_p[j] - ext_del
                            );
                        }
                    } else {
                        for (long j = 1; j < score_cols; ++j) {
                            del_c[j] = sc_p[j] - open_deletion;
                        }
                    }

                    auto compute_cell_general = [&](long j) {
                        // 1. Deletion (Move 1: HY_111_000)
                        const cawlign_fp del_score = del_c[j];

                        // 2. Insertion (Move 2: HY_000_111)
                        cawlign_fp ins_score = -INFINITY;
                        if (j >= 3) {
                            if (do_affine) {
                                ins_score = MAX_OP(
                                    sc_c[j - 3] - open_insertion,
                                    ins_c[j - 3] - (j > 3 ? extend_insertion : open_insertion)
                                );
                                ins_c[j] = ins_score;
                            } else {
                                ins_score = sc_c[j - 3] - open_insertion;
                            }
                        }

                        // 3. Match (Move 0: HY_111_111)
                        cawlign_fp match_score = -INFINITY;
                        if (r_codon >= 0 && j >= 3) {
                            const long qc = q_codons[j];
                            if (qc >= 0) {
                                match_score = sc_p[j - 3] + r_cost_row[qc];
                            } else {
                                const auto& qr = q_res[j];
                                cawlign_fp best_res = -1e100;
                                for (uint8_t k = 0; k < qr.count; ++k) {
                                    const cawlign_fp sc = r_cost_row[qr.codons[k]];
                                    if (sc > best_res) best_res = sc;
                                }
                                match_score = sc_p[j - 3] + best_res;
                            }
                        }

                        // Tie-breaking: 0 (match) > 1 (del) > 2 (ins)
                        cawlign_fp max_score = match_score;
                        uint8_t best_choice = HY_111_111;
                        if (del_score > max_score) {
                            max_score = del_score;
                            best_choice = HY_111_000;
                        }
                        if (ins_score > max_score) {
                            max_score = ins_score;
                            best_choice = HY_000_111;
                        }

                        if (r_codon >= 0) {
                            // Check 3x1
                            if (j >= 1) {
                                const cawlign_fp min_p = (j == 1 || j == score_cols - 1) ? 0.0f : miscall_cost;
                                if (sc_p[j - 1] - min_p + m1 > max_score) {
                                    const long p0 = q_enc[j - 1];
                                    if (p0 >= 0) {
                                        const cawlign_fp prev_sc = sc_p[j - 1];
                                        for (long k = 0; k < HY_3X1_COUNT; ++k) {
                                            cawlign_fp penalty;
                                            const long choice = HY_3X1_START + k;
                                            if ((j == 1 && choice == HY_111_001) || (j == score_cols - 1 && choice == HY_111_100))
                                                penalty = 0.0f;
                                            else if ((j == 1 && choice == HY_111_010) || (j == score_cols - 1 && choice == HY_111_010))
                                                penalty = miscall_cost;
                                            else
                                                penalty = 2.0f * miscall_cost;
                                            const cawlign_fp sc = prev_sc - penalty + r_3x1_row[HY_3X1_COUNT * p0 + k];
                                            if (sc > max_score) {
                                                max_score = sc;
                                                best_choice = (uint8_t)choice;
                                            }
                                        }
                                    }
                                }
                            }

                            // Check 3x2
                            if (j >= 2) {
                                const cawlign_fp min_p = (j == 2 || j == score_cols - 1) ? 0.0f : miscall_cost;
                                if (sc_p[j - 2] - min_p + m2 > max_score) {
                                    const long p0 = q_enc[j - 2] * char_count + q_enc[j - 1];
                                    if (p0 >= 0) {
                                        const cawlign_fp prev_sc = sc_p[j - 2];
                                        for (long k = 0; k < HY_3X2_COUNT; ++k) {
                                            cawlign_fp penalty;
                                            const long choice = HY_3X2_START + k;
                                            if ((j == 2 && choice == HY_111_011) || (j == score_cols - 1 && choice == HY_111_110))
                                                penalty = 0.0f;
                                            else
                                                penalty = miscall_cost;
                                            const cawlign_fp sc = prev_sc - penalty + r_3x2_row[HY_3X2_COUNT * p0 + k];
                                            if (sc > max_score) {
                                                max_score = sc;
                                                best_choice = (uint8_t)choice;
                                            }
                                        }
                                    }
                                }
                            }

                            // Check 3x4
                            if (j >= 4) {
                                const cawlign_fp min_p = (j == 4 || j == score_cols - 1) ? 0.0f : miscall_cost;
                                if (sc_p[j - 4] - min_p + m4 > max_score) {
                                    const cawlign_fp prev_sc = sc_p[j - 4];
                                    for (long k = 0; k < HY_3X4_COUNT; ++k) {
                                        const long c1 = q_enc[j - codon_spec_3x4[k][0]];
                                        const long c2 = q_enc[j - codon_spec_3x4[k][1]];
                                        const long c3 = q_enc[j - codon_spec_3x4[k][2]];
                                        const long p = (c1 * char_count + c2) * char_count + c3;
                                        if (p >= 0) {
                                            const long choice = HY_3X4_START + k;
                                            cawlign_fp penalty;
                                            if ((j == 4 && choice == HY_0111_1111) || (j == score_cols - 1 && choice == HY_1110_1111))
                                                penalty = 0.0f;
                                            else
                                                penalty = miscall_cost;
                                            const cawlign_fp sc = prev_sc - penalty + r_3x4_row[HY_3X4_COUNT * p + k];
                                            if (sc > max_score) {
                                                max_score = sc;
                                                best_choice = (uint8_t)choice;
                                            }
                                        }
                                    }
                                }
                            }

                            // Check 3x5
                            if (j >= 5) {
                                const cawlign_fp min_p = (j == 5 || j == score_cols - 1) ? 0.0f : miscall_cost;
                                if (sc_p[j - 5] - min_p + m5 > max_score) {
                                    const cawlign_fp prev_sc = sc_p[j - 5];
                                    for (long k = 0; k < HY_3X5_COUNT; ++k) {
                                        const long c1 = q_enc[j - codon_spec_3x5[k][0]];
                                        const long c2 = q_enc[j - codon_spec_3x5[k][1]];
                                        const long c3 = q_enc[j - codon_spec_3x5[k][2]];
                                        const long p = (c1 * char_count + c2) * char_count + c3;
                                        if (p >= 0) {
                                            const long choice = HY_3X5_START + k;
                                            cawlign_fp penalty;
                                            if ((j == 5 && choice == HY_00111_11111) || (j == score_cols - 1 && choice == HY_11100_11111))
                                                penalty = 0.0f;
                                            else if (j == 5 && choice >= HY_01110_11111)
                                                penalty = miscall_cost;
                                            else if (j == score_cols - 1 && (choice == HY_11010_11111 || choice == HY_10110_11111 || choice == HY_01110_11111))
                                                penalty = miscall_cost;
                                            else
                                                penalty = 2.0f * miscall_cost;
                                            const cawlign_fp sc = prev_sc - penalty + r_3x5_row[HY_3X5_COUNT * p + k];
                                            if (sc > max_score) {
                                                max_score = sc;
                                                best_choice = (uint8_t)choice;
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        sc_c[j] = max_score;

                        uint8_t move_byte = best_choice;
                        if (do_affine) {
                            if (del_score - extend_deletion >= max_score - open_deletion) move_byte |= 0x80;
                            if (ins_score - extend_insertion >= max_score - open_insertion) move_byte |= 0x40;
                        }
                        move_matrix[row_curr + j] = move_byte;

                        if (do_true_local && max_score > best_true_local_score) {
                            best_true_local_score = max_score;
                            best_true_local_R = ref_stride * i;
                            best_true_local_Q = j;
                        }
                    };

                    if (score_cols <= 7) {
                        for (long j = 1; j < score_cols; ++j) {
                            compute_cell_general(j);
                        }
                    } else {
                        for (long j = 1; j <= 5; ++j) {
                            compute_cell_general(j);
                        }

                        const cawlign_fp pen_2x = 2.0f * miscall_cost;
                        const cawlign_fp b1 = m1 - pen_2x;
                        const cawlign_fp b2 = m2 - miscall_cost;
                        const cawlign_fp b4 = m4 - miscall_cost;
                        const cawlign_fp b5 = m5 - pen_2x;
                        const cawlign_fp max_b = MAX_OP(MAX_OP(b1, b2), MAX_OP(b4, b5));
                        const cawlign_fp del_thresh = extend_deletion - open_deletion;
                        const cawlign_fp ins_thresh = extend_insertion - open_insertion;

                        for (long j = 6; j < score_cols - 1; ++j) {
                            // 1. Deletion
                            const cawlign_fp del_score = del_c[j];

                            // 2. Insertion
                            cawlign_fp ins_score = -INFINITY;
                            if (do_affine) {
                                ins_score = MAX_OP(sc_c[j - 3] - open_insertion, ins_c[j - 3] - extend_insertion);
                                ins_c[j] = ins_score;
                            } else {
                                ins_score = sc_c[j - 3] - open_insertion;
                            }

                            // 3. Match
                            cawlign_fp match_score;
                            const long qc = q_codons[j];
                            if (__builtin_expect(qc >= 0, 1)) {
                                match_score = sc_p[j - 3] + r_cost_row[qc];
                            } else {
                                const auto& qr = q_res[j];
                                cawlign_fp best_res = -1e100;
                                for (uint8_t k = 0; k < qr.count; ++k) {
                                    const cawlign_fp sc = r_cost_row[qr.codons[k]];
                                    if (sc > best_res) best_res = sc;
                                }
                                match_score = sc_p[j - 3] + best_res;
                            }

                            // Tie-breaking: 0 (match) > 1 (del) > 2 (ins)
                            const bool del_wins = (del_score > match_score);
                            cawlign_fp max_score = del_wins ? del_score : match_score;
                            uint8_t best_choice = del_wins ? HY_111_000 : HY_111_111;
                            if (__builtin_expect(ins_score > max_score, 0)) {
                                max_score = ins_score;
                                best_choice = HY_000_111;
                            }

                            const cawlign_fp bound1 = (sc_p[j - 1] - pen_2x) + m1;
                            const cawlign_fp bound2 = (sc_p[j - 2] - miscall_cost) + m2;
                            const cawlign_fp bound4 = (sc_p[j - 4] - miscall_cost) + m4;
                            const cawlign_fp bound5 = (sc_p[j - 5] - pen_2x) + m5;
                            const cawlign_fp max_bound = MAX_OP(MAX_OP(bound1, bound2), MAX_OP(bound4, bound5));
                            if (__builtin_expect(max_bound > max_score, 0)) {
                                // 3x1
                                if (bound1 > max_score) {
                                    const long p0 = q_enc[j - 1];
                                    if (p0 >= 0) {
                                        const cawlign_fp base_sc = sc_p[j - 1] - pen_2x;
                                        for (long k = 0; k < HY_3X1_COUNT; ++k) {
                                            const cawlign_fp sc = base_sc + r_3x1_row[HY_3X1_COUNT * p0 + k];
                                            if (sc > max_score) {
                                                max_score = sc;
                                                best_choice = (uint8_t)(HY_3X1_START + k);
                                            }
                                        }
                                    }
                                }

                                // 3x2
                                if (bound2 > max_score) {
                                    const long p0 = q_enc[j - 2] * char_count + q_enc[j - 1];
                                    if (p0 >= 0) {
                                        const cawlign_fp base_sc = sc_p[j - 2] - miscall_cost;
                                        for (long k = 0; k < HY_3X2_COUNT; ++k) {
                                            const cawlign_fp sc = base_sc + r_3x2_row[HY_3X2_COUNT * p0 + k];
                                            if (sc > max_score) {
                                                max_score = sc;
                                                best_choice = (uint8_t)(HY_3X2_START + k);
                                            }
                                        }
                                    }
                                }

                                // 3x4
                                if (bound4 > max_score) {
                                    const cawlign_fp base_sc = sc_p[j - 4] - miscall_cost;
                                    for (long k = 0; k < HY_3X4_COUNT; ++k) {
                                        const long c1 = q_enc[j - codon_spec_3x4[k][0]];
                                        const long c2 = q_enc[j - codon_spec_3x4[k][1]];
                                        const long c3 = q_enc[j - codon_spec_3x4[k][2]];
                                        const long p = (c1 * char_count + c2) * char_count + c3;
                                        if (p >= 0) {
                                            const cawlign_fp sc = base_sc + r_3x4_row[HY_3X4_COUNT * p + k];
                                            if (sc > max_score) {
                                                max_score = sc;
                                                best_choice = (uint8_t)(HY_3X4_START + k);
                                            }
                                        }
                                    }
                                }

                                // 3x5
                                if (bound5 > max_score) {
                                    const cawlign_fp base_sc = sc_p[j - 5] - pen_2x;
                                    for (long k = 0; k < HY_3X5_COUNT; ++k) {
                                        const long c1 = q_enc[j - codon_spec_3x5[k][0]];
                                        const long c2 = q_enc[j - codon_spec_3x5[k][1]];
                                        const long c3 = q_enc[j - codon_spec_3x5[k][2]];
                                        const long p = (c1 * char_count + c2) * char_count + c3;
                                        if (p >= 0) {
                                            const cawlign_fp sc = base_sc + r_3x5_row[HY_3X5_COUNT * p + k];
                                            if (sc > max_score) {
                                                max_score = sc;
                                                best_choice = (uint8_t)(HY_3X5_START + k);
                                            }
                                        }
                                    }
                                }
                            }

                            sc_c[j] = max_score;

                            uint8_t move_byte = best_choice;
                            if (do_affine) {
                                move_byte |= ((del_score - max_score >= del_thresh) ? 0x80 : 0)
                                           | ((ins_score - max_score >= ins_thresh) ? 0x40 : 0);
                            }
                            move_matrix[row_curr + j] = move_byte;

                            if (do_true_local && max_score > best_true_local_score) {
                                best_true_local_score = max_score;
                                best_true_local_R = ref_stride * i;
                                best_true_local_Q = j;
                            }
                        }

                        compute_cell_general(score_cols - 1);
                    }

                    last_col[i] = sc_c[score_cols - 1];

                    if (i + 1 < score_rows) {
                        std::swap(sc_p, sc_c);
                        std::swap(del_p, del_c);
                    }
                }

                // Determine terminal score and starting indices
                score = sc_c[score_cols - 1];
                if (do_true_local) {
                    if (best_true_local_score > score) {
                        score = best_true_local_score;
                        index_R = best_true_local_R;
                        index_Q = best_true_local_Q;
                    }
                } else if (do_local) {
                    for (long m = 0; m < score_rows - 1; ++m) {
                        if (last_col[m] > score) {
                            score = last_col[m];
                            index_R = ref_stride * m;
                        }
                    }

                    for (long j = 0; j < score_cols - 1; ++j) {
                        if (sc_c[j] > score) {
                            score = sc_c[j];
                            index_R = r_len;
                            index_Q = j;
                        }
                    }

                    for (long k = index_R; k < r_len; ++k) {
                        edit_ops[edit_ptr++] = -1;
                    }
                    for (long k = index_Q; k < q_len; ++k) {
                        edit_ops[edit_ptr++] = 1;
                    }
                }

                // Backtrack using move_matrix directly
                while (index_R && index_Q && (index_R >= 3 || index_Q >= 3) && !took_local_shortcut) {
                    const uint8_t move_byte = move_matrix[(index_R / 3) * score_cols + index_Q];
                    const long code = move_byte & 0x1F;

                    BacktrackAlignCodon(edit_ops, edit_ptr, index_R, index_Q, code);

                    if (index_R < 0 || index_Q < 0) {
                        return -INFINITY;
                    }

                    if (do_affine) {
                        long k = (index_R / 3) * score_cols + index_Q;
                        if (code == HY_111_000) {
                            while (index_R >= 3 && (move_matrix[k] & 0x80)) {
                                index_R -= 3;
                                edit_ops[edit_ptr++] = -1;
                                edit_ops[edit_ptr++] = -1;
                                edit_ops[edit_ptr++] = -1;
                                k -= score_cols;
                            }
                        } else if (code == HY_000_111) {
                            while (index_Q >= 3 && (move_matrix[k] & 0x40)) {
                                index_Q -= 3;
                                edit_ops[edit_ptr++] = 1;
                                edit_ops[edit_ptr++] = 1;
                                edit_ops[edit_ptr++] = 1;
                                k -= 3;
                            }
                        }
                    }
                }
            }
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
                                index_R = ref_stride * ( k / score_cols );
                            }

                        // grab the best score from the last row of the score matrix,
                        // skipping the very last entry ( we already checked it )
                        for (long k = ( score_rows - 1 ) * score_cols; k < score_rows * score_cols - 1; ++k )
                            if ( score_matrix[ k ] > score ) {
                                score = score_matrix[ k ];
                                index_R = r_len;
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

            if (!do_codon) {
                if (score_matrix && score_matrix != score_matrix_cache) {
                    delete [] score_matrix;
                }

                if ( do_affine ) {
                    if (insertion_matrix && insertion_matrix != insertion_matrix_cache) {
                        delete [] insertion_matrix;
                    }
                    if (deletion_matrix && deletion_matrix != deletion_matrix_cache) {
                        delete [] deletion_matrix;
                    }
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
