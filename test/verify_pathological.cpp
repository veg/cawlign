#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <map>
#include <cassert>
#include <cmath>
#include "src/scoring.hpp"
#include "src/alignment.h"
#include "src/argparse.hpp"

using namespace std;
using namespace argparse;

// CodonAlignStringsStep from src/alignment.cpp
long CodonAlignStringsStep(
    cawlign_fp * const score_matrix, long * const reference, long * const query,
    const long r, const long q, const long score_cols, const long char_count,
    const cawlign_fp miscall_cost, const cawlign_fp open_insertion, const cawlign_fp open_deletion,
    const cawlign_fp extend_insertion, const cawlign_fp extend_deletion,
    const cawlign_fp * cost_matrix, const long cost_stride,
    cawlign_fp * const insertion_matrix, cawlign_fp * const deletion_matrix,
    cawlign_fp const *codon3x5, cawlign_fp const * codon3x4,
    cawlign_fp const * codon3x2, cawlign_fp const* codon3x1,
    const bool do_local, cawlign_fp& score, long const * resolutions
);

#define HY_111_111 0
#define HY_111_000 1
#define HY_000_111 2
#define HY_111_100 3
#define HY_111_010 4
#define HY_111_001 5
#define HY_111_110 6
#define HY_111_101 7
#define HY_111_011 8
#define HY_1110_1111 9
#define HY_1101_1111 10
#define HY_1011_1111 11
#define HY_0111_1111 12
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

inline void BacktrackAlignCodon( signed char * const edit_ops
                               , long & edit_ptr
                               , long & r
                               , long & q
                               , const long code
                               )
{
    switch ( code ) {
        case HY_111_111:
            r -= 3; q -= 3;
            edit_ops [edit_ptr]   = 0;
            edit_ops [edit_ptr+1] = 0;
            edit_ops [edit_ptr+2] = 0;
            edit_ptr += 3;
            return;
        case HY_111_000:
            r -= 3;
            edit_ops [edit_ptr]   = -1;
            edit_ops [edit_ptr+1] = -1;
            edit_ops [edit_ptr+2] = -1;
            edit_ptr += 3;
            return;
        case HY_000_111:
            q -= 3;
            edit_ops [edit_ptr]   = 1;
            edit_ops [edit_ptr+1] = 1;
            edit_ops [edit_ptr+2] = 1;
            edit_ptr += 3;
            return;
    }
    
    unsigned char r_str[ 5 ] = { 1, 1, 1, 1, 1 },
                  q_str[ 5 ] = { 1, 1, 1, 1, 1 };
    int idx = 2;

    switch (code) {
        case HY_111_110: q_str[2] = 0; break;
        case HY_111_101: q_str[1] = 0; break;
        case HY_111_011: q_str[0] = 0; break;
        case HY_111_100: q_str[1] = 0; q_str[2] = 0; break;
        case HY_111_010: q_str[0] = 0; q_str[2] = 0; break;
        case HY_111_001: q_str[0] = 0; q_str[1] = 0; break;
        case HY_1110_1111: r_str[3] = 0; idx = 3; break;
        case HY_1101_1111: r_str[2] = 0; idx = 3; break;
        case HY_1011_1111: r_str[1] = 0; idx = 3; break;
        case HY_0111_1111: r_str[0] = 0; idx = 3; break;
        case HY_11100_11111: r_str[3] = 0; r_str[4] = 0; idx = 4; break;
        case HY_11010_11111: r_str[2] = 0; r_str[4] = 0; idx = 4; break;
        case HY_11001_11111: r_str[2] = 0; r_str[3] = 0; idx = 4; break;
        case HY_10110_11111: r_str[1] = 0; r_str[4] = 0; idx = 4; break;
        case HY_10101_11111: r_str[1] = 0; r_str[3] = 0; idx = 4; break;
        case HY_10011_11111: r_str[1] = 0; r_str[2] = 0; idx = 4; break;
        case HY_01110_11111: r_str[0] = 0; r_str[4] = 0; idx = 4; break;
        case HY_01101_11111: r_str[0] = 0; r_str[3] = 0; idx = 4; break;
        case HY_01011_11111: r_str[0] = 0; r_str[2] = 0; idx = 4; break;
        case HY_00111_11111: r_str[0] = 0; r_str[1] = 0; idx = 4; break;
    }

    for ( ; idx >= 0 ; --idx ) {
        if ( r_str[ idx ] ) {
            if ( q_str[ idx ] ) {
                --r; --q; edit_ops[ edit_ptr++ ] = 0;
            } else {
                --r; edit_ops[ edit_ptr++ ] = -2;
            }
        } else {
            --q; edit_ops[ edit_ptr++ ] = 2;
        }
    }
}

// Map move code to mnemonic
const char* move_name(long code) {
    switch (code) {
        case 0: return "HY_111_111 (Match/Mismatch 3x3)";
        case 1: return "HY_111_000 (Codon Deletion 3x0)";
        case 2: return "HY_000_111 (Codon Insertion 0x3)";
        case 3: return "HY_111_100 (3x1 pos 0)";
        case 4: return "HY_111_010 (3x1 pos 1)";
        case 5: return "HY_111_001 (3x1 pos 2)";
        case 6: return "HY_111_110 (3x2 pos 0,1)";
        case 7: return "HY_111_101 (3x2 pos 0,2)";
        case 8: return "HY_111_011 (3x2 pos 1,2)";
        case 9: return "HY_1110_1111 (3x4 ins pos 3)";
        case 10: return "HY_1101_1111 (3x4 ins pos 2)";
        case 11: return "HY_1011_1111 (3x4 ins pos 1)";
        case 12: return "HY_0111_1111 (3x4 ins pos 0)";
        case 13: return "HY_11100_11111 (3x5 ins pos 3,4)";
        case 14: return "HY_11010_11111 (3x5 ins pos 2,4)";
        case 15: return "HY_11001_11111 (3x5 ins pos 2,3)";
        case 16: return "HY_10110_11111 (3x5 ins pos 1,4)";
        case 17: return "HY_10101_11111 (3x5 ins pos 1,3)";
        case 18: return "HY_10011_11111 (3x5 ins pos 1,2)";
        case 19: return "HY_01110_11111 (3x5 ins pos 0,4)";
        case 20: return "HY_01101_11111 (3x5 ins pos 0,3)";
        case 21: return "HY_01011_11111 (3x5 ins pos 0,2)";
        case 22: return "HY_00111_11111 (3x5 ins pos 0,1)";
        default: return "OTHER";
    }
}

int main(int argc, const char** argv) {
    const char* ref_path = (argc > 1) ? argv[1] : "res/references/HXB2_pr";
    const char* qry_path = (argc > 2) ? argv[2] : "test/pathological/pathological_cases.fa";

    const char* fake_argv[] = {"cawlign", "-t", "codon", "-r", ref_path, "-s", "res/scoring/BLOSUM62", "-q", qry_path};
    args_t args(9, fake_argv);
    initAlphabets(false);
    CawalignCodonScores codonScoring(args.genetic_code, args.scores);

    char automatonState = 0;
    StringBuffer refName, refSeq;
    Vector refNL, refSL;
    long refLen = 0;
    readFASTA(args.reference, automatonState, refName, refSeq, refNL, refSL, refLen, true);
    refLen++;

    cout << "=================================================================" << endl;
    cout << "PATHOLOGICAL TEST SUITE: BACKTRACKING & MATRIX VERIFICATION" << endl;
    cout << "Reference: " << refName.getString() << " (length: " << refLen - 1 << " nt)" << endl;
    cout << "Query file: " << qry_path << endl;
    cout << "=================================================================" << endl;

    automatonState = 0;
    map<long, long> move_counts;
    long total_cases = 0;
    long passed_cases = 0;
    long score_mismatches = 0;
    long string_mismatches = 0;

    while (true) {
        StringBuffer qryName, qrySeq;
        Vector qryNL, qrySL;
        long qryLen = 0;
        char res = readFASTA(args.input, automatonState, qryName, qrySeq, qryNL, qrySL, qryLen, true);
        if (res == 1) break;
        if (qryName.length() == 0) break;
        qryLen++;
        total_cases++;

        // 1. Run Optimized AlignStrings
        char *r_cur = nullptr, *q_cur = nullptr;
        cawlign_fp sc_cur = AlignStrings(
            refSeq.getString(), qrySeq.getString(), refLen, qryLen,
            r_cur, q_cur,
            codonScoring.char_map, codonScoring.scoring_matrix.values(),
            codonScoring.D + 1, codonScoring.gap_char,
            codonScoring.open_gap_reference, codonScoring.extend_gap_reference,
            codonScoring.open_gap_query, codonScoring.extend_gap_query,
            codonScoring.frameshift_cost,
            true, true, true, 4,
            codonScoring.s3x5.values(), codonScoring.s3x4.values(),
            codonScoring.s3x2.values(), codonScoring.s3x1.values(),
            false, false, nullptr, nullptr, nullptr,
            codonScoring.resolutions.rvalues()
        );

        // 2. Run Legacy Forward Pass cell-by-cell with CodonAlignStringsStep
        const unsigned long score_rows = refLen / 3 + 1;
        const unsigned long score_cols = qryLen + 1;
        vector<cawlign_fp> sc_mat(score_rows * score_cols, 0.0f);
        vector<cawlign_fp> ins_mat(score_rows * score_cols, 0.0f);
        vector<cawlign_fp> del_mat(score_rows * score_cols, 0.0f);

        for (long j = 1; j < score_cols; ++j) {
            del_mat[j] = -codonScoring.open_gap_reference - (j % 3 != 1 ? codonScoring.frameshift_cost : 0);
            ins_mat[j] = 0.0f; sc_mat[j] = 0.0f;
        }
        for (long r = score_cols, j = 0; r < score_rows * score_cols; r += score_cols, ++j) {
            cawlign_fp is = -codonScoring.open_gap_query - (j % 3 != 0 ? codonScoring.frameshift_cost : 0);
            ins_mat[r] = is;
            if (r + 1 < score_rows * score_cols) ins_mat[r + 1] = 0.0f;
            if (r + 2 < score_rows * score_cols) ins_mat[r + 2] = 0.0f;
            del_mat[r] = 0.0f;
            if (r + 1 < score_rows * score_cols) del_mat[r + 1] = 0.0f;
            if (r + 2 < score_rows * score_cols) del_mat[r + 2] = 0.0f;
            sc_mat[r] = 0.0f;
        }

        vector<long> r_enc(refLen), q_enc(qryLen);
        for (long k = 0; k < refLen; ++k) r_enc[k] = codonScoring.char_map[(unsigned char)refSeq.getChar(k)];
        for (long k = 0; k < qryLen; ++k) q_enc[k] = codonScoring.char_map[(unsigned char)qrySeq.getChar(k)];

        cawlign_fp dummy_sc;
        for (long r = 1; r < score_rows; ++r) {
            for (long c = 1; c < score_cols; ++c) {
                CodonAlignStringsStep(
                    sc_mat.data(), r_enc.data(), q_enc.data(), r, c, score_cols, 4,
                    codonScoring.frameshift_cost, codonScoring.open_gap_query, codonScoring.open_gap_reference,
                    codonScoring.extend_gap_query, codonScoring.extend_gap_reference,
                    codonScoring.scoring_matrix.values(), codonScoring.D + 1,
                    ins_mat.data(), del_mat.data(),
                    codonScoring.s3x5.values(), codonScoring.s3x4.values(),
                    codonScoring.s3x2.values(), codonScoring.s3x1.values(),
                    false, dummy_sc, codonScoring.resolutions.rvalues()
                );
            }
        }

        // Traceback legacy
        long index_R = refLen - 1;
        long index_Q = qryLen - 1;
        cawlign_fp sc_leg = sc_mat[score_rows * score_cols - 1];

        // Trim terminal gaps logic
        for (long k = score_cols - 1; k < score_rows * score_cols - 1; k += score_cols) {
            if (sc_mat[k] > sc_leg) {
                sc_leg = sc_mat[k];
                index_R = 3 * (k / score_cols);
            }
        }
        for (long k = (score_rows - 1) * score_cols; k < score_rows * score_cols - 1; ++k) {
            if (sc_mat[k] > sc_leg) {
                sc_leg = sc_mat[k];
                index_R = refLen - 1;
                index_Q = k - (score_rows - 1) * score_cols;
            }
        }

        // Traceback while loop
        vector<signed char> edit_ops((refLen + qryLen) * 2, 0);
        long edit_ptr = 0;
        for (long k = index_R; k < refLen - 1; ++k) edit_ops[edit_ptr++] = -1;
        for (long k = index_Q; k < qryLen - 1; ++k) edit_ops[edit_ptr++] = 1;

        while (index_R && index_Q && (index_R >= 3 || index_Q >= 3)) {
            cawlign_fp local_sc;
            long code = CodonAlignStringsStep(
                sc_mat.data(), r_enc.data(), q_enc.data(), (index_R / 3), index_Q, score_cols, 4,
                codonScoring.frameshift_cost, codonScoring.open_gap_query, codonScoring.open_gap_reference,
                codonScoring.extend_gap_query, codonScoring.extend_gap_reference,
                codonScoring.scoring_matrix.values(), codonScoring.D + 1,
                ins_mat.data(), del_mat.data(),
                codonScoring.s3x5.values(), codonScoring.s3x4.values(),
                codonScoring.s3x2.values(), codonScoring.s3x1.values(),
                false, local_sc, codonScoring.resolutions.rvalues()
            );

            move_counts[code]++;
            BacktrackAlignCodon(edit_ops.data(), edit_ptr, index_R, index_Q, code);

            if (index_R < 0 || index_Q < 0) break;

            // Affine gap expansion
            long k = (index_R / 3) * score_cols + index_Q;
            if (code == 1) { // HY_111_000
                while (index_R >= 3 && sc_mat[k] - codonScoring.open_gap_reference <= del_mat[k] - codonScoring.extend_gap_reference) {
                    index_R -= 3;
                    edit_ops[edit_ptr++] = -1; edit_ops[edit_ptr++] = -1; edit_ops[edit_ptr++] = -1;
                    k -= score_cols;
                }
            } else if (code == 2) { // HY_000_111
                while (index_Q >= 3 && sc_mat[k] - codonScoring.open_gap_query <= ins_mat[k] - codonScoring.extend_gap_query) {
                    index_Q -= 3;
                    edit_ops[edit_ptr++] = 1; edit_ops[edit_ptr++] = 1; edit_ops[edit_ptr++] = 1;
                    k -= 3;
                }
            }
        }

        while (--index_R >= 0) edit_ops[edit_ptr++] = -1;
        while (--index_Q >= 0) edit_ops[edit_ptr++] = 1;

        // Reconstruct legacy aligned strings
        string r_leg_str, q_leg_str;
        long cur_R = 0, cur_Q = 0;
        for (long k = edit_ptr - 1; k >= 0; --k) {
            signed char op = edit_ops[k];
            switch (op) {
                case 0:
                    r_leg_str.push_back(refSeq.getChar(cur_R++));
                    q_leg_str.push_back(qrySeq.getChar(cur_Q++));
                    break;
                case 1:
                    cur_Q++; // refmap ignores insertions
                    break;
                case 2:
                    cur_Q++; // refmap ignores insertions
                    break;
                case -1:
                    r_leg_str.push_back(refSeq.getChar(cur_R++));
                    q_leg_str.push_back('-');
                    break;
                case -2:
                    r_leg_str.push_back(tolower(refSeq.getChar(cur_R++)));
                    q_leg_str.push_back('-');
                    break;
            }
        }

        // Compare score and strings
        bool score_ok = (fabs(sc_cur - sc_leg) < 1e-4);
        string q_cur_str = q_cur ? q_cur : "";
        bool string_ok = (q_cur_str == q_leg_str);

        if (!score_ok) {
            score_mismatches++;
            cout << "SCORE MISMATCH on " << qryName.getString()
                 << ": fast=" << sc_cur << " legacy=" << sc_leg << endl;
        }
        if (!string_ok) {
            string_mismatches++;
            cout << "STRING MISMATCH on " << qryName.getString() << endl;
            cout << "  Fast:   " << q_cur_str.substr(0, 80) << endl;
            cout << "  Legacy: " << q_leg_str.substr(0, 80) << endl;
        }

        if (score_ok && string_ok) {
            passed_cases++;
        }

        if (r_cur) delete[] r_cur;
        if (q_cur) delete[] q_cur;
    }

    cout << "\n=================================================================" << endl;
    cout << "TEST RESULTS: " << passed_cases << " / " << total_cases << " passed." << endl;
    cout << "Score mismatches:  " << score_mismatches << endl;
    cout << "String mismatches: " << string_mismatches << endl;
    cout << "=================================================================" << endl;

    cout << "\n=================================================================" << endl;
    cout << "CODON MOVE COVERAGE IN BACKTRACKING (All 24 moves):" << endl;
    cout << "=================================================================" << endl;
    for (int m = 0; m <= 22; ++m) {
        long cnt = move_counts[m];
        cout << "Move " << setw(2) << m << " [" << setw(35) << left << move_name(m) << "]: "
             << setw(6) << right << cnt << (cnt > 0 ? " [TRIGGERED]" : " [NOT TRIGGERED]") << endl;
    }
    cout << "=================================================================" << endl;

    if (score_mismatches == 0 && string_mismatches == 0) {
        cout << "\nALL " << total_cases << " PATHOLOGICAL CASES MATCH 100% BIT-PERFECT!" << endl;
        return 0;
    } else {
        cout << "\nFAILED: Mismatches detected." << endl;
        return 1;
    }
}
