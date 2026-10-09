#!/usr/bin/env python3
"""
Generate pathological test cases to stress-test cawlign codon-aware alignment.
Specifically targets rare branches in forward DP matrix calculations and backtracking:
1. All 24 moves:
   - 3x3 match (0)
   - 3x0 codon deletion (1)
   - 0x3 codon insertion (2)
   - 3x1 moves (3, 4, 5) [-2 frameshift]
   - 3x2 moves (6, 7, 8) [-1 frameshift]
   - 3x4 moves (9, 10, 11, 12) [+1 frameshift]
   - 3x5 moves (13..22) [+2 frameshift]
2. Ragged edge conditions (start and end of sequence: q=1,2,4,5, q=N-1, etc.)
3. Ambiguity codes (R, Y, S, W, K, M, B, D, H, V, N) inside codons & partial codons
4. In-frame stop codons vs frameshift bypasses
5. Compensatory frameshifts (+1/-1, +2/-2, etc.)
6. Long affine gaps vs frameshifts
7. Ultra-short query edge cases
"""

import sys
import os

def load_fasta(path):
    name = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith('>'):
                if name: break
                name = line[1:]
            else:
                seq.append(line)
    return ''.join(seq)

def make_cases():
    cases = {}
    
    # Use HXB2_pr as primary base (297 bp, 99 codons)
    ref = load_fasta("res/references/HXB2_pr")
    assert len(ref) % 3 == 0

    # 1. Base identical
    cases["exact_match"] = ref

    # 2. Single codon deletion at various positions (start, mid, end)
    cases["codon_del_first"] = ref[3:]
    cases["codon_del_mid"] = ref[:60] + ref[63:]
    cases["codon_del_last"] = ref[:-3]

    # 3. Single codon insertion at various positions
    cases["codon_ins_first"] = "GCG" + ref
    cases["codon_ins_mid"] = ref[:60] + "TGG" + ref[60:]
    cases["codon_ins_last"] = ref + "CAT"

    # 4. Long affine indels (in-frame)
    cases["long_del_inframe"] = ref[:30] + ref[75:] # 45 nt / 15 codons deleted
    cases["long_ins_inframe"] = ref[:45] + ("GCTGCT" * 6) + ref[45:] # 36 nt inserted

    # 5. Out-of-frame long indels (affine + frameshift)
    cases["long_del_fs_minus_1"] = ref[:30] + ref[76:] # 46 nt deleted (net -1)
    cases["long_del_fs_minus_2"] = ref[:30] + ref[77:] # 47 nt deleted (net -2)
    cases["long_ins_fs_plus_1"] = ref[:45] + ("GCTGCT" * 6 + "A") + ref[45:] # 37 nt inserted (+1)
    cases["long_ins_fs_plus_2"] = ref[:45] + ("GCTGCT" * 6 + "AA") + ref[45:] # 38 nt inserted (+2)

    # 6. Specific 3x1 moves (-2 net deletion in query within a codon)
    # Ref codon at pos 30 is ref[30:33]. Query keeps 1 nt:
    cases["fs_3x1_keep_first"] = ref[:30] + ref[30:31] + ref[33:] # keep pos 0 (HY_111_100)
    cases["fs_3x1_keep_mid"] = ref[:30] + ref[31:32] + ref[33:]   # keep pos 1 (HY_111_010)
    cases["fs_3x1_keep_last"] = ref[:30] + ref[32:33] + ref[33:]  # keep pos 2 (HY_111_001)

    # 7. Specific 3x2 moves (-1 net deletion in query within a codon)
    # Ref codon at pos 30 is ref[30:33]. Query keeps 2 nt:
    cases["fs_3x2_keep_01"] = ref[:30] + ref[30:32] + ref[33:] # drop pos 2 (HY_111_110)
    cases["fs_3x2_keep_02"] = ref[:30] + ref[30:31] + ref[32:33] + ref[33:] # drop pos 1 (HY_111_101)
    cases["fs_3x2_keep_12"] = ref[:30] + ref[31:33] + ref[33:] # drop pos 0 (HY_111_011)

    # 8. Specific 3x4 moves (+1 net insertion in query within a codon, homopolymer & non-homopolymer)
    # 4 possible positions of insertion relative to 3 ref characters:
    cases["fs_3x4_ins_pos0"] = ref[:30] + "A" + ref[30:33] + ref[33:]
    cases["fs_3x4_ins_pos1"] = ref[:30] + ref[30:31] + "A" + ref[31:33] + ref[33:]
    cases["fs_3x4_ins_pos2"] = ref[:30] + ref[30:32] + "A" + ref[32:33] + ref[33:]
    cases["fs_3x4_ins_pos3"] = ref[:30] + ref[30:33] + "A" + ref[33:]

    # 9. Specific 3x5 moves (+2 net insertion in query within a codon)
    # 10 choices of placing 2 inserted chars among the 3 reference chars:
    for i, ins in enumerate(["AA", "TT", "GG", "CC", "AT", "GC", "TA", "CG", "TG", "CA"]):
        cases[f"fs_3x5_var_{i}"] = ref[:42] + ref[42:43] + ins[:1] + ref[43:44] + ins[1:] + ref[44:45] + ref[45:]

    # 10. Ragged edge moves at start of query (q=1, 2, 4, 5)
    cases["ragged_start_q1"] = ref[2:] # 1st codon has only 1 nt
    cases["ragged_start_q2"] = ref[1:] # 1st codon has only 2 nt
    cases["ragged_start_q4"] = "A" + ref # 1st codon has 4 nt
    cases["ragged_start_q5"] = "AA" + ref # 1st codon has 5 nt

    # 11. Ragged edge moves at end of query (q = score_cols - 1 special rules)
    cases["ragged_end_minus_1"] = ref[:-1] # end is missing 1 nt
    cases["ragged_end_minus_2"] = ref[:-2] # end is missing 2 nt
    cases["ragged_end_plus_1"] = ref + "A" # end has +1 nt
    cases["ragged_end_plus_2"] = ref + "AT" # end has +2 nt
    cases["ragged_end_plus_4"] = ref + "ATCG" # end has +4 nt

    # 12. Ambiguity characters in match and frameshift positions
    # Single and multiple degenerate bases
    ambigs = ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N']
    for a in ambigs:
        # Codon with ambiguity
        cases[f"ambig_codon_{a}"] = ref[:15] + a + ref[16:]
        # 3x4 with ambiguity in the inserted base
        cases[f"ambig_3x4_{a}"] = ref[:21] + ref[21:22] + a + ref[22:24] + ref[24:]
        # 3x2 with ambiguity in remaining base
        cases[f"ambig_3x2_{a}"] = ref[:27] + a + ref[28:29] + ref[30:]
        # 3x1 with ambiguity
        cases[f"ambig_3x1_{a}"] = ref[:33] + a + ref[36:]

    # 13. Triple ambiguous codon (NNN)
    cases["ambig_nnn"] = ref[:30] + "NNN" + ref[33:]
    cases["ambig_multiple"] = ref[:30] + "RYSWKMBDHVN"[0:3] + ref[33:60] + "RYSWKMBDHVN"[3:6] + ref[63:]

    # 14. In-frame stop codons in query vs frameshift bypass
    # Introduce TAA, TAG, TGA at codon 20
    cases["stop_TAA_inframe"] = ref[:60] + "TAA" + ref[63:]
    cases["stop_TAG_inframe"] = ref[:60] + "TAG" + ref[63:]
    cases["stop_TGA_inframe"] = ref[:60] + "TGA" + ref[63:]

    # Stop codon created if in-frame, but frameshift avoids it
    # Ref codon is replaced by T-A-A with an extra base nearby
    cases["stop_with_adjacent_fs"] = ref[:60] + "TA" + "A" + "C" + ref[63:]

    # 15. Compensatory frameshifts (closely spaced +1 then -1, +2 then -2)
    # +1 insertion followed by -1 deletion 6 nt later
    cases["comp_plus1_minus1"] = ref[:60] + "A" + ref[60:66] + ref[67:]
    # -1 deletion followed by +1 insertion 6 nt later
    cases["comp_minus1_plus1"] = ref[:60] + ref[61:66] + "C" + ref[66:]
    # +2 insertion followed by -2 deletion
    cases["comp_plus2_minus2"] = ref[:60] + "AA" + ref[60:66] + ref[68:]
    # Three +1 insertions in a row (restoring reading frame after 3 codons)
    cases["comp_three_plus1"] = ref[:30] + "A" + ref[30:33] + "C" + ref[33:36] + "G" + ref[36:]

    # 16. Multiple homopolymer run errors (454 characteristic errors)
    # Create poly-A, poly-T, poly-G, poly-C runs with expansions/contractions
    poly_ref = ref[:30] + "AAAAAA" + "TTTTTT" + "GGGGGG" + "CCCCCC" + ref[54:]
    cases["poly_expand_all"] = ref[:30] + "AAAAAAA" + "TTTTTTT" + "GGGGGGG" + "CCCCCCC" + ref[54:]
    cases["poly_contract_all"] = ref[:30] + "AAAAA" + "TTTTT" + "GGGGG" + "CCCCC" + ref[54:]

    # 17. Extreme length sequences (ultra short query)
    cases["micro_1nt"] = "A"
    cases["micro_2nt"] = "AT"
    cases["micro_3nt"] = "ATG"
    cases["micro_4nt"] = "ATGC"
    cases["micro_5nt"] = "ATGCA"
    cases["micro_6nt"] = "ATGCAT"
    cases["micro_7nt"] = "ATGCATG"

    # 18. Large terminal non-homologous extensions (testing trim/local boundaries)
    cases["flank_5prime"] = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC" + ref
    cases["flank_3prime"] = ref + "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
    cases["flank_both"] = "CCCCCCCCCCCCCCCCCCCC" + ref + "GGGGGGGGGGGGGGGGGGGG"

    # 19. Highly mutated sequence (~30% divergence across the gene + interspersed indels)
    import random
    random.seed(42)
    mut_seq = list(ref)
    nucs = ['A', 'C', 'G', 'T']
    for pos in range(0, len(mut_seq), 4):
        mut_seq[pos] = random.choice([n for n in nucs if n != mut_seq[pos]])
    # Add a few frameshifts
    mut_seq.insert(35, 'T')
    mut_seq.insert(90, 'A')
    del mut_seq[150]
    cases["high_divergence_with_fs"] = ''.join(mut_seq)

    return cases

if __name__ == "__main__":
    out_dir = "test/pathological"
    os.makedirs(out_dir, exist_ok=True)
    out_fasta = os.path.join(out_dir, "pathological_cases.fa")
    
    cases = make_cases()
    with open(out_fasta, "w") as f:
        for name, seq in cases.items():
            f.write(f">{name}\n{seq}\n")
            
    print(f"Generated {len(cases)} pathological test cases in {out_fasta}")
