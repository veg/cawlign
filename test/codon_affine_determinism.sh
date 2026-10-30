#!/bin/bash
# Regression test: codon alignment with -l local / -l global must not depend on
# uninitialised memory (record order, thread count, malloc fill pattern).
# usage: codon_affine_determinism.sh <cawlign binary> <source dir> <scratch dir>
set -u
BIN=$1; SRC=$2; TMP=$3
mkdir -p "$TMP"
CMP="python3 $SRC/test/compare_fasta.py"
ARGS=(-t codon -s "$SRC/res/scoring/BLOSUM62" -c "$SRC/res/genetic_codes/universal" -r "$SRC/res/references/HXB2_pol")

# first 40 records, plus the same records in reverse order
python3 - "$SRC/test/HIV1-pol-326.fa" "$TMP/fwd.fa" "$TMP/rev.fa" "$TMP/one.fa" <<'P'
import sys
recs=[]
for l in open(sys.argv[1]):
    if l.startswith('>'): recs.append(l)
    elif recs: recs[-1]+=l
recs=recs[:40]
open(sys.argv[2],'w').write(''.join(recs))
open(sys.argv[3],'w').write(''.join(reversed(recs)))
open(sys.argv[4],'w').write(recs[0])
P

fail=0
for mode in local global; do
    OMP_NUM_THREADS=1 "$BIN" "${ARGS[@]}" -l $mode "$TMP/fwd.fa" > "$TMP/$mode.t1.fwd" 2>/dev/null
    OMP_NUM_THREADS=1 "$BIN" "${ARGS[@]}" -l $mode "$TMP/rev.fa" > "$TMP/$mode.t1.rev" 2>/dev/null
    OMP_NUM_THREADS=8 "$BIN" "${ARGS[@]}" -l $mode "$TMP/fwd.fa" > "$TMP/$mode.t8.fwd" 2>/dev/null
    OMP_NUM_THREADS=8 "$BIN" "${ARGS[@]}" -l $mode "$TMP/rev.fa" > "$TMP/$mode.t8.rev" 2>/dev/null
    for o in t1.rev t8.fwd t8.rev; do
        if ! $CMP "$TMP/$mode.t1.fwd" "$TMP/$mode.$o" >/dev/null; then
            echo "FAIL: -l $mode differs between t1.fwd and $o"; fail=1
        fi
    done
    # single record under different malloc fill patterns
    MALLOC_PERTURB_=0   OMP_NUM_THREADS=1 "$BIN" "${ARGS[@]}" -l $mode "$TMP/one.fa" > "$TMP/$mode.p0" 2>/dev/null
    MALLOC_PERTURB_=255 OMP_NUM_THREADS=1 "$BIN" "${ARGS[@]}" -l $mode "$TMP/one.fa" > "$TMP/$mode.p255" 2>/dev/null
    if ! $CMP "$TMP/$mode.p0" "$TMP/$mode.p255" >/dev/null; then
        echo "FAIL: -l $mode differs under MALLOC_PERTURB_ 0 vs 255"; fail=1
    fi
    [ -s "$TMP/$mode.t1.fwd" ] || { echo "FAIL: -l $mode produced no output"; fail=1; }
done
exit $fail
