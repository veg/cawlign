#!/bin/bash

CAWLIGN=./build/cawlign
RES_DIR=./res
TEST_DIR=./test
BASELINES_DIR=$TEST_DIR/baselines

mkdir -p $BASELINES_DIR

# Nucleotide default
$CAWLIGN -r $RES_DIR/references/HXB2_pol -q $TEST_DIR/cases/KX988404.fa > $BASELINES_DIR/nuc_default.fa

# Nucleotide linear
$CAWLIGN -S linear -r $RES_DIR/references/HXB2_pol -q $TEST_DIR/cases/KX988404.fa > $BASELINES_DIR/nuc_linear.fa

# Nucleotide local
$CAWLIGN -l local -r $RES_DIR/references/HXB2_pol -q $TEST_DIR/cases/KX988404.fa > $BASELINES_DIR/nuc_local.fa

# Nucleotide custom scoring (NUC4.4)
$CAWLIGN -s $RES_DIR/scoring/NUC4.4 -r $RES_DIR/references/HXB2_pol -q $TEST_DIR/cases/KX988404.fa > $BASELINES_DIR/nuc_custom_scoring.fa

# Create a compatible protein scoring file by prepending [ALPHABET] section
# We'll use the aminoacids from BLOSUM62 as the alphabet
ALPHABET=$(grep "aminoacids =" $RES_DIR/scoring/BLOSUM62 | cut -d'=' -f2 | tr -d ' ')
cat > $TEST_DIR/prot.scoring <<EOF
[ALPHABET]
alphabet = $ALPHABET
EOF
# Append the rest of BLOSUM62 (MATRIX and PARAMETERS sections)
sed -n '/\[MATRIX\]/,$p' $RES_DIR/scoring/BLOSUM62 >> $TEST_DIR/prot.scoring

# Protein
$CAWLIGN -t protein -r $RES_DIR/references/HXB2_pr -s $TEST_DIR/prot.scoring -q $TEST_DIR/cases/OQ948666-pr.fa > $BASELINES_DIR/prot_default.fa

# Codon (universal)
$CAWLIGN -t codon -r $RES_DIR/references/HXB2_pr -s $RES_DIR/scoring/BLOSUM62 -q $TEST_DIR/cases/OQ948666-pr.fa > $BASELINES_DIR/codon_default.fa

# Codon (vertebrate-mtdna)
$CAWLIGN -t codon -c vertebrate-mtdna -r $RES_DIR/references/COXI-mouse -s $RES_DIR/scoring/BLOSUM62 -q $TEST_DIR/cases/R_norvegicus_COXI.fa > $BASELINES_DIR/codon_mtdna.fa

# Reverse complement annotated
$CAWLIGN -R annotated -r $RES_DIR/references/HXB2_pol -q $TEST_DIR/cases/KX988404-rc.fa > $BASELINES_DIR/nuc_rc.fa

# Output formats
$CAWLIGN -f refalign -r $RES_DIR/references/HXB2_pol -q $TEST_DIR/cases/KX988404.fa > $BASELINES_DIR/nuc_refalign.fa
$CAWLIGN -f pairwise -r $RES_DIR/references/HXB2_pol -q $TEST_DIR/cases/KX988404.fa > $BASELINES_DIR/nuc_pairwise.fa
$CAWLIGN -I -f refmap -r $RES_DIR/references/HXB2_pol -q $TEST_DIR/cases/KX988404.fa > $BASELINES_DIR/nuc_refmap_with_ref.fa
