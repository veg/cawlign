SEQS=$1
RES_DIR=$2

SEQ_FILE=$(basename ${SEQS})

CAWLIGN_MSA=${RES_DIR}/${SEQ_FILE}.cawlign.msa
BEALIGN_BAM=${RES_DIR}/${SEQ_FILE}.bam
BEALIGN_MSA=${RES_DIR}/${SEQ_FILE}.bealign.msa
DIFFS=${RES_DIR}/${SEQ_FILE}.json
DIFFA=${RES_DIR}/${SEQ_FILE}.diff.fa

echo "Running CAWLIGN"
echo ""
echo ""

/usr/bin/time -l -h -p  cawlign -r HXB2_pol -t codon -I -s HIV_BETWEEN_F -o $CAWLIGN_MSA $SEQS

echo ""
echo ""
echo "Running BEALIGN"
echo ""
echo ""

/usr/bin/time -l -h -p  bealign -r HXB2_pol -m HIV_BETWEEN_F -K $SEQS $BEALIGN_BAM
bam2msa $BEALIGN_BAM $BEALIGN_MSA

echo ""
echo ""
echo "Running fasta_diff"
echo ""
echo 

fasta_diff -m $BEALIGN_MSA -t id_sequence $CAWLIGN_MSA > /dev/null 2>${DIFFS}
sed -iE 's/'\''/"/g' ${DIFFS}

BASEDIR=$(dirname "$0")

python3.9 ${BASEDIR}/compare-overlap.py -b $BEALIGN_MSA -c $CAWLIGN_MSA -d $DIFFS > $DIFFA