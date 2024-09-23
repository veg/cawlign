
#include <iostream>
#include "argparse.hpp"
#include "tn93_shared.h"
#include "alignment.h"
#include "scoring.hpp"

#ifdef _OPENMP
    #include <omp.h>
#endif

using namespace std;
using namespace argparse;


//---------------------------------------------------------------

/**
 * The main entry point for the cawlign tool.
 *
 * This function initializes command-line arguments, parses the reference sequence from a FASTA file,
 * and aligns query sequences against the reference using either codon or nucleotide scoring.
 * Supports parallel processing using OpenMP.
 *
 * @param argc The number of command-line arguments.
 * @param argv The array of command-line arguments.
 * @return 0 on successful execution, non-zero on error.
 */
int main (int argc, const char * argv[]) {

    args_t args = args_t (argc, argv);
    
    initAlphabets(args.data_type == protein);
    
    if (args.out_format == refalign) {
        ERROR_NO_USAGE ("This output mode is currently not implemented.");
    }

    char automatonState = 0,
         fasta_result = 2;
    // 0 - between sequences
    // 1 - reading sequence name
    // 2 - reading sequence
    
    CawalignSimpleScores* alignmentScoring = nullptr;
    
    if (args.reference == nullptr) {
        ERROR_NO_USAGE ("No reference sequence has been found.");
    }
    
    if (args.scores == nullptr) {
        if (args.data_type != nucleotide) {
            ERROR_NO_USAGE ("Default scoring is only available for nucleotide data. Please provide a suitable scoring file as a -s argument.");
        }
        alignmentScoring = new CawalignSimpleScores (kNucleotideAlphabet, kNucScoring, 10., 10., 0.5, 0.5);
    } else {
        if (args.data_type == codon) {
            alignmentScoring = new CawalignCodonScores (args.scores);
        } else {
            alignmentScoring = new CawalignSimpleScores (args.scores);
        }
    }
    
    StringBuffer refName,
                 refSequence;

    Vector       refNameLengths,
                 refSeqLengths;
    
    long referenceSequenceLength = 0;
    fasta_result = readFASTA (args.reference, automatonState, refName, refSequence, refNameLengths, refSeqLengths, referenceSequenceLength, true);
    if (fasta_result == 1) {
        ERROR_NO_USAGE ("The FASTA reference sequence could not be parsed.");
    }
    referenceSequenceLength++;
    
    if (args.data_type == codon) {
        CawalignCodonScores* scores = (CawalignCodonScores*)alignmentScoring;
        if (referenceSequenceLength % 3 != 0) {
            ERROR_NO_USAGE ("The reference sequence must have length divisible by 3 (data_type is codon).");
        }
        for (long i = 0; i < referenceSequenceLength; i+=3) {
            const long code = (validFlags[refSequence.getChar(i)]<<4) + (validFlags[refSequence.getChar(i+1)]<<2) + validFlags[refSequence.getChar(i+2)];
            if (code >= 0 && code < scores->translation_table.length()) {
                const char translation = scores->translation_table.value(code);
                if (translation == scores->stop_codon_index) {
                    ERROR_NO_USAGE ("The reference sequence must not have stop codons in it (data_type is codon).");
                }
            }
        }
    }
 
    long sequences_read    = 0,
         sequences_written = 0;
        
    
    automatonState = 0;
    fasta_result   = 2;
    
    VectorFP   scoreCache,
                   insertCache,
                   deleteCache;
    
    #pragma omp parallel shared (automatonState, fasta_result, sequences_read, sequences_written, args, refName, refSequence, alignmentScoring) private (nameLengths, seqLengths, names, sequences, scoreCache,insertCache,deleteCache)
    while (fasta_result == 2) {
        
        StringBuffer names,
                     sequences;

        Vector       nameLengths,
                     seqLengths;
            
        long         sequenceLength = 0;
        #pragma omp critical
        {
            fasta_result = readFASTA (args.input, automatonState, names, sequences, nameLengths, seqLengths, sequenceLength, true);
            sequenceLength++;
            if (fasta_result == 1) {
                ERROR_NO_USAGE ("Error reading the input FASTA file.");
            }
        }
        
        if (fasta_result == 2 || (fasta_result == 3 && names.length() > 0)) {
            // read a non-trivial sequence
            char * alignedRefSeq = nullptr,
            * alignedQrySeq = nullptr;
            
            cawlign_fp score;
            
            if (args.data_type != data_t::codon) {
                if (args.space_type == quadratic) {
                    AlignStrings(
                                         refSequence.getString(),
                                         sequences.getString(),
                                         referenceSequenceLength,
                                         sequenceLength,
                                         alignedRefSeq,
                                         alignedQrySeq,
                                         alignmentScoring->char_map,
                                         alignmentScoring->scoring_matrix.values(),
                                         alignmentScoring->D+1,
                                         alignmentScoring->gap_char,
                                         alignmentScoring->open_gap_reference,
                                         alignmentScoring->extend_gap_reference,
                                         alignmentScoring->open_gap_query,
                                         alignmentScoring->extend_gap_query,
                                         0.,
                                         args.local_option == trim,
                                         args.affine,
                                         false,
                                         alignmentScoring->D,
                                         nullptr,
                                         nullptr,
                                         nullptr,
                                         nullptr,
                                         args.local_option == local,
                                         args.out_format != refmap
                                         );
                    
                } else { // linear space
                    const unsigned long   size_allocation = sequenceLength+1;
                    
                    cawlign_fp          *data_buffers[6] = {nullptr};
                    
                    for (int i = 0; i < 6; i++) {
                        data_buffers[i] = new cawlign_fp [size_allocation] {0};
                    }
                    
                    char          *alignment_route = new char[2*size_allocation] {0};
                    long          *ops = new long [referenceSequenceLength + 2];
                    ops [0] = -1;
                    ops [referenceSequenceLength + 1] = sequenceLength;
                    
                    for (long i = 1L; i <= referenceSequenceLength; i++) {
                        ops[i] = -2;
                    }
                    
                    LinearSpaceAlign (refSequence.getString(),
                                              sequences.getString(),
                                              referenceSequenceLength,
                                              sequenceLength,
                                              alignmentScoring->char_map,
                                              alignmentScoring->scoring_matrix.values(),
                                              alignmentScoring->D+1,
                                              alignmentScoring->open_gap_reference,
                                              alignmentScoring->extend_gap_reference,
                                              alignmentScoring->open_gap_query,
                                              alignmentScoring->extend_gap_query,
                                              args.local_option == trim,
                                              args.affine,
                                              ops,
                                              score,
                                              0,
                                              referenceSequenceLength,
                                              0,
                                              sequenceLength,
                                              data_buffers,
                                              0,
                                              alignment_route);
                    
                    
                    
                    
                    StringBuffer     result1,
                    result2;
                    
                    long             last_column     = ops[referenceSequenceLength + 1];
                    
                    for (long position = referenceSequenceLength - 1L; position>=0; position--) {
                        
                        long current_column     = ops[position+1];
                        
                        if (current_column<0) {
                            if (current_column == -2) {
                                current_column = last_column;
                            } else if (current_column == -3) {
                                // find the next matched char or a -1
                                long    p   = position, s2p;
                                while (ops[p+1] < -1) {
                                    p--;
                                }
                                
                                s2p = ops[p+1];
                                
                                for (long j = last_column-1; j>s2p;) {
                                    if (args.out_format != refmap) {
                                        result1.appendChar(alignmentScoring->gap_char);
                                        result2.appendChar(sequences.getChar (j--));
                                    }
                                }
                                
                                last_column     = s2p+1;
                                
                                for (; position>p; position--) {
                                    result2.appendChar(alignmentScoring->gap_char);
                                    result1.appendChar(refSequence.getChar (position));
                                }
                                position ++;
                                continue;
                            } else {
                                for (last_column--; last_column >=0L; last_column--) {
                                    if (args.out_format != refmap) {
                                        result1.appendChar(alignmentScoring->gap_char);
                                        result2.appendChar(sequences.getChar (last_column));
                                    }
                                }
                                while (position>=0) {
                                    result2.appendChar(alignmentScoring->gap_char);
                                    result1.appendChar(refSequence.getChar (position--));
                                }
                                break;
                            }
                        }
                        
                        if (current_column == last_column) { // insert in sequence 2
                            result2.appendChar(alignmentScoring->gap_char);
                            result1.appendChar(refSequence.getChar (position));
                        } else {
                            last_column--;
                            
                            for (; last_column > current_column; last_column--) { // insert in column 1
                                if (args.out_format != refmap) {
                                    result1.appendChar(alignmentScoring->gap_char);
                                    result2.appendChar(sequences.getChar (last_column));
                                }
                            }
                            
                            result1.appendChar(refSequence.getChar (position));
                            result2.appendChar(sequences.getChar (current_column));
                        }
                    }
                    
                    for (last_column--; last_column >=0; last_column--) {
                        if (args.out_format != refmap) {
                            result1.appendChar(alignmentScoring->gap_char);
                            result2.appendChar(sequences.getChar (last_column));
                        }
                    }
                    
                    
                    result2.flip();
                    alignedQrySeq = result2.getString();
                    result2.detach();
                    
                    if (args.out_format == pairwise) {
                        result1.flip();
                        alignedRefSeq = result1.getString();
                        result1.detach();
                    }
                    
                    delete[]    alignment_route;
                    delete[]    ops;
                    for (int i = 0; i < 6; i++) {
                        delete [] data_buffers[i];
                    }
                    
                }
            } else {
 
                CawalignCodonScores* codonScoring = ( CawalignCodonScores*)alignmentScoring;
                
                long score_size = (referenceSequenceLength / 3 + 1) * (sequenceLength + 1);
                scoreCache.storeValue(0.,score_size-1);
                if (args.affine) {
                    insertCache.storeValue(0.,score_size-1);
                    deleteCache.storeValue(0.,score_size-1);
                }
                
                AlignStrings(
                                     refSequence.getString(),
                                     sequences.getString(),
                                     referenceSequenceLength,
                                     sequenceLength,
                                     alignedRefSeq,
                                     alignedQrySeq,
                                     alignmentScoring->char_map,
                                     alignmentScoring->scoring_matrix.values(),
                                     alignmentScoring->D+1,
                                     alignmentScoring->gap_char,
                                     alignmentScoring->open_gap_reference,
                                     alignmentScoring->extend_gap_reference,
                                     alignmentScoring->open_gap_query,
                                     alignmentScoring->extend_gap_query,
                                     codonScoring->frameshift_cost,
                                     args.local_option == trim,
                                     args.affine,
                                     true,
                                     4,
                                     codonScoring->s3x5.values(),
                                     codonScoring->s3x4.values(),
                                     codonScoring->s3x2.values(),
                                     codonScoring->s3x1.values(),
                                     args.local_option == local,
                                     args.out_format != refmap,
                                     scoreCache.rvalues(),
                                     insertCache.rvalues(),
                                     deleteCache.rvalues()
                                     );
                
            }
            
            if (alignedQrySeq) {
                if (args.out_format == pairwise) {
#pragma omp critical
                    {
                        fprintf (args.output, ">%s\n%s\n>%s\n%s\n", refName.getString(), alignedRefSeq, names.getString(), alignedQrySeq);
                    }
                    
                } else {
#pragma omp critical
                    {
                        if (args.include_reference) {
                            if (sequences_written == 0) {
                                fprintf (args.output, ">%s\n%s\n", refName.getString(), refSequence.getString());
                            }
                        }
                        fprintf (args.output, ">%s\n%s\n", names.getString(), alignedQrySeq);
                    }
                    
                    
#pragma omp atomic
                    sequences_written ++;
                }
                if (alignedRefSeq) {
                    delete [] (alignedRefSeq);
                }
                delete [] (alignedQrySeq);
            }
            
#pragma omp critical
            {
                sequences_read++;
                
                if (args.quiet == false && sequences_read % 100 == 0) {
                    cerr << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << setw(8) << sequences_read << " sequences";
                }
            }
            
        }
    }
    
    if (args.quiet == false) {
      cerr << endl;
    }
    
    if (alignmentScoring) {
        delete alignmentScoring;
    }
    return 0;

}

