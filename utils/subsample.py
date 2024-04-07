import os, argparse, subprocess, re, sys, json, math, itertools, csv, shutil, random, itertools
from Bio import SeqIO

    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='subset FASTA'
    )
    
    parser.add_argument(
        '-f', '--fasta',
        type=str,
        help='FASTA',
        required=True,
    )

    
    parser.add_argument(
        '-r', '--sampling_rate',
        type=float,
        help='CAWLIGN MSA',
        required=True,
    )
    
    

   
    args = parser.parse_args()
    
    with open(args.fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if random.random () < args.sampling_rate:
                print (">%s\n%s\n" % (record.id.split (' ')[0], str (record.seq)))
   
            
                
            
            
    
        
                
                            


