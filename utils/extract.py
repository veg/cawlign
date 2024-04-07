import os, argparse, subprocess, re, sys, json, math, itertools, csv, shutil, random, itertools
from Bio import SeqIO

    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='extract ID from FASTA'
    )
    
    parser.add_argument(
        '-f', '--fasta',
        type=str,
        help='FASTA',
        required=True,
    )

    
    parser.add_argument(
        '-i', '--id',
        type=str,
        help='ID',
        required=True,
    )
    

   
    args = parser.parse_args()
    
    with open(args.fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            id = record.id.split (' ')[0]
            if id == args.id:
               print (">%s\n%s\n" % (id, str (record.seq)))
               sys.exit (0)
   
            
                
            
            
    
        
                
                            


