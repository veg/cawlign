import os, argparse, subprocess, re, sys, json, math, itertools, csv, shutil, random, itertools
from Bio import SeqIO

    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='compare network degree distros'
    )
    
    parser.add_argument(
        '-b', '--bealign',
        type=str,
        help='BEALIGN MSA',
        required=True,
    )

    
    parser.add_argument(
        '-c', '--cawlign',
        type=str,
        help='CAWLIGN MSA',
        required=True,
    )
    
    parser.add_argument(
        '-d', '--diff',
        type=str,
        help='FASTA_DIFF json',
        required=True,
    )
    

   
    args = parser.parse_args()
    
    b_seq = {}
    
    with open(args.bealign) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            b_seq [record.id] = record.seq

    c_seq = {}
    
    with open(args.cawlign) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            c_seq [record.id] = record.seq

    with open (args.diff) as handle:
        json = json.load (handle)
        
    for s_id in json["updated_sequences"]:
        print (">[bealign] %s\n%s\n" % (s_id, str(b_seq[s_id])))
        print (">[cawalign] %s\n%s\n\n" % (s_id, str(c_seq[s_id])))

   
            
                
            
            
    
        
                
                            


