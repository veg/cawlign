import argparse
import os
from Bio.Align import substitution_matrices
from Bio.Data import CodonTable
from Bio.Seq import Seq


def load_biopython_matrix(matrix_name):
    """Loads a Biopython substitution matrix by name."""
    try:
        return substitution_matrices.load(matrix_name)
    except ValueError:
        print(f"Error: Biopython matrix '{matrix_name}' not found. "
              "Please check available matrices in Bio.Align.substitution_matrices.")
        exit(1)

def get_all_codons():
    """Generates a list of all 64 DNA codons in alphabetical order."""
    bases = ['A', 'C', 'G', 'T']
    codons = sorted([b1 + b2 + b3 for b1 in bases for b2 in bases for b3 in bases])
    return codons

def generate_protein_scoring_file(biopython_matrix, matrix_name, output_file, open_gap, extend_gap):
    """Generates a cawlign protein scoring file from a Biopython matrix."""
    amino_acids = sorted(list(biopython_matrix.alphabet))
    # Remove '*' and 'X' if present, then add them at the end as per cawlign example
    if '*' in amino_acids:
        amino_acids.remove('*')
    if 'X' in amino_acids:
        amino_acids.remove('X')
    
    # Add common cawlign specific characters if not present
    cawlign_amino_acids = "ARNDCQEGHILKMFPSTWYVBZX*"
    final_amino_acids = []
    for aa in cawlign_amino_acids:
        if aa in amino_acids or aa in ['X', '*']: # X and * are handled explicitly
            final_amino_acids.append(aa)
        elif aa not in final_amino_acids:
            final_amino_acids.append(aa)
    
    # Ensure all original matrix amino acids are included
    for aa in sorted(list(biopython_matrix.alphabet)):
        if aa not in final_amino_acids:
            final_amino_acids.append(aa)

    # Reorder to match cawlign example if possible, otherwise keep sorted
    # For simplicity, let's just use the sorted list of unique amino acids from the matrix
    # and ensure X and * are at the end if they exist in the matrix.
    
    # Let's try to match the BLOSUM62 example order as much as possible
    # ARNDCQEGHILKMFPSTWYVBZX*
    reference_order = "ARNDCQEGHILKMFPSTWYVBZX*"
    final_amino_acids_ordered = []
    for char in reference_order:
        if char in biopython_matrix.alphabet or char in ['X', '*']:
            final_amino_acids_ordered.append(char)
    
    # Add any remaining amino acids from the matrix that weren't in the reference order
    for char in sorted(list(biopython_matrix.alphabet)):
        if char not in final_amino_acids_ordered:
            final_amino_acids_ordered.append(char)

    amino_acids_str = "".join(final_amino_acids_ordered)
    
    # Generate cost matrix
    matrix_size = len(amino_acids_str)
    cost_matrix = []
    for i in range(matrix_size):
        row = []
        for j in range(matrix_size):
            aa1 = amino_acids_str[i]
            aa2 = amino_acids_str[j]
            try:
                score = int(biopython_matrix[aa1, aa2])
            except KeyError:
                # Handle cases where an amino acid might not be in the Biopython matrix
                # Assign a default low score, e.g., -4 or 0
                score = -4 
            row.append(str(score))
        cost_matrix.extend(row)
    
    # Add the (D+1) x (D+1) extension for cawlign. For protein, this is usually for ambiguities.
    # For simplicity, we'll fill the last row/column with a default low score.
    # The cawlign example uses -4 for the last row/column.
    for _ in range(matrix_size + 1): # Add D+1 elements for the last row
        cost_matrix.append(str(-4))
    
    # Now, for each existing row, add one element for the last column
    # This requires reconstructing the matrix or carefully inserting.
    # Let's just append the last column elements to the end of each row's scores
    # and then add the last row.
    
    # Re-think: The cawlign format is a flattened (D+1)x(D+1) matrix.
    # So we need to generate D rows of D elements, then add D elements for the last column,
    # then add D+1 elements for the last row.
    
    # Let's regenerate the cost_matrix with the (D+1) x (D+1) structure directly.
    cost_values = []
    for i in range(matrix_size):
        for j in range(matrix_size):
            aa1 = amino_acids_str[i]
            aa2 = amino_acids_str[j]
            try:
                score = int(biopython_matrix[aa1, aa2])
            except KeyError:
                score = -4 # Default score for missing pairs
            cost_values.append(str(score))
        # Add the (i, D) element (score against ambiguity/unknown)
        cost_values.append(str(-4)) # Default for ambiguity
    
    # Add the (D, j) elements (ambiguity/unknown against j) and (D, D)
    for _ in range(matrix_size + 1):
        cost_values.append(str(-4)) # Default for ambiguity row/column
    
    cost_str = ",".join(cost_values)

    content = f"""[CODE]
aminoacids = {amino_acids_str}
; translations and resolutions are placeholders for protein matrices
translations = 
resolutions = 

[MATRIX]
cost = {cost_str}

[PARAMETERS]
open_deletion    = {open_gap}.
open_insertion   = {open_gap}.
extend_insertion = {extend_gap}.
extend_deletion  = {extend_gap}.
frameshift_cost  = -1.
"""
    with open(output_file, 'w') as f:
        f.write(content)
    print(f"Successfully generated cawlign protein scoring file: {output_file}")


def generate_codon_scoring_file(biopython_matrix, matrix_name, output_file, genetic_code_number, open_gap, extend_gap, frameshift_cost):
    """Generates a cawlign codon scoring file from a Biopython matrix and genetic code."""
    try:
        codon_table = CodonTable.unambiguous_dna_by_id[genetic_code_number]
    except KeyError:
        print(f"Error: Genetic code {genetic_code_number} not found. "
              "Please check valid NCBI genetic code table numbers.")
        exit(1)

    all_codons = get_all_codons()
    
    # Use the amino acid alphabet from the Biopython matrix, ensuring X and * are present if needed
    amino_acids_from_matrix = sorted(list(biopython_matrix.alphabet))
    
    # Ensure the cawlign example order for amino acids if possible
    reference_order = "ARNDCQEGHILKMFPSTWYVBZX*"
    final_amino_acids_ordered = []
    for char in reference_order:
        if char in amino_acids_from_matrix or char in ['X', '*']:
            final_amino_acids_ordered.append(char)
    
    for char in amino_acids_from_matrix:
        if char not in final_amino_acids_ordered:
            final_amino_acids_ordered.append(char)

    amino_acids_str = "".join(final_amino_acids_ordered)

    # Generate translations string
    translations = []
    for codon in all_codons:
        aa = codon_table.forward_table.get(codon, 'X') # Use 'X' for stop codons
        translations.append(aa)
    translations_str = ",".join(translations)

    # Generate cost matrix (codon-codon scores)
    codon_matrix_size = len(all_codons)
    cost_values = []
    for i in range(codon_matrix_size):
        for j in range(codon_matrix_size):
            codon1 = all_codons[i]
            codon2 = all_codons[j]
            
            aa1 = codon_table.forward_table.get(codon1, 'X')
            aa2 = codon_table.forward_table.get(codon2, 'X')
            
            try:
                score = int(biopython_matrix[aa1, aa2])
            except KeyError:
                score = -4 # Default score for missing amino acid pairs
            cost_values.append(str(score))
        
        # Add the (i, D) element (score against ambiguity/unknown)
        cost_values.append(str(-4)) # Default for ambiguity
    
    # Add the (D, j) elements (ambiguity/unknown against j) and (D, D)
    for _ in range(codon_matrix_size + 1):
        cost_values.append(str(-4)) # Default for ambiguity row/column
    
    cost_str = ",".join(cost_values)

    content = f"""[CODE]
aminoacids = {amino_acids_str}
translations = {translations_str}
; resolutions are placeholders for now
resolutions = 

[MATRIX]
cost = {cost_str}

[PARAMETERS]
open_deletion    = {open_gap}.
open_insertion   = {open_gap}.
extend_insertion = {extend_gap}.
extend_deletion  = {extend_gap}.
frameshift_cost  = {frameshift_cost}.
"""
    with open(output_file, 'w') as f:
        f.write(content)
    print(f"Successfully generated cawlign codon scoring file: {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Convert Biopython substitution matrix to Cawlign scoring file.")
    parser.add_argument("--matrix_name", required=True, 
                        help="Name of the Biopython substitution matrix (e.g., 'BLOSUM62', 'PAM250').")
    parser.add_argument("--output_file", required=True, 
                        help="Path to the output cawlign scoring file.")
    parser.add_argument("--type", choices=['protein', 'codon'], required=True,
                        help="Type of scoring file to generate: 'protein' or 'codon'.")
    parser.add_argument("--genetic_code", type=int, 
                        help="NCBI genetic code table number (required for 'codon' type).")
    parser.add_argument("--open_gap_penalty", type=float, default=-1.0,
                        help="Penalty for opening a gap (deletion or insertion). Default: -1.0.")
    parser.add_argument("--extend_gap_penalty", type=float, default=-1.0,
                        help="Penalty for extending a gap (deletion or insertion). Default: -1.0.")
    parser.add_argument("--frameshift_cost", type=float, default=-1.0,
                        help="Frameshift cost (only for 'codon' type). Default: -1.0.")

    args = parser.parse_args()

    if args.type == 'codon' and args.genetic_code is None:
        parser.error("--genetic_code is required for 'codon' type.")

    biopython_matrix = load_biopython_matrix(args.matrix_name)

    if args.type == 'protein':
        generate_protein_scoring_file(biopython_matrix, args.matrix_name, args.output_file, 
                                      args.open_gap_penalty, args.extend_gap_penalty)
    elif args.type == 'codon':
        generate_codon_scoring_file(biopython_matrix, args.matrix_name, args.output_file, 
                                    args.genetic_code, args.open_gap_penalty, 
                                    args.extend_gap_penalty, args.frameshift_cost)

if __name__ == "__main__":
    main()