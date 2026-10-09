# Cawlign Scoring File Format

Cawlign uses plain text files to define scoring matrices for sequence alignment. These files are structured into distinct sections, each serving a specific purpose in defining the scoring parameters.

## General Structure

*   **Sections:** The file is organized into sections, indicated by a name enclosed in square brackets (e.g., `[ALPHABET]`, `[CODE]`, `[MATRIX]`, `[PARAMETERS]`).
*   **Comments:** Lines beginning with a semicolon (`;`) are treated as comments and are ignored by the parser.
*   **Key-Value Pairs:** Within each section, parameters are defined using a `key = value` syntax.

## Sections Details

### `[ALPHABET]` (For Nucleotide/DNA/RNA Scoring Matrices)

This section defines the alphabet of characters used in the scoring matrix. This is typically used for nucleotide sequences.

*   **`alphabet = <string>`**: A string containing all valid characters in the alphabet. The order of characters in this string is crucial as it defines the indexing for the `cost` matrix in the `[MATRIX]` section.

**Example (from NUC4.4):**
```
[ALPHABET]
alphabet = ATGCSWRYKMBVHD
```

### `[CODE]` (For Amino Acid/Codon Scoring Matrices)

This section is used for amino acid or codon-based scoring matrices. It defines the amino acid alphabet, codon translations, and how ambiguous characters are resolved.

*   **`aminoacids = <string>`**: A string of amino acids that defines the alphabet for the scoring matrix. This string typically includes standard amino acids, a stop codon symbol (e.g., `X`), and an unresolved character (e.g., `*`). The order is used for indexing the `cost` matrix.
*   **`translations = <comma-separated-list>`**: A comma-separated list of amino acid characters. This list corresponds to alphabetically sorted codons (e.g., AAA, AAC, AAG, AAT, ACA, etc.). Each entry in the list is the amino acid that the corresponding codon translates to. This is essential for codon-based alignment where codon scores are derived from amino acid scores.
*   **`resolutions = <comma-separated-list-of-triplets/quadruplets>`**: Defines how ambiguous characters are resolved. Each entry is a comma-separated list where the first element is the ambiguous character, followed by the characters it can resolve to.

**Example (from BLOSUM62):**
```
[CODE]
aminoacids = ARNDCQEGHILKMFPSTWYVBZX*
translations = K,N,K,N,T,T,T,T,R,S,R,S,I,I,M,I,Q,H,Q,H,P,P,P,P,R,R,R,R,L,L,L,L,E,D,E,D,A,A,A,A,G,G,G,G,V,V,V,V,X,Y,X,Y,S,S,S,S,X,C,W,C,L,F,L,F
resolutions = R,A,G,,Y,C,T,,S,C,G,,W,A,T,,M,A,C,,K,G,T,,H,A,C,T,,B,C,G,T,,V,A,C,G,,D,A,G,T,,N,A,C,G,T
```

### `[MATRIX]`

This section contains the actual scoring matrix values.

*   **`cost = <comma-separated-list-of-numbers>`**: A flattened list of numbers representing the scoring matrix. If the alphabet (from `[ALPHABET]` or `aminoacids` in `[CODE]`) has `D` characters, this list represents a `(D+1) x (D+1)` matrix in row-major order.
    *   `element (i,j)` (where `0 <= i < D`, `0 <= j < D`) is the cost for changing character `i` into character `j` (based on their order in the `alphabet` or `aminoacids` string).
    *   The `(i, D)` and `(D, j)` elements represent the cost of aligning character `i` (or `j`) against a non-gap but allowed character (e.g., an ambiguity).

**Example (snippet from BLOSUM62):**
```
[MATRIX]
cost = 4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0,-2,-1,0,-4,-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3,-1,0,-1,-4,-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3,3,0,-1,-4,-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3,4,1,-1,-4,0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4,-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2,0,3,-1,-4,-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4,0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3,-1,-2,-1,-4,-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3,0,0,-1,-4,-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3,-3,-3,-1,-4,-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1,-4,-3,-1,-4,-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2,0,1,-1,-4,-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1,-3,-1,-1,-4,-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1,-3,-3,-1,-4,-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2,-2,-1,-2,-4,1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2,0,0,0,-4,0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0,-1,-1,0,-4,-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3,-4,-3,-2,-4,-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1,-3,-2,-1,-4,0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4,-3,-2,-1,-4,-2,-1,3,4,-3,0,1,-1,0,-3,-4,0,-3,-3,-2,0,-1,-4,-3,-3,4,1,-1,-4,-1,0,0,1,-3,3,4,-2,0,-3,-3,1,-1,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4,0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2,0,0,-2,-1,-1,-1,-1,-1,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,1
```

### `[PARAMETERS]`

This section defines various parameters, primarily gap penalties. If a parameter is not explicitly defined, cawlign may derive it automatically.

*   **`open_deletion = <number>`**: The penalty for opening a deletion gap.
*   **`open_insertion = <number>`**: The penalty for opening an insertion gap.
*   **`extend_deletion = <number>`**: The penalty for extending a deletion gap by one character.
*   **`extend_insertion = <number>`**: The penalty for extending an insertion gap by one character.
*   **`frameshift_cost = <number>`**: (Specific to codon alignments) The penalty for a frameshift.

**Example (from BLOSUM62):**
```
[PARAMETERS]
open_deletion    = -1.
open_insertion   = -1.
extend_insertion = -1.
extend_deletion  = -1.
frameshift_cost  = -1.
```

**Example (from NUC4.4):**
```
[PARAMETERS]
open_deletion = 10.
extend_deletion = 0.5
open_insertion = 10.
extend_insertion = 0.5
```

## Python Script for Generating Cawlign Scoring Files

A Python script `utils/biopython_to_cawlign_scoring.py` is available to convert Biopython substitution matrices into cawlign scoring files.

### Usage

```bash
python utils/biopython_to_cawlign_scoring.py --help
```

### Examples

**To generate a protein scoring file (e.g., from BLOSUM62):**

```bash
python /Users/sergei/Development/cawlign/utils/biopython_to_cawlign_scoring.py \
    --matrix_name BLOSUM62 \
    --output_file res/scoring/my_blosum62.txt \
    --type protein \
    --open_gap_penalty -10.0 \
    --extend_gap_penalty -0.5
```

*   `--matrix_name`: The name of the Biopython substitution matrix (e.g., `BLOSUM62`, `PAM250`).
*   `--output_file`: The desired path for the output cawlign scoring file.
*   `--type protein`: Specifies that a protein scoring file should be generated.
*   `--open_gap_penalty`: (Optional) Penalty for opening a gap.
*   `--extend_gap_penalty`: (Optional) Penalty for extending a gap.

**To generate a codon scoring file (e.g., from BLOSUM62 using Genetic Code 1):**

```bash
python /Users/sergei/Development/cawlign/utils/biopython_to_cawlign_scoring.py \
    --matrix_name BLOSUM62 \
    --output_file res/scoring/my_codon_blosum62.txt \
    --type codon \
    --genetic_code 1 \
    --open_gap_penalty -10.0 \
    --extend_gap_penalty -0.5 \
    --frameshift_cost -15.0
```

*   `--matrix_name`: The name of the Biopython substitution matrix (e.g., `BLOSUM62`, `PAM250`).
*   `--output_file`: The desired path for the output cawlign scoring file.
*   `--type codon`: Specifies that a codon scoring file should be generated.
*   `--genetic_code`: The NCBI genetic code table number (e.g., `1` for Standard, `4` for Mycoplasma, etc.). This is **required** for `codon` type.
*   `--open_gap_penalty`: (Optional) Penalty for opening a gap.
*   `--extend_gap_penalty`: (Optional) Penalty for extending a gap.
*   `--frameshift_cost`: (Optional) Penalty for a frameshift (only applicable for `codon` type).

**Note on `translations` and `resolutions` in the generated files:**
*   For protein matrices, the `translations` and `resolutions` fields in the `[CODE]` section are left empty as they are not directly applicable in the same way as for codon matrices.
*   For codon matrices, `translations` are generated based on the provided genetic code. `resolutions` are currently left empty as handling ambiguous codons comprehensively would add significant complexity to this initial script.