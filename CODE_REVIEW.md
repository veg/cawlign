# Cawlign Code Review & Test Suite Design

## 1. Architecture Review

### Core Alignment Engine
`cawlign` implements a specialized pairwise alignment algorithm optimized for viral and mitochondrial sequences where frameshifts (insertions/deletions of 1, 2, 4, or 5 nucleotides) are common.

*   **Codon Alignment (`AlignStrings`)**: Uses a 24-move state machine. This is the "secret sauce" of the tool, allowing it to recover from frameshifts while maintaining a high-quality protein-level scoring.
*   **Linear Space Alignment (`LinearSpaceAlign`)**: Implements Hirschberg's algorithm for $O(N+M)$ memory usage. Currently restricted to nucleotide and protein data.
*   **Parallelization**: Uses OpenMP `#pragma omp parallel` to process sequences in a FASTA file concurrently, making it highly efficient for large datasets.

### Findings & Fixes
*   **Scoring Matrix Mismatch**: The `SimpleScores` constructor expected a `(D+1)x(D+1)` matrix, but the bundled `BLOSUM62` resource provided a `DxD` matrix.
    *   *Fix*: Updated `src/scoring.cpp` to detect `DxD` matrices and automatically pad them with a default value (0.0) for the $(D+1)$-th (unresolved) character.
*   **Memory Management**: The code relies heavily on raw pointers (`new`/`delete`). While no immediate leaks were found in the hot path, a transition to `std::vector` or `std::unique_ptr` would improve safety.
*   **Linear/Quadratic Consistency**: Minor differences in gap placement were observed between linear and quadratic modes due to different tie-breaking rules in backtracking.

---

## 2. Test Suite Design Plan

### Goal
To verify all major run modes using real-world data from the `test/cases/` directory and ensure no regressions in alignment quality.

### Scenarios Covered
1.  **Nuc Default**: Quadratic space, global alignment, default nucleotide scoring.
2.  **Nuc Linear**: Hirschberg linear space algorithm.
3.  **Nuc Local**: Smith-Waterman local alignment.
4.  **Nuc Custom**: Loading custom scoring files (NUC4.4).
5.  **Protein**: BLOSUM62 scoring with the new matrix padding fix.
6.  **Codon Universal**: Standard genetic code with HXB2 reference.
7.  **Codon MtDNA**: Vertebrate mitochondrial genetic code.
8.  **Reverse Complement**: Annotated RC detection (`-R annotated`).
9.  **Output Formats**: `refalign`, `pairwise`, and reference inclusion (`-I`).

---

## 3. Implementation Details

### Testing Infrastructure
*   **`test/test_runner.py`**: A robust FASTA comparison tool that ignores whitespace and handles multi-line sequences to ensure stable cross-platform verification.
*   **`CMakeLists.txt` Integration**: Uses `add_test()` to register scenarios with `ctest`.
*   **`test/baselines/`**: Contains "gold standard" outputs generated after manual verification of the alignment quality.

### Code Coverage
Integrated using GCC/Clang `--coverage` flags.
*   **Flags**: `-O0 --coverage` is applied to the `cawlign` target when `ENABLE_COVERAGE=ON`.
*   **Target**: A `make coverage` target was added to automate `lcov` collection and HTML report generation.

---

## 4. Log of Changes
- **Modified**: `src/scoring.cpp` (Fixed scoring matrix dimension check).
- **Modified**: `CMakeLists.txt` (Added CTest and Coverage support).
- **Added**: `test/test_runner.py` (Custom verification logic).
- **Added**: `test/baselines/*.fa` (Gold standard results).
- **Added**: `test/baselines/generate_baselines.sh` (Automation for baseline updates).
- **Modified**: `.gitignore` (Added coverage artifacts).
