# COMP 5970 Project 1
Contains implementation for Needleman-Wunsch and Smith-Waterman algorithms, using the BLOSUM 62 scoring matrix. **NOTE: The BLOSUM 62 matrix file has the comment at the top removed and the first line begins with the Protein header**


### How to Run

##### If inputting two FASTA sequences
`python project1.py path_to_sequence1 path_to_sequence2`

##### If choosing one of the three pairs provided in the directory
`python project1.py [pair1|pair2|pair3]`

##### Default: Selects the pair1 sequences by using the paths `pair1/1k4rA_dengue_virus.fasta` and `pair1/5ire_zika_virus.fasta`.
`python project1.py`

### Output
After running the code, the global alignment for the two sequences will be printed out using Needleman-Wunsch, followed by local alignment using Smith-Waterman.   For amino acids (characters) whose matches have positive scores in BLOSUM62, '|' is used to indicate a positive match; otherwise  '*' is used. No symbol is used to denote the relationship between an amino acid and a gap

