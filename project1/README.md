# COMP 5970 Project 1
Contains implementation for Needleman-Wunsch and Smith-Waterman algorithms, using the BLOSUM 62 scoring matrix


###How to run

`python project1.py path_to_sequence1 path_to_sequence2`

Alternatively, if no arguments are supplied, it defaults to selecting the pair1 sequences.
`python project1.py`
###Output
After running the code, the global alignment for the two sequences will be printed out using Needleman-Wunsch, followed by local alignment using Smith-Waterman.   For amino acids (characters) whose matches have positive scores in BLOSUM62, '|' is used to indicate a positive match; otherwise  '*' is used. No symbol is used to denote the relationship between an amino acid and a gap

