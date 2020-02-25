# Comp 5970 Project 2
Contains implementation for ID3 decision tree offline learning algorithm on a fasta and corrosponding sa sequence. Ten attributes are chosen to create the tree with, and they are stored in `attributes.py` Running `id3.py` can be done through either a `train` or `test` mode. The `train` mode will create a decision tree based on 75% of the files given and store it in a binary file for later use. The `test` mode will test the remaining 25% of the data using the binary tree file, and print out statistics on how well it performed.

## How To Run
This project was created and run on the Auburn tux machines using Python 2.7.5
### Training:
##### If supplying a fasta and sa path:
`python id3.py path_to_fasta_directory path_to_sa_directory`
OR
`python id3.py train path_to_fasta_directory path_to_sa_directory`
##### Defaults to `fasta/` and `sa/` paths:
`python id3.py train`

### Testing:
`python id3.py test`

### Default
trains using paths `fasta/` and `sa/` directories in the same directory: 
`python id3.py`

## Output:
There is no output for `train` other than an indication that a tree has been created and stored in `tree.bin`. For `test`, accuracy, precision, recall, accuracy, and f1 were computing using formulas given in class. Sample output for `test` is given below, using the `fasta` and `sa` files given with the first 75% used as training and the other 25% being used as test:

```
>python id3.py test
Running Test mode
Precision: 0.685489
Recall: 0.719986
Accuracy: 0.707241
F1: 0.702314
```
