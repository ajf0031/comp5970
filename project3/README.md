# COMP 5970 Project 3
Contains implementation for Gaussian Naive Bayes offline learning algorithm on a directory of fasta, ss, and pssm files. One hundred features are chosen using the first 20 columns of the pssm matrix, with a sliding window of 2 on either side of the central amino acid. Dummy values of -1 are used for edge values that do not have prior or post amino acids in each sequence. The three classifiers are found in the ss directory and are either H, E, or C based on the secondary protein strctures. Running `project3.py` can be done through either a `train` or `test` mode. The `train` mode will generate the conditional means and variances of each feature divided into the three classifiers based on the first 75% of the files, and stores the values in `train.bin`. The `test` mode will test the remaining 25% of the data using the means, variances, and class probabilities generated in `train` to find the conditional probabilities of each amino acid and make a prediction on which secondary protein structure based on the Gaussian Naive Bayes model.

## How To Run
This project was created and run on the Auburn tux machines using Python 2.7.5
### Training:
##### If supplying a fasta, ss, and pssm path:
```python project3.py train path_to_fasta_directory path_to_sa_directory path_to_pssm_directory```
##### Defaults to `fasta/`, `ss/`, and `pssm/` paths:
```python project3.py train```

### Testing:
##### If supplying a fasta, ss, and pssm path:
```python project3.py test path_to_fasta_directory path_to_ss_directory path_to_pssm_directory ```
### Default
trains using paths `fasta/`, `ss/`, and `pssm/` directories in the current directory: 
```python project3.py test```

## Output:
There is no output for `train` other than an indication that the means and variances have been created and stored in `train.bin`. For `test`, accuracy is computing using Q3 accuracy. Sample output for `test` is given below, using the `fasta`, `ss`, and `pssm` files given with the first 75% used as training and the other 25% being used as test:

```
>python project3.py test
Running Test mode
q3 accuracy = 0.599138
```
