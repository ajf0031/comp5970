 # COMP 5970 Project 5
Contains implementation for Lindear Regression using stochastic gradient descent using fasta, pssm, and tmalign files. Twenty five features for each protein were selected, and run for each unique pair of proteins, making a total of 50 features. The first 20 features each protein were the average proportion of the percentages in the pssm file. The next two were the predicted exposed and buried proportions using the ID3 training data from project 2. The final 3 features for each protein were the proportions of secondary structures (H, E, C) given from the GNB training data implemented in project 3. The labels for each unique pair is given by the mean of the two TM align scores. The `train` mode will generate the weights as defined by linear regression using stochastic gradient descent, and `test` will make predictions using those weights and determine the average squared error from the predicted TM score. The first 75% of unique pairs are training and the last 25% are testing data.

## How To Run
This project was created and run on the Auburn tux machines using Python 2.7.5
### Training:
##### If supplying a fasta, pssm, and tm_align path:
```python project5.py train  path_to_fasta_directory path_to_pssm_directory path_to_tm_align_directory```
##### Defaults to `fasta/`, `pssm/`, and `tmalign/` paths:
```python project5.py train```

### Testing:
##### If supplying a fasta, pssm, and tm_align path:
```python project5.py test path_to_fasta_directory path_to_pssm_directory path_to_tm_align_directory ```
##### Defaults to `fasta/`, `pssm/`, and `tmalign/` paths:
```python project5.py test```

## Output:
##### Train
```
>python project5.py train
Running Training mode
Running stochastic gradient descent with 1000000 steps and step size 0.0000100000
Iteration 0
Iteration 10000
Iteration 20000
Iteration 30000
Iteration 40000
Iteration 50000
...
Iteration 10000000
Finished training. Data stored in weights.bin
```
##### Test
```
>python project5.py test
Running Test mode
Average Squared error: 0.0022008080416
```

