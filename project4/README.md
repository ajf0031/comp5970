# COMP 5970 Project 4
Contains implementation for Logistic Regression using stochastic gradient ascent using rr, and pssm files. Two hundred features are chosen using the first 20 columns of the pssm matrix, with a sliding window of 2 on either side of the central amino acid, paired with another amino acid. Dummy values of -1 are used for edge values that do not have prior or post amino acids in each sequence. The rr files list each contact pair in a given protein sequence.. The `train` mode will generate the weights as defined by logistic regression using stochastic gradient ascent, and `test` will make predictions using the sigmoid function and calculate L10, L5, and L2 accuracy.

## How To Run
This project was created and run on the Auburn tux machines using Python 2.7.5
### Training:
##### If supplying a rr, and pssm path:
```python project4.py train  path_to_rr_directory path_to_pssm_directory```
##### Defaults to `rr/`, and `pssm/` paths:
```python project4.py train```

### Testing:
##### If supplying an rr, and pssm path:
```python project4.py test path_to_rr_directory path_to_pssm_directory ```
##### Default
trains using paths`rr/`, and `pssm/` directories in the current directory: 

```python project4.py test```

## Output:
##### Train
```
>python project4.py train
Running Training mode
Running stochastic gradient ascent with 1000000 steps and step size 0.0000100000
Iteration 0
Iteration 10000
Iteration 20000
Iteration 30000
Iteration 40000
Iteration 50000
...
Iteration 10000000
Finished training. Data stored in train.bin
```
##### Test
```
>python project4.py test
Running Test mode
L10 Accuracy: 0.65302
L5 Accuracy: 0.65502
L2 Accuracy: 0.62548
```
