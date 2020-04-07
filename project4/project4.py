from __future__ import division
import sys
from os import listdir
from os.path import isfile, join
from math import log, exp 
import pickle
import re
import random
TRAINING_PROPORTION = .75

def main():
	mode, rr_path, pssm_path = command_line_parser()
		
	if mode == "train":
		print "Running Training mode"
		rr_list, pssm_list = read_files(mode, rr_path, pssm_path)
		pssm_list = pssm_transform(pssm_list)
		features, labels = feature_label_generation(rr_list, pssm_list)
		weights = stochastic_logistic_regression(features, labels, 10000000, .000001)	
		with open("train.bin", "wb") as f:
			pickle.dump(weights, f)
		print "Finished training. Data stored in train.bin"
	else:
		print "Running Test mode"
		with open("train.bin", "rb") as f:
			weights = pickle.load(f)
		rr_list, pssm_list = read_files(mode, rr_path, pssm_path)
		pssm_list = pssm_transform(pssm_list)
		calculate_accuracy(rr_list, pssm_list, weights)


#Transforms pssm list to contain a sliding window of the previous and post 2 sequences, -1 values are included for edge cases
def pssm_transform(pssm_list):
	dummy_values = [-1] * 20
	
	for i in range(len(pssm_list)):
		# first value in file
		j = 0
		pssm_list[i][j].extend(dummy_values)
		pssm_list[i][j].extend(dummy_values)
		pssm_list[i][j].extend(pssm_list[i][j+1])
		pssm_list[i][j].extend(pssm_list[i][j+2])
		#2nd value in file
		j = 1
		pssm_list[i][1].extend(dummy_values)
		pssm_list[i][1].extend(pssm_list[i][j-1][0:20])
		pssm_list[i][1].extend(pssm_list[i][j+1][0:20])
		pssm_list[i][1].extend(pssm_list[i][j+2][0:20])

		for j in range(2, len(pssm_list[i]) - 2): 

			pssm_list[i][j].extend(pssm_list[i][j-2][0:20])
			pssm_list[i][j].extend(pssm_list[i][j-1][0:20])
			pssm_list[i][j].extend(pssm_list[i][j+1][0:20])
			pssm_list[i][j].extend(pssm_list[i][j+2][0:20])
			
		#j = n-1
		j = len(pssm_list[i]) - 2	
		pssm_list[i][j].extend(pssm_list[i][j-2][0:20])
		pssm_list[i][j].extend(pssm_list[i][j-1][0:20])
		pssm_list[i][j].extend(pssm_list[i][j+1][0:20])
		pssm_list[i][j].extend(dummy_values)

		#j = n
		j = len(pssm_list[i]) - 1	
		pssm_list[i][j].extend(pssm_list[i][j-2][0:20])
		pssm_list[i][j].extend(pssm_list[i][j-1][0:20])
		pssm_list[i][j].extend(dummy_values)
		pssm_list[i][j].extend(dummy_values)
	return pssm_list


def command_line_parser():
	if len(sys.argv) == 1:
		mode = "train"
	else:
		mode = sys.argv[1]
		if mode != "train" and mode != "test" and len(sys.argv) == 2:
			print "Invalid mode: please use either \"train\" or \"test\""
			exit()
	
	if len(sys.argv) == 4:
		rr_path = sys.argv[2]
		pssm_path = sys.argv[3]
	else:
		rr_path = "rr"	
		pssm_path = "pssm"

	return mode, rr_path, pssm_path


def read_files(mode, rr_path, pssm_path):

	files = [f for f in listdir(pssm_path) if isfile(join(pssm_path, f))]
	files = sorted(files)
	
	if mode == "train":
		start = 0
		end = int(TRAINING_PROPORTION*len(files))

	else:
		start = int(TRAINING_PROPORTION*len(files))
		end = len(files)
			
	
	#pssm intake
	pssm_list = []

	for i in range(start, end):
		fp = open("%s/%s" % (pssm_path, files[i]), "r")
		lines = fp.readlines()
		lines = lines[3:-6]
		for j in range(len(lines)):
			lines[j].strip()
			all_numbers = map(int, re.findall('-?\d+', lines[j]))
			lines[j] = all_numbers[1:21]

		pssm_list.append(lines)
	

	#rr intake
	rr_list = []

	files = [f for f in listdir(rr_path) if isfile(join(rr_path, f))]
	files = sorted(files)
	for i in range(start + 1, end + 1):
		fp = open("%s/%s" % (rr_path, files[i]), "r")
		lines = fp.readlines()
		lines = lines[1:]
		for j in range(len(lines)):
			lines[j].strip()
			pairs = map(int, re.findall('\d+', lines[j]))
			distance = map(float, re.findall('\d*\.\d+', lines[j]))
			lines[j] = pairs[0:2]
			lines[j].extend(distance)
					
		rr_list.append(lines)
		fp.close()


	return rr_list, pssm_list 



def log_likelihood(features, target, weights):
	scores = 0
	for i in range(len(features)):
		scores = features[i] + weights[i]

	log_likelihood = 0
	for i in range(len(features)):
		log_likelhood += target * scores - log(1 + exp(scores))
	
	return log_likelihood


def stochastic_logistic_regression(features, target, steps, step_size):
	print "Running stochastic gradient ascent with %d steps and step size %.10f" % (steps, step_size)
	weights = [0 for x in range(201)]
	
	for step in xrange(steps):
		gradient = [0 for x in range(201)]
		rand_seed = random.randrange(len(features))
		scores = get_scores(features[rand_seed], weights)
		sigmoid = calculate_sigmoid(scores)
		for j in range(len(gradient)-1):
			gradient[j] += features[rand_seed][j] * (target[rand_seed] - sigmoid)
		#account for w0 which assumes feature = 1
		gradient[-1] = 1 * (target[rand_seed] - sigmoid)
		for j in range(len(gradient)):
			weights[j] += step_size * gradient[j]  
		if step % 10000 == 0:
			test = "incorrect"
			if sigmoid > .5 and target[rand_seed] == 1: 
				test = "correct"	
			elif sigmoid < .5 and target[rand_seed] == 0: 
				test = "correct"
			print "Iteration %d" % step
	return weights

def calculate_sigmoid(scores):
	return exp(scores)/(1 + exp(scores))

def get_scores(features, weights):
	scores = 0
	for i in range(len(features)):
		scores += features[i] * weights[i]
	scores += weights[-1] * 1
	return scores


def feature_label_generation(rr_list, pssm_list):
	balanced_features = []
	balanced_labels = []
	unbalanced_features = []

	#add all contacts to balanced data
	for i in range(len(rr_list)):
		for j in range(len(rr_list[i])):
			feature = pssm_list[i][rr_list[i][j][0]-1] + pssm_list[i][rr_list[i][j][1]-1]
			balanced_features.append(feature)
			balanced_labels.append(1)

	#find all non contacts in the data 
	for i in range(len(pssm_list)):
		index = 0
		for j in range(len(pssm_list[i])):
			for k in range(j + 5, len(pssm_list[i])):
				feature = []
				if k < len(pssm_list[i]):
					if j == rr_list[i][index][0] - 1 and k == rr_list[i][index][1] - 1:
						continue;
					else:
					
						feature = pssm_list[i][j] + pssm_list[i][k]
						unbalanced_features.append(feature)
						
				else:
					break
	#Add in an equal number of non-contacts
	for i in range(len(balanced_features)):
		x = random.choice(unbalanced_features)
		balanced_features.append(x)
		balanced_labels.append(0)
	
	return balanced_features, balanced_labels	

	
def calculate_accuracy(rr_list, pssm_list, weights):
	#L10
	correct = 0
	total = 0
	for i in range(len(rr_list)):
		for j in range(int(len(pssm_list[i]) / 10)):
			prediction =  - 1
			feature = pssm_list[i][rr_list[i][j][0]-1] + pssm_list[i][rr_list[i][j][1]-1]
			scores = get_scores(feature, weights)
			probability = calculate_sigmoid(scores)
			if probability < .5:
				total += 1
			else:
				correct += 1
				total += 1
	print "L10 Accuracy: %.5f " % (correct/total)
	
	#L5
	correct = 0
	total = 0
	for i in range(len(rr_list)):
		for j in range(int(len(pssm_list[i]) / 5)):
			prediction =  - 1
			feature = pssm_list[i][rr_list[i][j][0]-1] + pssm_list[i][rr_list[i][j][1]-1]
			scores = get_scores(feature, weights)
			probability = calculate_sigmoid(scores)
			if probability < .5:
				total += 1
			else:
				correct += 1
				total += 1
	print "L5 Accuracy: %.5f " % (correct/total)
	
	#L2
	correct = 0
	total = 0
	for i in range(len(rr_list)):
		for j in range(int(len(pssm_list[i]) / 2)):
			prediction =  - 1
			feature = pssm_list[i][rr_list[i][j][0]-1] + pssm_list[i][rr_list[i][j][1]-1]
			scores = get_scores(feature, weights)
			probability = calculate_sigmoid(scores)
			if probability < .5:
				total += 1
			else:
				correct += 1
				total += 1
	print "L2 Accuracy: %.5f " % (correct/total)


if  __name__ == "__main__":
	main()

