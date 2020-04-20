from __future__ import division
import sys
from os import listdir
from os import path
from os.path import isfile, join
from math import log, exp 
import pickle
import re
import random
from id3 import tree_traversal
from GNB import GNB
TRAINING_PROPORTION = .75

def main():
	mode, fasta_path, pssm_path, tm_path = command_line_parser()
		
	if mode == "train":
		print "Running Training mode"
		protein_features = get_protein_features(mode, fasta_path, pssm_path)	
		features, labels = feature_label_generation(mode, protein_features, pssm_path, tm_path)
		weights = stochastic_linear_regression(features, labels, 100000000, .00001)	
		with open("weights.bin", "wb") as f:
			pickle.dump(weights, f)
		print "Finished training. Data stored in train.bin"
	else:
		print "Running Test mode"
		with open("weights.bin", "rb") as f:
			weights = pickle.load(f)
		protein_features = get_protein_features(mode, fasta_path, pssm_path)
		features, labels = feature_label_generation(mode, protein_features, pssm_path, tm_path)
		calculate_accuracy(weights, features, labels)


def calculate_accuracy(weights, features, labels):
	average_squared_error = 0

	for i in range(len(features)):
		prediction = 0
		for j in range(len(features[i])):
			prediction += features[i][j] * weights[j]
		prediction += weights[-1]
		average_squared_error += (labels[i] - prediction) ** 2

	average_squared_error /= len(features)
	print average_squared_error
	return average_squared_error


def get_protein_features(mode, fasta_path, pssm_path):
	pssm_proportions = get_pssm_proportions(pssm_path)
	EB_proportions = tree_traversal(fasta_path)
	HEC_proportions = GNB(pssm_path)
	
	protein_features = [ 0 for x in range(len(pssm_proportions))]
	for i in range(len(pssm_proportions)):
		protein_features[i] = pssm_proportions[i] + EB_proportions[i] + HEC_proportions[i]
	return protein_features


def get_pssm_proportions(pssm_path):
	files = [f for f in listdir(pssm_path) if isfile(join(pssm_path, f))]
	files = sorted(files)
	
	start = 0
	end = len(files)

	pssm_list = []
	pssm_proportions = []
	#iterate over all applicable pssm files
	for i in range(start, end):
		fp = open("%s/%s" % (pssm_path, files[i]), "r")
		lines = fp.readlines()
		lines = lines[3:-6]
		
		#get the raw contents (the second set of 20 features) of each file
		for j in range(len(lines)):
			lines[j].strip()
			all_numbers = map(int, re.findall('-?\d+', lines[j]))
			lines[j] = all_numbers[21:41]

		#average all values and then divide by 100 to get average proportion per file
		pssm_proportion = [0 for x in range(len(lines[0]))] 
		for k in range(len(lines[0])):

			for m in range(len(lines)):
				pssm_proportion[k] += lines[m][k]
		for k in range(len(lines[0])):
			pssm_proportion[k] /= len(lines)
			pssm_proportion[k] /= 100
	
		#print pssm_proportion
		pssm_proportions.append(pssm_proportion)
	return pssm_proportions


#Transforms pssm list to contain a sliding window of the previous and post 2 sequences, -1 values are included for edge cases
def command_line_parser():
	if len(sys.argv) == 1:
		mode = "train"
	else:
		mode = sys.argv[1]
		if mode != "train" and mode != "test" and len(sys.argv) == 2:
			print "Invalid mode: please use either \"train\" or \"test\""
			exit()
	
	if len(sys.argv) == 4:
		fasta = sys.argv[2]
		pssm_path = sys.argv[3]
		tm_path = sys.argv[4]
	else:
		fasta_path = "fasta"	
		pssm_path = "pssm"
		tm_path = "tmalign"
	return mode, fasta_path, pssm_path, tm_path


def stochastic_linear_regression(features, labels, steps, step_size):
	print "Running stochastic gradient descent with %d steps and step size %.10f" % (steps, step_size)
	weights = [0 for x in range(51)]
	
	for step in xrange(steps):
		gradients = [0 for x in range(51)]
		random_seed = random.randrange(len(features))
		feature = features[random_seed]
		label = labels[random_seed]		
		misassignment = 0
		for i in range(len(weights) - 1):
			misassignment += weights[i] * feature[i]
		misassignment = (label - (weights[-1] + misassignment))

		for i in range(len(gradients) - 1):
			gradients[i] = -2 * feature[i] * misassignment

		#get value for weight 0, assum x0 = 1
		weights[-1] = weights[-1] + 2 * step_size * misassignment
		
		#get values for rest of weights
		for i in range(len(weights) - 1):
			weights[i] = weights[i] - step_size * gradients[i]

		if step % 10000 == 0:
			print "Iteration %d" % (step)

	return weights
def feature_label_generation(mode, protein_features, pssm_path, tm_path):
	files = [path.splitext(f)[0] for f in listdir(pssm_path) if isfile(join(pssm_path, f))]
	files = sorted(files)
	protein_pair_features = []
	protein_pair_labels = []

	count = 0
	for i in range(len(files)):
		for j in range(i + 1, len(files)):
			#get feature set for each pair
			pair = protein_features[i] + protein_features[j]
			protein_pair_features.append(pair)

			#get tm scores for each unique pair
			tm_file = files[i] + "_" + files[j] + "_tmalign"
			fp = open("%s/%s" % (tm_path, tm_file), "r")
			lines = fp.readlines()
			score_chain_1 = map(float, re.findall('\d*\.\d+', lines[17]))
			score_chain_2 = map(float, re.findall('\d*\.\d+', lines[18]))
			average_score = (score_chain_1[0] + score_chain_2[0]) / 2
			protein_pair_labels.append(average_score)

	if mode == "train":
		start = 0
		end = int(TRAINING_PROPORTION*len(protein_pair_features))

	else:
		start = int(TRAINING_PROPORTION*len(protein_pair_features))
		end = len(protein_pair_features)

	return protein_pair_features[start:end], protein_pair_labels[start:end]

if  __name__ == "__main__":
	main()

