from __future__ import division
import sys
from os import listdir
from os.path import isfile, join
from math import sqrt, pi, exp 
import pickle
import re

TRAINING_PROPORTION = .75

def main():
	mode, fasta_path, sa_path, pssm_path = command_line_parser()
		
	if mode == "train":
		print "Running Training mode"
		fasta_list, sa_list, pssm_list = read_files(mode, fasta_path, sa_path, pssm_path)
		pssm_list = pssm_transform(pssm_list)

		class_probabilities = calculate_class_probabilities(sa_list)
		mean, var = calculate_mean(sa_list, pssm_list)
		data = [class_probabilities, mean, var]
		with open("train.bin", "wb") as f:
			pickle.dump(data, f)		
		print "Finished training. Data stored in train.bin"

	else:
		print "Running Test mode"

		with open("train.bin", "rb") as f:
			data = pickle.load(f)
		class_probabilities = data[0]
		mean = data[1]
		var = data[2]
		fasta_list, sa_list, pssm_list = read_files(mode, fasta_path, sa_path, pssm_path)
		pssm_list = pssm_transform(pssm_list)
		prediction = GNB(pssm_list, mean, var, class_probabilities)

		q3(sa_list, prediction)

def q3(sa_list, prediction):
	num_correct = 0
	num_total = 0

	for i in range(len(sa_list)):
		for j in range(len(sa_list[i])):
			if sa_list[i][j] == "H" and prediction[i][j] == 0:
				num_correct += 1

			elif sa_list[i][j] == "E" and prediction[i][j] == 1:
				num_correct += 1
			elif sa_list[i][j] == "C" and prediction[i][j] == 2:
				num_correct += 1
			num_total += 1

	print "q3 accuracy = %f" % (num_correct/num_total)


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
	
	if len(sys.argv) == 5:
		fasta_path = sys.argv[2]
		sa_path = sys.argv[3]
		pssm_path = sys.argv[4]
	else:
		fasta_path = "fasta"
		sa_path = "ss"	
		pssm_path = "pssm"

	return mode, fasta_path, sa_path, pssm_path


def read_files(mode, fasta_path, sa_path, pssm_path):
	sequences_list = []
	fasta_sequence = ""	

	files = [f for f in listdir(fasta_path) if isfile(join(fasta_path, f))]
	files = sorted(files)
	
	if mode == "train":
		start = 0
		end = int(TRAINING_PROPORTION*len(files))
	else:
		start = int(TRAINING_PROPORTION*len(files))
		end = len(files)
	
	#fasta intake
	for i in range(start, end):
		fp = open("%s/%s" % (fasta_path, files[i]), "r")
		sequence_as_list = fp.readlines()
		fasta_sequence = ""
		for line in range(1, len(sequence_as_list)):
			fasta_sequence += sequence_as_list[line].strip()
		sequences_list.append(fasta_sequence)
		fp.close()
	
	#sa intake
	sa_list = []
	sa_sequence = ""	
	files = [f for f in listdir(sa_path) if isfile(join(sa_path, f))]
	files = sorted(files)

	for i in range(start, end):
		fp = open("%s/%s" % (sa_path, files[i]), "r")
		sequence_as_list = fp.readlines()
		sa_sequence = ""
		for line in range(1, len(sequence_as_list)):
			sa_sequence += sequence_as_list[line].strip()
		sa_list.append(sa_sequence)
		fp.close()


	#pssm intake
	pssm_list = []
	files = [f for f in listdir(pssm_path) if isfile(join(pssm_path, f))]
	files = sorted(files)

	for i in range(start, end):
		fp = open("%s/%s" % (pssm_path, files[i]), "r")
		lines = fp.readlines()
		lines = lines[3:-6]
		for j in range(len(lines)):
			lines[j].strip()
			all_numbers = map(int, re.findall('-?\d+', lines[j]))
			lines[j] = all_numbers[1:21]

		pssm_list.append(lines)
	return sequences_list, sa_list, pssm_list 


def calculate_class_probabilities(sa_list):
	#in H, E, C order
	class_probabilities = [0, 0, 0]
	
	for i in range(len(sa_list)):
		for j in range(len(sa_list[i])):
			if sa_list[i][j] == "H":
				class_probabilities[0] += 1
			elif sa_list[i][j] == "E":
				class_probabilities[1] += 1
			elif sa_list[i][j] == "C":
				class_probabilities[2] += 1
			else:
				print "Unknown shape in ss file. Exiting..."
				quit()

	total = class_probabilities[0] + class_probabilities[1] + class_probabilities[2]
	class_probabilities[0] /= total
	class_probabilities[1] /= total
	class_probabilities[2] /= total


	return class_probabilities

def calculate_mean(sa_list, pssm_list):
	
	#conditional mean and var listed in H, E, C order
	mean = [[0 for col in range(3)] for row in range(100)]
	var = [[0 for col in range(3)] for row in range(100)]
	sa_map = {"H":0, "E":1, "C":2}

	total = [0, 0, 0]
	for i in range(len(pssm_list)):
		for j in range(len(pssm_list[i])):	
			for k in range(len(pssm_list[i][j])):
				mean[k][sa_map[sa_list[i][j]]] += pssm_list[i][j][k]
			
			total[sa_map[sa_list[i][j]]] += 1
					
	for k in range(len(mean)):
		mean[k][0] /= total[0]
		mean[k][1] /= total[1]
		mean[k][2] /= total[2]
	

	for i in range(len(pssm_list)):
		for j in range(len(pssm_list[i])):
			for k in range(len(pssm_list[i][j])):
				var[k][sa_map[sa_list[i][j]]] += (pssm_list[i][j][k] - mean[k][sa_map[sa_list[i][j]]]) ** 2
	
			
	for k in range(len(mean)):
		var[k][0] = var[k][0] / total[0]
		var[k][1] = var[k][1] / total[1]
		var[k][2] = var[k][2] / total[2]
		
	return mean, var


def calculate_pdf(x, mean, var):
	pdf = [1, 1, 1]

	for i in range(0, 3):
		for j in range(len(mean)):
			pdf[i] *= (1 / sqrt(2 * pi * var[j][i])) * exp(-0.5 * ((x[j] - mean[j][i]) ** 2) / var[j][i])
	return pdf


def GNB(pssm_list, mean, var, class_prob):
	expected_sa_list = []
	prediction = []	 

	for i in range(len(pssm_list)):
		file_prediction = []	
		for j in range(len(pssm_list[i])):
			conditional_prob = [1, 1, 1]
			pdf = calculate_pdf(pssm_list[i][j], mean, var)	
			total_prob = 0
		
			for k in range(0,3):
				total_prob += (pdf[k] * class_prob[k])
			for k in range(0,3):	
				conditional_prob[k] = (pdf[k] * class_prob[k])/total_prob
			
			file_prediction.append(conditional_prob.index(max(conditional_prob))) 

		prediction.append(file_prediction)

	return prediction


if  __name__ == "__main__":
	main()

