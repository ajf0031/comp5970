from __future__ import division
import sys
from os import listdir
from os.path import isfile, join
from math import sqrt, pi, exp 
import pickle
import re


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


def read_files(pssm_path):
	pssm_list = []
	files = [f for f in listdir(pssm_path) if isfile(join(pssm_path, f))]
	files = sorted(files)
	start = 0
	end = len(files)

	for i in range(start, end):
		fp = open("%s/%s" % (pssm_path, files[i]), "r")
		lines = fp.readlines()
		lines = lines[3:-6]
		for j in range(len(lines)):
			lines[j].strip()
			all_numbers = map(int, re.findall('-?\d+', lines[j]))
			lines[j] = all_numbers[1:21]

		pssm_list.append(lines)
	return pssm_list 


def calculate_pdf(x, mean, var):
	pdf = [1, 1, 1]

	for i in range(0, 3):
		for j in range(len(mean)):
			pdf[i] *= (1 / sqrt(2 * pi * var[j][i])) * exp(-0.5 * ((x[j] - mean[j][i]) ** 2) / var[j][i])
	return pdf


def GNB(pssm_path):
	pssm_list = read_files(pssm_path)
	pssm_list = pssm_transform(pssm_list)
	with open("GNB.bin", "rb") as f:
		data = pickle.load(f)
	class_prob = data[0]
	mean = data[1]
	var = data[2]
	
	proportion_HEC = []
	for i in range(len(pssm_list)):
		for j in range(len(pssm_list[i])):
			conditional_prob = [1, 1, 1]
			pdf = calculate_pdf(pssm_list[i][j], mean, var)	
			total_prob = 0
		
			for k in range(0,3):
				total_prob += (pdf[k] * class_prob[k])
			for k in range(0,3):	
				conditional_prob[k] = (pdf[k] * class_prob[k])/total_prob
			
		proportion_HEC.append(conditional_prob)
		
		

	return proportion_HEC


