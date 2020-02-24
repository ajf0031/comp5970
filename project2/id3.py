from __future__ import division
import sys
from os import listdir
from os.path import isfile, join
from math import log
import get_attributes
from node import Node
import pickle

def main():
	if sys.argv[1] == "train":
		start = 0
		end = 112
		fasta_list, sa_list = read_sequences(start, end)
		buried_exposed = calculate_buried_exposed(fasta_list, sa_list)
		attributes_in = get_attributes.attributes
		root = choose_decision(buried_exposed, attributes_in)

		with open("tree.bin", 'wb') as f:	
			pickle.dump(root, f)
	elif sys.argv[1] == "test":
		expected_buried_exposed = {}
		with open("tree.bin", "rb") as f:
			root = pickle.load(f)
		start = 112
		end = 150
		fasta_list, sa_list = read_sequences(start, end)
		actual_buried_exposed = calculate_buried_exposed(fasta_list, sa_list)
		attributes = get_attributes.attributes
		for amino_acid in actual_buried_exposed:
			iterator = root
			while iterator != None:
				if iterator.data == "buried":
					expected_buried_exposed[amino_acid] = "buried"
					iterator = None
				elif iterator.data == "exposed":
					expected_buried_exposed[amino_acid] = "exposed"
					iterator = None
				else:
					if amino_acid in attributes[iterator.data]:
						iterator = iterator.right
					else:
						iterator = iterator.left
		print expected_buried_exposed	
		true_positive = 0
		true_negative = 0
		false_positive = 0
		false_negative = 0
		for amino_acid in actual_buried_exposed:
			exposed =  actual_buried_exposed[amino_acid][0]
			buried = actual_buried_exposed[amino_acid][1]
			if expected_buried_exposed[amino_acid] == "buried":
				true_negative += actual_buried_exposed[amino_acid][1]
				false_negative += actual_buried_exposed[amino_acid][0] 
			else:
				true_positive += actual_buried_exposed[amino_acid][0]
				false_positive += actual_buried_exposed[amino_acid][1] 
		print [true_negative, false_negative, true_positive, false_positive]
		print "Precision: %f" % (true_positive / (true_positive + false_positive))
		print "Recall %f" % (true_positive / (true_positive + false_negative))
		print "Accuracy %f" % ((true_positive + true_negative)/(true_positive + true_negative + false_positive + false_negative))
	else:
		print "Improper use command line arguments. Please use train or test"

def read_sequences(start, end):
	fasta_path = "fasta/"
	sa_path = "sa/"

	sequences_list = []
	sequence = ""	
	files = [f for f in listdir(fasta_path) if isfile(join(fasta_path, f))]
	files = sorted(files)
	for i in range(start, end):
		fp = open("%s%s" % (fasta_path, files[i]), "r")
		sequence_as_list = fp.readlines()
		sequence = ""
		for line in range(1, len(sequence_as_list)):
			sequence += sequence_as_list[line].strip()
		sequences_list.append(sequence)
		fp.close()

	sa_list = []
	sequence = ""	
	files = [f for f in listdir(sa_path) if isfile(join(sa_path, f))]
	files = sorted(files)

	for i in range(start, end):
		fp = open("%s%s" % (sa_path, files[i]), "r")
		sequence_as_list = fp.readlines()
		sequence = ""
		for line in range(1, len(sequence_as_list)):
			sequence += sequence_as_list[line].strip()
		sa_list.append(sequence)
		fp.close()
	return sequences_list, sa_list 

#returns list of tuples of each amino acid in the form of (exposed_count, buried_count)
def calculate_buried_exposed(fasta_list, sa_list):
	amino_acid_exposed = {}
	amino_acid_total = {}
	amino_acid_sa_tuple = {}
	for i in range(len(fasta_list)):
		for j, v in enumerate(fasta_list[i]):	
			if v not in amino_acid_total:
				if sa_list[i][j] == "E":
					amino_acid_exposed[v] = 1
				else:
					amino_acid_exposed[v] = 0
				amino_acid_total[v] = 1
			else:
				if sa_list[i][j] == "E":
					amino_acid_exposed[v] += 1
				amino_acid_total[v] += 1

	for value in amino_acid_exposed:
		amino_acid_sa_tuple[value] = (amino_acid_exposed[value], amino_acid_total[value]- amino_acid_exposed[value])

	return amino_acid_sa_tuple


def calculate_entropy(exposed_count, buried_count):
	probability_exposed = exposed_count / (exposed_count + buried_count)
	probability_buried = buried_count / (exposed_count + buried_count)
	entropy = -probability_exposed*log(probability_exposed, 2) - probability_buried*log(probability_buried, 2)
	return entropy

def choose_decision(amino_acid_tuple, attributes_in):
	attributes = attributes_in
	cur_node = Node()
	attributes_gain = {}
	
	total = [0, 0] # exposed, buried
	for amino_acid in amino_acid_tuple:
		total[0] += amino_acid_tuple[amino_acid][0]
		total[1] += amino_acid_tuple[amino_acid][1]

	if total[0] == 0:
		cur_node.data = "buried"
		return cur_node 
	if total[1] == 0:
		cur_node.data = "exposed"
		return cur_node
	if len(attributes) == 0:
		if total[0] > total[1]:
			cur_node.data = "exposed"
		else:
			cur_node.data = "buried"
		return cur_node
	total_entropy = calculate_entropy(total[0], total[1])
	for attribute in attributes:
		attribute_total = [0, 0]
		
		for amino_acid in amino_acid_tuple:
			if amino_acid in attributes[attribute]: 
				attribute_total[0] += amino_acid_tuple[amino_acid][0]					
				attribute_total[1] += amino_acid_tuple[amino_acid][1]
				
		if attribute_total == [0, 0] or attribute_total == total:
			continue 
		conditional_entropy = (attribute_total[0] + attribute_total[1])/(total[0] + total[1]) * calculate_entropy(attribute_total[0], attribute_total[1])
		conditional_entropy_not = (total[0] + total[1] - attribute_total[0] - attribute_total[1])/(total[0] + total[1])	* calculate_entropy(total[0]-attribute_total[0], total[1] - attribute_total[1])
		attributes_gain[attribute] = total_entropy - conditional_entropy - conditional_entropy_not

	if len(attributes_gain) == 0:
		if total[0] > total[1]:
			cur_node.data = "exposed"
		else:
			cur_node.data = "buried"
		return cur_node
	best_attribute =  max(attributes_gain, key=attributes_gain.get)
	not_attribute_amino_acid_tuple = {}
	attribute_amino_acid_tuple = {}
	for amino_acid in amino_acid_tuple:
		if amino_acid in attributes[best_attribute]:
			attribute_amino_acid_tuple[amino_acid] = amino_acid_tuple[amino_acid]
		else:
			not_attribute_amino_acid_tuple[amino_acid] = amino_acid_tuple[amino_acid]
	
	
	print best_attribute	

	#print "attribute amino acid set: %s" % attribute_amino_acid_tuple
	#print "not attribute amino acid set: %s" % not_attribute_amino_acid_tuple

	element = attributes.pop(best_attribute)

	#print "remaining attributes: %s" % attributes	
	cur_node.data = best_attribute
	cur_node.left = choose_decision(not_attribute_amino_acid_tuple, attributes)

	print "Finished left, moving right"
	cur_node.right = choose_decision(attribute_amino_acid_tuple, attributes)
	attributes[best_attribute] = element
	return cur_node


if  __name__ == "__main__":
	main()

