from __future__ import division
import sys
from os import listdir
from os.path import isfile, join
from math import log
import attributes
from node import Node
import pickle

TRAINING_PROPORTION = .75

def main():
	mode, fasta_path, sa_path= command_line_parser()
		
	if mode == "train":
		print "Running Training mode"
		fasta_list, sa_list = read_sequences(mode, fasta_path, sa_path)
		buried_exposed = calculate_buried_exposed(fasta_list, sa_list)
		attributes_in = attributes.attributes
		root = id3(buried_exposed, attributes_in)
		with open("tree.bin", 'wb') as f:	
			pickle.dump(root, f)
		print "Finished training, decision tree stored in tree.bin"

	else:
		print "Running Test mode"
		with open("tree.bin", "rb") as f:
			root = pickle.load(f)
		start = 112
		end = 150
		fasta_list, sa_list = read_sequences(mode, fasta_path, sa_path)
		actual_buried_exposed = calculate_buried_exposed(fasta_list, sa_list)
		expected_buried_exposed = tree_traversal(root, actual_buried_exposed)
		calculate_statistics(actual_buried_exposed, expected_buried_exposed)
		
def command_line_parser():
	if len(sys.argv) != 2:
		mode = "train"
	else:
		mode = sys.argv[1]
		if mode != "train" and mode != "test" and len(sys.argv) == 2:
			print "Invalid mode: please use either \"train\" or \"test\""
			exit()
	
	if len(sys.argv) == 3:
		fasta_path = sys.argv[1]
		sa_path = sys.argv[2]
	else:
		fasta_path = "fasta"
		sa_path = "sa"	

	return mode, fasta_path, sa_path



def tree_traversal(root, actual_buried_exposed):
	expected_buried_exposed = {}
	attributes_dict = attributes.attributes
	
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
				if amino_acid in attributes_dict[iterator.data]:
					iterator = iterator.right
				else:
					iterator = iterator.left
	
	return expected_buried_exposed	


def calculate_statistics(actual_buried_exposed, expected_buried_exposed):
	true_positive, true_negative, false_positive, false_negative = (0, 0, 0, 0)

	for amino_acid in actual_buried_exposed:
		exposed =  actual_buried_exposed[amino_acid][0]
		buried = actual_buried_exposed[amino_acid][1]
		if expected_buried_exposed[amino_acid] == "buried":
			true_negative += actual_buried_exposed[amino_acid][1]
			false_negative += actual_buried_exposed[amino_acid][0] 
		else:
			true_positive += actual_buried_exposed[amino_acid][0]
			false_positive += actual_buried_exposed[amino_acid][1] 

	precision = (true_positive / (true_positive + false_positive))
	recall = (true_positive / (true_positive + false_negative))
	accuracy = ((true_positive + true_negative)/(true_positive + true_negative + false_positive + false_negative))
	f1 = 2*(precision*recall)/(precision + recall)

	print "Precision: %f" % precision
 	print "Recall: %f" % recall
	print "Accuracy: %f" % accuracy
	print "F1: %f" % f1


def read_sequences(mode, fasta_path, sa_path):
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

	for i in range(start, end):
		fp = open("%s/%s" % (fasta_path, files[i]), "r")
		sequence_as_list = fp.readlines()
		fasta_sequence = ""
		for line in range(1, len(sequence_as_list)):
			fasta_sequence += sequence_as_list[line].strip()
		sequences_list.append(fasta_sequence)
		fp.close()

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
		amino_acid_sa_tuple[value] = (amino_acid_exposed[value], \
					amino_acid_total[value]- amino_acid_exposed[value])

	return amino_acid_sa_tuple


# Recursive function that partitions the amino acids based on the attribute on the highest gain, 
# and removes that selected attribute from the next iteration. Takes in a dictionary of amino acids
# with [exposed, buried] and a dictionary of attribute dictionaries.
def id3(amino_acid_tuple, attributes_in):
	attributes = attributes_in
	cur_node = Node()
	
	total = [0, 0] # exposed, buried
	for amino_acid in amino_acid_tuple:
		total[0] += amino_acid_tuple[amino_acid][0]
		total[1] += amino_acid_tuple[amino_acid][1]
	
	#Edge cases- if there are no exposed or no buried, then it is a perfect match and that is a leaf node
	#If there no more attributes remaining, then gain cannot be run
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
	attributes_gain = calculate_gain(total_entropy, attributes, amino_acid_tuple, total)
	
	# If the only attributes remaining have no information gain, then we pick whichever sa count is higher
	if len(attributes_gain) == 0:
		if total[0] > total[1]:
			cur_node.data = "exposed"
		else:
			cur_node.data = "buried"
		return cur_node
	
	best_attribute =  max(attributes_gain, key=attributes_gain.get)

	# These dictionaries partition the amino acids into those classified by the presence or absnece
	# of an attribute	
	not_attribute_amino_acid_tuple = {}
	attribute_amino_acid_tuple = {}

	for amino_acid in amino_acid_tuple:
		if amino_acid in attributes[best_attribute]:
			attribute_amino_acid_tuple[amino_acid] = amino_acid_tuple[amino_acid]
		else:
			not_attribute_amino_acid_tuple[amino_acid] = amino_acid_tuple[amino_acid]

	cur_node.data = best_attribute
	element = attributes.pop(best_attribute)

	#Break up the tree into the left and right partition
	cur_node.left = id3(not_attribute_amino_acid_tuple, attributes)
	cur_node.right = id3(attribute_amino_acid_tuple, attributes)

	attributes[best_attribute] = element
	
	return cur_node


def calculate_entropy(exposed_count, buried_count):
	probability_exposed = exposed_count / (exposed_count + buried_count)
	probability_buried = buried_count / (exposed_count + buried_count)
	entropy = -probability_exposed*log(probability_exposed, 2) - probability_buried*log(probability_buried, 2)
	return entropy


def calculate_gain(total_entropy, attributes, amino_acid_tuple, total):
	attributes_gain = {}
	for attribute in attributes:
		attribute_total = [0, 0]	
	
		for amino_acid in amino_acid_tuple:
			if amino_acid in attributes[attribute]: 
				attribute_total[0] += amino_acid_tuple[amino_acid][0]					
				attribute_total[1] += amino_acid_tuple[amino_acid][1]
		
		#If the attribute explains none or all the data, it is useless	
		if attribute_total == [0, 0] or attribute_total == total:
			continue 

		entropy_of_attribute = (attribute_total[0] + attribute_total[1])/ \
					(total[0] + total[1]) * calculate_entropy( \
					attribute_total[0], attribute_total[1])

		entropy_of_not_attribute = (total[0] + total[1] - attribute_total[0] - attribute_total[1])/ \
					(total[0] + total[1])	* calculate_entropy( \
					total[0]-attribute_total[0], total[1] - attribute_total[1])
		
		attributes_gain[attribute] = total_entropy - entropy_of_attribute - entropy_of_not_attribute
	
	return attributes_gain	


if  __name__ == "__main__":
	main()

