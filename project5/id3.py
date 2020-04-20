from __future__ import division
import sys
from os import listdir
from os.path import isfile, join
import attributes
from node import Node
import pickle

def read_sequences(fasta_path):
	sequences_list = []
	fasta_sequence = ""	

	files = [f for f in listdir(fasta_path) if isfile(join(fasta_path, f))]
	files = sorted(files)

	start = 0
	end = len(files)

	for i in range(start, end):
		fp = open("%s/%s" % (fasta_path, files[i]), "r")
		sequence_as_list = fp.readlines()
		fasta_sequence = ""
		for line in range(1, len(sequence_as_list)):
			fasta_sequence += sequence_as_list[line].strip()
		sequences_list.append(fasta_sequence)
		fp.close()
	return sequences_list

def tree_traversal(fasta_path):
	fasta_list = read_sequences(fasta_path)
	

	buried_exposed_list = []
	attributes_dict = attributes.attributes
	with open("tree.bin", "rb") as f:
		root = pickle.load(f)
	for fasta_sequence in fasta_list:
		proportion_buried_exposed = [0, 0]
		
		for amino_acid in fasta_sequence:
			iterator = root
			while iterator != None:
				if iterator.data == "buried":
					proportion_buried_exposed[0] += 1
					iterator = None
				elif iterator.data == "exposed":
					proportion_buried_exposed[1] += 1
					iterator = None
				else:
					if amino_acid in attributes_dict[iterator.data]:
						iterator = iterator.right
					else:
						iterator = iterator.left

		proportion_buried_exposed[0] /= len(fasta_sequence)
		proportion_buried_exposed[1] /= len(fasta_sequence)
		buried_exposed_list.append(proportion_buried_exposed)	
	return buried_exposed_list


