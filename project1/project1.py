import sys
 
BLOSUM_MATRIX = []
BLOSUM_MAP = {} #maps the blosum amino acid characters to integers
MAX_CHARS_PER_LINE = 80 #max number of characters printed per line

def main():
	read_blosum()
	seq1, seq2 = read_sequences()
	needleman_wunsch(seq1, seq2)
	smith_waterman(seq1, seq2)


def read_blosum():
	try:
		blosum_file = open("blosum.txt", "r")
	except: 
		exit("Error opening blosum file, there should be a file called blosum.txt in this directory")
	index = 0
	blosum_file.readline()
	for line in blosum_file: 
		char = line.split()
		BLOSUM_MAP[char[0]] = index
		index += 1
		BLOSUM_MATRIX.append([int(x) for x in char[1:]])
	
	blosum_file.close()


def read_sequences():
	if len(sys.argv) == 2:
		if sys.argv[1] == "pair1":
			file1 = "pair1//1k4rA_dengue_virus.fasta"
			file2 = "pair1//5ire_zika_virus.fasta"
		elif sys.argv[1] == 'pair2':
			file1 = "pair2//1dzl_HPV.fasta"
			file2 = "pair2//3mge_HIV.fasta"
		elif sys.argv[1] == "pair3":
			file1 = "pair3//2plv1_polio_virus.fasta"
			file2 = "pair3//4rhv1_rhino_virus.fasta"
		else:
			exit("Please format using python project1.py [pair1|pair2|pair3]")
	elif len(sys.argv) == 3:
		file1 = sys.argv[1]
		flie2 = sys.argv[2]
	else:
		file1 = "pair1//1k4rA_dengue_virus.fasta"
		file2 = "pair1//5ire_zika_virus.fasta"

	try:	
		fp1 = open(file1, "r")
		fp2 = open(file2, "r")
	except:
		exit("Error opening file. Please make sure to include the proper path of the file")
	
	sequence1_list = fp1.readlines()
	sequence2_list = fp2.readlines()
	sequence1 = "" 
	sequence2 = ""

	for line in range(1,len(sequence1_list)):
		sequence1 += sequence1_list[line].strip()
	for line in range(1,len(sequence2_list)):
		sequence2 += sequence2_list[line].strip()

	fp1.close()
	fp2.close()	
		
	return sequence1, sequence2


def needleman_wunsch(seq1, seq2):	
	value_map = { 0: "diagonal", 1: "left", 2: "up"}
	
	#BASE CASE for value_matrix and direction_matrix
	value_matrix = [[0 for x in range(len(seq1)+1)] for y in range(len(seq2)+1)]
	direction_matrix = [["" for x in range(len(seq1)+1)] for y in range(len(seq2)+1)]
	
	for i in range(1, len(seq1) + 1):	
		blosum_index = BLOSUM_MAP[seq1[i-1]] #extracts the blosum index from the sequence	
		value_matrix[0][i] = value_matrix[0][i-1] + BLOSUM_MATRIX[-1][blosum_index]
		direction_matrix[0][i] = "left"	
	for j in range(1, len(seq2) + 1):
		blosum_index = BLOSUM_MAP[seq2[j-1]]
		value_matrix[j][0] = value_matrix[j-1][0] + BLOSUM_MATRIX[blosum_index][-1]
		direction_matrix[j][0] = "up"
	
	#general case
	for i in range(1, len(seq1) + 1):
		blosum_index_i = BLOSUM_MAP[seq1[i-1]]
		for j in range(1, len(seq2) + 1):
			blosum_index_j = BLOSUM_MAP[seq2[j-1]]
			
			#find the 3 potential values that V(i,j) can be
			diagonal = BLOSUM_MATRIX[blosum_index_j][blosum_index_i] + value_matrix[j-1][i-1]
			vertical = BLOSUM_MATRIX[blosum_index_j][-1] + value_matrix[j-1][i]
			horizontal = BLOSUM_MATRIX[-1][blosum_index_i] + value_matrix[j][i-1]

			value_choice = max([diagonal, horizontal, vertical])
			value_matrix[j][i] = value_choice
			
			#put in the direction of the max score- Note for duplicates, diagonal is chosen first
			direction_matrix[j][i] = value_map[[diagonal, horizontal, vertical].index(value_choice)]

	s1, match, s2 = backtrack(seq1, seq2, len(seq1), len(seq2), value_matrix, direction_matrix)
	
	print("Needleman Wunsch Global Alignment")	
	print("Score: %s" % value_matrix[-1][-1])
	print_alignment(s1, match, s2)


def smith_waterman(seq1, seq2):
	max_value = 0
	max_value_i = 0
	max_value_j = 0
	value_map = { 0: "diagonal", 1: "left", 2: "up"}
	
	#BASE CASE for value_matrix
	# Note we keep the default values of 0 for the Value matrix 
	# and "" for the direction matrix because local alignment ends at 0
	value_matrix = [[0 for x in range(len(seq1)+1)] for y in range(len(seq2)+1)]	
	direction_matrix = [["" for x in range(len(seq1)+1)] for y in range(len(seq2)+1)]
	
	#general case
	for i in range(1, len(seq1) + 1):
		blosum_index_i = BLOSUM_MAP[seq1[i-1]]
		for j in range(1, len(seq2) + 1):
			blosum_index_j = BLOSUM_MAP[seq2[j-1]]
			#find the 3 potential values that V(i,j) can be. The minimum they can be set to is 0.
			diagonal = max(BLOSUM_MATRIX[blosum_index_j][blosum_index_i] + value_matrix[j-1][i-1], 0)
			vertical = max(BLOSUM_MATRIX[blosum_index_j][-1] + value_matrix[j-1][i], 0)
			horizontal = max(BLOSUM_MATRIX[-1][blosum_index_i] + value_matrix[j][i-1], 0)

			value_choice = max([diagonal, horizontal, vertical])
			if value_choice > max_value:
				max_value_i = i
				max_value_j = j
				max_value = value_choice		
			value_matrix[j][i] = value_choice
			
			#put in the direction of the max score, but only if if is greater than 0, 
			# otherwise there is no direction
			if value_choice > 0:
				direction_matrix[j][i] = value_map[[diagonal, horizontal, vertical].index(value_choice)]

	s1, match, s2 = backtrack(seq1, seq2, max_value_i, max_value_j, value_matrix, direction_matrix)
	
	print("Smith Waterman Local Alignment")
	print("Score: %s" % max_value)
	print_alignment(s1, match, s2)


def backtrack(seq1, seq2, start_x, start_y, value_matrix, direction_matrix):
	#s1 and s2 are the char arrays with gaps
	s1 = []
	s2 = []
	match = []
	x = start_x
	y = start_y
	while (direction_matrix[y][x] != ""):
		if direction_matrix[y][x] == "diagonal":
			s1.append(seq1[x-1])
			s2.append(seq2[y-1])
			if value_matrix[y][x] > value_matrix[y-1][x-1]:
				match.append("|")
			else:
				match.append("*")
			x -= 1
			y -= 1
		elif direction_matrix[y][x] == "left":
			s1.append(seq1[x-1])
			s2.append("-")
			match.append(" ")
			x -= 1
		else:
			s1.append("-")
			s2.append(seq2[y-1])
			match.append(" ")
			y -= 1
	
	s1.reverse()
	match.reverse()
	s2.reverse()

	return s1, match, s2


def print_alignment(sequence1, match_characters, sequence2):
	index = 0
	remaining_list_size = len(sequence1)
	while MAX_CHARS_PER_LINE < remaining_list_size:
		print (''.join(sequence1[index:index+MAX_CHARS_PER_LINE]))
		print (''.join(match_characters[index:index+MAX_CHARS_PER_LINE]))
		print (''.join(sequence2[index:index+MAX_CHARS_PER_LINE]))
		print("")
		index += MAX_CHARS_PER_LINE
		remaining_list_size -= MAX_CHARS_PER_LINE
	print (''.join(sequence1[index:]))
	print (''.join(match_characters[index:]))
	print (''.join(sequence2[index:]))
	print ("")


if __name__ == "__main__":
	main()

