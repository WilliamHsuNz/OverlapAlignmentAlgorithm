import numpy as np
'''
Author: William Hsu
Date: 11 May 2015

This programme performs alignment of two sequences using the overlap alignment algorithm.  
The programme takes a score matrix and a gap penalty and constructs a score table.  The 
score table is backtacked to find the best alignment.  

Instructions 
Use the following commandline to run the programme.  -m is for match score -M is for mismatch
score and -d is for gap penalty.  

python -m Problem3a -m 2 -M -2 -d -3 Problem2c_descendant1.csv Problem2c_descendant2.csv


for 3 (c)
python -m Problem3a -m 2 -M -2 -d -4 Problem2c_descendant1.csv Problem2c_descendant2.csv
python -m Problem3a -m 2 -M -2 -d -2 Problem2c_descendant1.csv Problem2c_descendant2.csv
python -m Problem3a -m 2 -M -2 -d -1 Problem2c_descendant1.csv Problem2c_descendant2.csv
'''
def overlap_alignment(sequence_1, sequence_2, score_matrix, gapPenalty):
	
	cols = len(sequence_1)+1
	rows = len(sequence_2)+1 
	
	geneDict = {"A":0, "C":1, "G":2, "T":3}
				
	#creating score table			
	score_table = [[0 for x in range(cols)] for x in range(rows)]
	for i in range(rows):
		score_table[i][0] = 0 #gapPenalty * i
	for j in range(cols):
		score_table[0][j] = 0 #gapPenalty * j	
	
	for i in range(1,rows):
		for j in range(1, cols):
			match = score_table[i-1][j-1] + score_matrix[geneDict[sequence_2[i-1]]][geneDict[sequence_1[j-1]]]
			delete = score_table[i-1][j] + gapPenalty	
			insert = score_table[i][j-1] + gapPenalty
			score_table[i][j] = max(match, delete, insert) 	
	 		j = j +1
		i = i+1
		
	#print score_table
	print "Score Table"
	for i in range(rows):
		for j in range(cols):
			print "%4d" %score_table[i][j], 
		print

	max_row_index, max_col_index = find_position_to_backtrack(score_table)
	
	#backtrack and form alignment
	AlignmentA = ""
	AlignmentB = ""

	i = max_row_index 
	j = max_col_index 
	while (i > 0 and j > 0): #changed to and
		if (i > 0 and j > 0 and (score_table[i][j] == score_table[i-1][j-1] + score_matrix[geneDict[sequence_2[i-1]]][geneDict[sequence_1[j-1]]])):
			AlignmentA = sequence_1[j-1] + AlignmentA
			AlignmentB = sequence_2[i-1] + AlignmentB
			i = i -1
			j = j -1
		elif (i > 0 and score_table[i][j] == score_table[i-1][j] + gapPenalty):
			AlignmentB = sequence_2[i-1] + AlignmentB
			AlignmentA = "-" + AlignmentA
			i = i -1
		elif (j > 0 and score_table[i][j] == score_table[i][j-1] + gapPenalty):
			AlignmentB = "-" + AlignmentB
			AlignmentA = sequence_1[j-1] + AlignmentA
			j = j - 1
	print AlignmentA
	print AlignmentB
	
# create a score matrix	
def create_score_matrix(m, M):
	score_matrix = [[0 for x in range(4)] for x in range(4)]
	for i in range(4):
		for j in range(4):
			if i == j:
				score_matrix [i][j] = m
			else: 
				score_matrix [i][j] = M	
	return score_matrix

#find position to backtrack from	
def find_position_to_backtrack(score_table):
	height = len(score_table)
	width = len(score_table[0])
	max = score_table[height-1][width-1]
	max_row_index = height-1
	max_col_index = width-1
	
	for i in range(height):
		if score_table[i][width-1] > max:
			max = score_table[i][width - 1]
			max_row_index = i
			max_col_index = width - 1
	for j in range(width):
		if score_table[height - 1][j] > max:
			max_row_index = height - 1
			max_col_index = j
	return max_row_index, max_col_index		

if __name__ == '__main__':
	from optparse import OptionParser
	import csv
	
	sequence_1 = []
	sequence_2 = []
	p = OptionParser(usage='%prog data_file')
	p.add_option('-m', '--match', dest='match', type='int')
	p.add_option('-M', '--Mismatch', dest='Mismatch', type='int')
	p.add_option('-d', '--gapPenalty', dest='gapPenalty', type='int')
	options, args = p.parse_args()
	if len(args) < 2:
		p.error('must provide the path to two CSV file to read')
			
	s1 = open(args[0])
	s2 = open(args[1])
	try:
		s1_iterable = csv.reader(s1, skipinitialspace=True, delimiter=',', quoting=csv.QUOTE_NONE)
		s2_iterable = csv.reader(s2, skipinitialspace=True, delimiter=',', quoting=csv.QUOTE_NONE)
		for item in s1_iterable:
			sequence_1.extend(item)
		for item in s2_iterable:
			sequence_2.extend(item)	
		sequence_1 = sequence_1[:35]
		sequence_2 = sequence_2[len(sequence_2)-35:len(sequence_2)]
		
		score_matrix = create_score_matrix(options.match, options.Mismatch)
		overlap_alignment(sequence_1, sequence_2, score_matrix, options.gapPenalty)
		
	finally:
		s1.close()
		s2.close()	 
