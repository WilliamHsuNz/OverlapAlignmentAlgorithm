# OverlapAlignmentAlgorithm
This programme implements an overlap alignment algorithm to align two molecular sequences. 
The programme takes as input two sequences, a score matrix and a gap penalty.  The programme
constructs a score table and the score table is backtracked to find the best alignment.  

Instructions 
Use the following commandline to run the programme.  -m is for match score -M is for mismatch
score and -d is for gap penalty.  

python -m Problem3a -m 2 -M -2 -d -3 Problem2c_descendant1.csv Problem2c_descendant2.csv
