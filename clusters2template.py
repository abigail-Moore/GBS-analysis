#! /usr/bin/env python

#This script goes through the smaller clusters of sequences from the second rtd_run.py run.
#(These are the clusters of full sequences that are made from each cluster of 30bp sequences.)
#The program chooses which of the smaller clusters to take (according to an input value),
#chooses the first sequence of each cluster, labels them ">largecluster_smallcluster" [sequence],
#and saves them to a fasta file. (So there will be one fasta file for each large cluster.)
#If there is only one cluster, then it takes the first sequence from that cluster.

import sys #We want to be able to output error messages to the screen and get input from the command line.
import gzip #We want to be able to open zipped files.

Usage = '''
clusters2template.py version 1.0, 24 Dec. 2013
This program analyzes results of the rtd_run.py analysis of each cluster of 
sequences, chooses which smaller clusters to analyze further, and saves one 
sequence from each of these smaller clusters to a fasta file.
clusters2template.py [lower limit of group size] [group name] [.clstats file] 
[.cluni file] [out file folder]
The .clstats and .cluni files can just be written "[path to that directory]/*.clstats", 
since the file names for those files are unreasonably long and difficult to remember.
'''

if len(sys.argv) <5:
	print Usage
else:
	IndLowLim = int(sys.argv[1])
	GrpName = sys.argv[2]
	InFileCLSt = sys.argv[3]
	InFileCLU = sys.argv[4]
	OutFilePre = sys.argv[5]

print("clusters2template.py")

OutFileName = OutFilePre+"template_"+GrpName+".fa"

'''
The .clstats file is in the following format:
0	12	8	0.089552238806

0: Cluster Number
12: number of sequences in the cluster
8: number of individuals in the cluster
0.089552238806: "dirt", a measure of sequence dissimilarity
'''

'''The .cluni file is in the following format (with sequences and phred scores truncated):
0	10.try1_998uniqed	GATCCTCATT	6	DDADBFHHHG	220_1	6	.	.
0: cluster number
10.try1_998uniqed: presumably sequence 10 in this file
GATCCTCATT: sequence
6: number of times it is present
DDADBFHHHG: phred score
220_1: individual(s) in which the sequence is present
6: number of times the sequence is present in each individual
.	.: something for the reverse sequences (but in ours they are already concatenated with the
forward sequences)
'''

#First to parse the .clstats file to get the group numbers for the clusters we want:
SubGrpList=[ ]
InFile = open(InFileCLSt, 'rU')
for Line in InFile:
	Line = Line.strip('\r').strip('\n').split('\t')
	if int(Line[2]) >= IndLowLim:
		SubGrpList.append(int(Line[0]))
InFile.close()

print("There are %d groups of sequences that are present in more than %d individuals." \
% (len(SubGrpList), IndLowLim))
print("These groups of sequences will be used for further analyses.")
sys.stderr.write("There are %d groups of sequences that are present in more than %d individuals." \
% (len(SubGrpList), IndLowLim))
sys.stderr.write("These groups of sequences will be used for further analyses.\n")


#Now to make a dictionary with one sequence per cluster:

SeqDict = { }
OldNum = 999
InFile = open(InFileCLU, 'rU')
for Line in InFile:
	Line = Line.strip('\r').strip('\n').split('\t')
	SubGrpNum = int(Line[0])
	if SubGrpNum != OldNum:
		SeqDict[SubGrpNum] = Line[2]
	OldNum = SubGrpNum
InFile.close()

print("There are %d groups of sequences in total.\n" % len(SeqDict))
sys.stderr.write("There are %d groups of sequences in total.\n" % len(SeqDict))

#Now to make a dictionary of the sequences we want:

TemplateDict = { }
if len(SeqDict) == 1: #If we have a .cluni file with the sequences in one group
	SeqName = ">"+GrpName+"_0"
	TemplateDict[SeqName] = SeqDict[SubGrpNum]
elif len(SeqDict) == 0: #If we don't have a .cluni file at all, just take the first sequence
	#from the filed of uniqued sequences.
	SeqName = ">"+GrpName+"_0"
	SeqFileName = OutFilePre+"_"+GrpName+"uniqued.gz"
	LineNum = 0
	InFile = gzip.open(SeqFileName)
	for Line in InFile:
		if LineNum < 1:
			Line = Line.strip('\r').strip('\n').split('\t')
			TemplateDict[SeqName] = Line[0]
			LineNum += 1
else: #If we have .cluni and .clstats files with sequences in multiple groups
	for SubGrpNum in SubGrpList:
		SeqName = ">"+GrpName+"_"+str(SubGrpNum)
		TemplateDict[SeqName] = SeqDict[SubGrpNum]

#Now to write this dictionary to a file:

OutFile = open(OutFileName,'w')
for SeqName in TemplateDict.keys():
	OutFile.write("%s\n%s\n" % (SeqName, TemplateDict[SeqName]))
OutFile.close()
	

print("%d template sequences (one per smaller cluster) were written to the file %s.\n" \
	% (len(TemplateDict),OutFileName))
sys.stderr.write("%d template sequences (one per smaller cluster) were written to the file %s.\n" \
	% (len(TemplateDict),OutFileName))
