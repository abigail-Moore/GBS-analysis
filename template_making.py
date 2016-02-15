#! /usr/bin/env python

#template_making.py version 1.0 26 May 2014 Abby Moore

#This script is supposed to read the list of sequences from the groups of multiple sequences
#after you have gone through the list and marked the ones you want to keep.  There are two ways
#to mark them:
#	1) You can add a column (tab-separated), and mark the ones you want with * and the ones
#		you don't want with x in this first column.
#	2) You can simply delete the ones you don't want and do nothing to the ones you want.
#The two styles of marking can also be combined.

#Two possible formats of sequence list file:
'''
* [0]: x for don't want, * for want
2213_0 [1]: alignment name (this is just for choosing the sequences; we don't need it here)
748_0 [2]: sequence name (we need this)
extra bit at start [3]: any random notes written about the sequence

or (if the sequences we didn't want have just been deleted)

2213_0 [0]: alignment name (this is just for choosing the sequences; we don't need it here)
748_0 [1]: sequence name (we need this)
extra bit at start [2]: any random notes written about the sequence
'''


import sys #We want to be able to get information from the command line.

Usage = '''
template_making.py version 1.0
This script reads lists of sequences from groups of multiple sequences created by 
blast_parse_template_4.py once they have been annotated according to which sequences to
keep and which to discard.
template_making.py [annotated list of sequences] [file with sequences in it]
[file of unique sequences to which the desired sequences should be added]
'''

print("template_making.py")

if len(sys.argv) < 4:
	print(Usage)
else:
	ListFileName = sys.argv[1]
	SeqFileName = sys.argv[2]
	FinalSeqsName = sys.argv[3]

#The dictionaries and lists we will fill out:
SeqsWanted = [ ]
SeqDict = { }
SeqsOut = [ ]

#First, read in the file of the annotated list of sequences.
InFile = open(ListFileName, 'rU')
for Line in InFile:
	Line = Line.strip('\r').strip('\n').split('\t')
	if Line[0] == 'x': #This means you don't want the sequence.
		"do nothing"
	elif Line[0] == '*': #This means you do want the sequence.
		SeqsWanted.append(Line[2])
	elif Line[0] == "": #Skip blank lines.
		"do nothing"
	else: #We assume that if none of the above are true, only the lines with sequences we want
		#are still in the file.
		SeqsWanted.append(Line[1])
InFile.close()

print("The names of %d sequences that we want to add to the data matrix were read from the file %s.\n" \
	% (len(SeqsWanted), ListFileName))
sys.stderr.write("The names of %d sequences that we want to add to the data matrix were read from the file %s.\n" \
	% (len(SeqsWanted), ListFileName))

#Then read the file of sequences.
InFile = open(SeqFileName, 'rU')
for Line in InFile:
	Line = Line.strip('\n').strip('\r')
	#Each sequence is two lines, first the name, with a > symbol first,
	if Line[0] == ">":
		SeqName = Line[1:]
	#then the sequence itself.
	else:
		Seq = Line
		SeqDict[SeqName] = Seq
InFile.close()

print ("%d sequences were read from the file %s.\n" % (len(SeqDict), SeqFileName))
sys.stderr.write("%d sequences were read from the file %s.\n" % (len(SeqDict), SeqFileName))

#Now to make the list of formatted sequences to be added to the template sequence file:
for SeqName in SeqsWanted:
	Line = ">"+SeqName+"\n"+SeqDict[SeqName]+"\n"
	SeqsOut.append(Line)

#Now to open the outfile (_without_ writing over it) and add the new sequences.
OutFile = open(FinalSeqsName, 'a')
for Line in SeqsOut:
	OutFile.write(Line)
OutFile.close()

print("%d sequences were added to the file %s.\n" % (len(SeqsOut), FinalSeqsName))
sys.stderr.write("%d sequences were added to the file %s.\n" % (len(SeqsOut), FinalSeqsName))
