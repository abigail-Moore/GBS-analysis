#! /usr/bin/env python

#This script takes the set of template files (one for each cluster from the rtd_run.py run)
#and combines them into a single fasta file.  It also writes (but does not execute) a shell script
#for getting rid of the original template files.

import sys #We want to be able to write messages on the screen.

Usage = '''
template_concat.py version 1.0 24 Dec. 2013
This script reads a list of template fasta files and concatenates them.  It also writes a script
to get red of the original template files.
template_concat.py [output file prefix] [list of template files]
'''

if len(sys.argv) <3:
	print Usage
else:
	OutFilePre = sys.argv[1]
	TemplateList = sys.argv[2:]

print("template_concat.py")

OutFileNameScript = OutFilePre+"_del.sh"
OutFileNameConcat = OutFilePre+".fa"

#Making a list of the sequences in the various files:

SeqList = [ ]

for InFileName in TemplateList:
	InFile = open(InFileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n')
		SeqList.append(Line)

print("A total of %d template sequences were read from %d input files.\n" \
	% (len(SeqList), len(TemplateList))) 
sys.stderr.write("A total of %d template sequences were read from %d input files.\n" \
	% (len(SeqList), len(TemplateList)))

ScriptList = [ ]

for InFileName in TemplateList:
	Line = "rm "+InFileName
	ScriptList.append(Line)

OutFile = open(OutFileNameConcat, 'w')
for Seq in SeqList:
	OutFile.write(Seq+"\n")
OutFile.close()

print("The sequences were written to the file %s.\n" % (OutFileNameConcat))
sys.stderr.write("The sequences were written to the file %s.\n" % (OutFileNameConcat))

OutFile = open(OutFileNameScript, 'w')
for Line in ScriptList:
	OutFile.write(Line+"\n")
OutFile.close()

print("The script for removing the intermediate files was written to %s.\n" % (OutFileNameScript))
print("In order to remove these files, type the following:\n\tchmod u+x %s\n\t%s\n" % \
	(OutFileNameScript, OutFileNameScript))
print("However, it is good to open the file first, to make sure you will not delete other files that you want to keep!!\n")
sys.stderr.write("The script for removing the intermediate files was written to %s.\n" % (OutFileNameScript))
sys.stderr.write("In order to remove these files, type the following:\n\tchmod u+x %s\n\t%s\n" % \
	(OutFileNameScript, OutFileNameScript))
sys.stderr.write("However, it is good to open the file first, to make sure you will not delete other files that you want to keep!!\n")
