#! /usr/bin/env python

#This file is supposed to read the .clstats output from the first rtd_run.py run, which is only for the first 30 bases.
#It goes through that file to figure out which groups of sequences are present in most of the individuals (with a range that the
#user inputs) and are present fewer than a certain number of times.
#Then it looks through the .cluni file to find the sequences belonging to these groups and makes a dictionary
#linking the sequence to the group
#Then it looks through the file of uniqued sequences and finds the sequences whose first 30 bases correspond to one of
#the sequences in the selected groups and makes a new uniqued sequences file with the sequences labeled.
#Then it makes separate outfiles for each group number.
#Finally, it writes a shell script to analyze these various files and parse the output of
#the analyses to remove duplicate sequences and make a new file for use as a reference in
#the subsequent analyses.
#Version 5 is changed so that the final comparison of the uniqued sequences takes place with BLAST instead of
#with BLAT.
#Version 5.0, 18 March 2014

import sys #We want to be able to output error messages to the screen and get input from the command line.
import gzip #We want to be able to open zipped files.
import os #We want to be able to talk to the shell.

from collections import defaultdict #We need defaultdict to be able to create dictionaries with multiple levels.

Usage = '''
bam30_uniq_split_5.py is supposed to take the output from an rtd_run.py analysis of the first 30 bases
of the sequence data, look for the groups of sequences that are present in a large proportion of the individuals,
retreive those whole sequences according to group from the uniqued sequence files, run blat on those sequences,
and parse the output to remove duplicate sequences.
Input should be in the following format:
bam30_uniq_split_5.py [lower limit of group size] [upper limit of group size] [lower limit of group size for smaller clusters] \
[.clstats file name from rtd_run.py] [.cluni file name from rtd_run.py] [file containing uniqued (full) sequences] \
[outfile path] [outfile prefix--without the path]
'''

if len(sys.argv) <9:
	print Usage
else:
	IndLowLim = int(sys.argv[1])
	IndUpLim = int(sys.argv[2])
	IndLowLimSm = sys.argv[3]
	InFileCLSt = sys.argv[4]
	InFileCLU = sys.argv[5]
	InFileUSeqs = sys.argv[6]
	OutFilePath = sys.argv[7]
	OutFilePre = sys.argv[8]

#IndLowLim = 130
#IndUpLim = 154
#InFileCLSt = "/mnt/disk2/trial/rtd_blat-minScore20-minMatch2-repMatch1000000_l8_l8_mcl_I2.0_0.05dirt_100indiv_200seq.clstats"
#InFileCLU = "rtd_blat-minScore20-minMatch2-repMatch1000000_l8_l8_mcl_I2.0.cluni"
#InFileUSeqs = "/home/moorea/trial/seqsbyind1/228uniqued.uniqued.gz"
#OutFilePre = "/mnt/disk2/trial/2013_06_13"

#Reading the clstats file to find which groups of sequences are present in between IndLowLim and IndUpLim 
#individuals (inclusive), include fewer than 1000 unique sequences, and have a "dirt" value lower than 0.05
#Only those individuals with a dirt value below 0.05 will be in the bam file, but the minimum number of individuals
#varies according to the settings specified with rtd_run.py and the default is 100.

print("bam30_uniq_split_5.py")

GrpList=[ ]
MaxSeqNum = 1000
MaxDirt = 0.05
InFile = open(InFileCLSt, 'rU')
for Line in InFile:
	Line = Line.strip('\r').strip('\n').split('\t')
	if int(Line[1]) < MaxSeqNum:
		if int(Line[2]) >= IndLowLim and int(Line[2]) <= IndUpLim:
			if float(Line[3]) <= MaxDirt:
				GrpList.append(Line[0])
InFile.close()

print("There are %d groups of sequences that are present in between %d and %d individuals, include fewer than \
%d different sequences, and have a percentage \"dirt\" less than %.2f." \
% (len(GrpList), IndLowLim, IndUpLim, MaxSeqNum, MaxDirt))
print("These groups of sequences will be used for further analyses.")
sys.stderr.write("There are %d groups of sequences that are present in between %d and %d individuals, include fewer than \
%d different sequences, and have a percentage \"dirt\" less than %.2f." \
% (len(GrpList), IndLowLim, IndUpLim, MaxSeqNum, MaxDirt))
sys.stderr.write("These groups of sequences will be used for further analyses.\n")
#GrpList = ['569']


#Now to retrieve the sequences for each group and make a dictionary:

SeqDict = { }

InFile = open(InFileCLU,'rU')
for Line in InFile:
	Line = Line.strip('\r').strip('\n').split('\t')
	GrpNum = Line[0]
	if (GrpNum in GrpList) == True:
		SeqDict[Line[2]] = GrpNum

for GrpNum in GrpList:
	#Using some shell commands to retrieve a list of unique sequences in that group.
	OutLine = "cut -f 1,3 "+InFileCLU+" | grep "+GrpNum+" | cut -f 2 | uniq"
	GrpSeqStr = os.popen(OutLine,'r').read()
	GrpSeqList = GrpSeqStr.split('\n')
	#append these sequences to the dictionary
	for Seq in GrpSeqList:
		if len(Seq) > 2:
			SeqDict[Seq] = GrpNum
			
print("There are a total of %d 30 bp sequences that belong to these groups.\n" % len(SeqDict))
sys.stderr.write("There are a total of %d 30 bp sequences that belong to these groups.\n" % len(SeqDict))

#Now to make a new file of uniqued sequences, containing only those which occur in these groups and labeled
#according to group.
FullSeqDict = { }
TotalSeqs = 0 
if InFileUSeqs.endswith('.gz'):
	InFile = gzip.open(InFileUSeqs)
else:
	InFile = open(InFileUSeqs,'rU')
for Line in InFile:
	Line = Line.strip('\n').strip('\r')
	Seq = Line[0:30]
	GrpNum = SeqDict.get(Seq, -99) #Seeing if the sequence belongs to a group.
	if GrpNum != -99:
		TotalSeqs += 1
		try:
			FullSeqDict[GrpNum] += Line+'\n' #Adding the sequence to the entry for that group,
			#if it exists.
		except KeyError:
			FullSeqDict[GrpNum] = Line+"\n" #If not, creating an entry for that group.
InFile.close()

print ("A total of %d full sequences were found that belong to these groups.\n" % TotalSeqs)
sys.stderr.write("A total of %d full sequences were found that belong to these groups.\n" % TotalSeqs)

#Now to write the files:
for GrpNum in GrpList: #For each group number
	OutFileName = OutFilePath+OutFilePre+"_"+GrpNum+"uniqued.gz"
	OutFile = gzip.open(OutFileName, 'w')
	OutFile.write(FullSeqDict[GrpNum])
	OutFile.close()
	print("Outfile %s written.\n" % OutFileName)
	sys.stderr.write("Outfile %s written.\n" % OutFileName)

#Now to write the shell script:
GrpScript = [ ]
#This says that it is a shell scrip
Line = "#! /bin/bash\n"
#This makes a subject file for the rtd pipeline so it will run
Line3 = "new_subj_rtd.py\n"
GrpScript.append(Line)
#Each group of sequences is analyzed separately:
for GrpNum in GrpList:
	OutFileName = OutFilePath+OutFilePre+"_"+GrpNum+"uniqued.gz"
	OutFileDir = OutFilePath+GrpNum
	#This runs the rtd pipeline for the first time, but it will crash:
	Line1 = "rtd_run.py -np 5 "+OutFileDir+" "+OutFileName+"\n"
	#Then we change into that directory so we can run new_subj_rtd.py
	Line2 = "cd "+OutFileDir+"\n"
	#After running the rtd pipeline for the 2nd time, we retrieve the template sequences
	#for the clusters we want to analyze further
	Line4 = "clusters2template.py "+IndLowLimSm+" "+GrpNum+" "+OutFileDir\
	+"/*.clstats "+OutFileDir+"/*.cluni "+OutFilePath+OutFilePre+"\n" 
	GrpScript.append(Line1)
	GrpScript.append(Line2)
	GrpScript.append(Line3)
	GrpScript.append(Line1)
	GrpScript.append(Line4)
#After analyzing all of the groups, we put all of the template sequences into one
#file for further filtering.
Line = "template_concat.py "+OutFilePath+OutFilePre+" "+OutFilePath+OutFilePre+"template_*.fa\n"
GrpScript.append(Line)
#We make that file into a blast database, so we can blast the sequences against themselves
#to see if any are the same or very similar.  (We need to do this, because the groups were
#originally created by looking at only the first 30bp of the sequences.)
Line = "makeblastdb -in "+OutFilePath+OutFilePre+".fa -out "+OutFilePath+OutFilePre+"db -dbtype nucl\n"
GrpScript.append(Line)
#Then we blast the sequences against themselves.
Line = "blastn -db "+OutFilePath+OutFilePre+"db -query "+OutFilePath+OutFilePre+".fa -out "+OutFilePath+OutFilePre+"out.txt -outfmt \'6 std qlen slen\'\n"
GrpScript.append(Line)
#And blast them against the primer sequences to remove any remaining primers.
Line = "blastn -db "+OutFilePath+OutFilePre+"db -query "+OutFilePath+"primers.fa -out "+OutFilePath+"primersout.txt -outfmt \'6 std qlen slen\'\n"
GrpScript.append(Line)
#Then we need to parse the blast output files and come up with a file of sequences that don't have matches and a set of alignments
#of the sequences that do have matches for further examination.
#But first, we need to make the folder in which to put the alignments
Line = "mkdir "+OutFilePath+"/altemp\n"
GrpScript.append(Line)
Line = "blast_parse_template_4.py "+OutFilePath+OutFilePre+"_cleaned.fa "+OutFilePath+"altemp/ "+OutFilePath+OutFilePre+"out.txt "+OutFilePath+"primersout.txt "+OutFilePath+OutFilePre+".fa\n"
GrpScript.append(Line)

OutFileName = OutFilePath+OutFilePre+"script.sh"
OutFile=open(OutFileName,'w')
for Line in GrpScript:
	OutFile.write(Line)
OutFile.close()
print("The next step is to run the following script: %s \n" % (OutFileName))
print("In order to do this, you must type the following:\n\tchmod u+x %s\n\t%s\n" %\
	(OutFileName,OutFileName))
print("This script will not delete anything, so you do not need to worry.\n")
sys.stderr.write("The next step is to run the following script: %s \n" % (OutFileName))
sys.stderr.write("In order to do this, you must type the following:\n\tchmod u+x %s\n\t%s\n" %\
	(OutFileName,OutFileName))
sys.stderr.write("This script will not delete anything, so you do not need to worry.\n")
