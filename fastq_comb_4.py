#! /usr/bin/env python

# This program takes forward and reverse sequences, reverse-complements the reverse
# sequences, reverses the Phred scores, and combines everything into a single fastq
# file.
#The start of the reverse sequence and the end of the forward sequence are compared
#and the overlapping bases are trimmed from the reverse sequence.
#At this point, there is nothing to mark the break between the forward and reverse
#sequences.  The forward sequence is just the first 88 bases (after the barcode) and the
#reverse sequence is the rest of the sequence.
#Version 3 is modified to also look for restriction sites later in the sequence, which would
#be from an insert that is shorter than 100bp (and thus part of the adapter was sequenced as well).
#Version4 is modified to check the sequences that do not appear to overlap using muscle.
#***Muscle will need to be installed for this to work and the line that calls muscle will need
#***to be modified to call the version of muscle that is installed on your computer.*******
#Version5 is modified to use Clustal Omega instead of muscle, in the hopes that it will be faster.
#This is called using clustalo.  This may need to be changed on other computers.

#The string from the first line of each sequence is:
# @HWI-ST558:60:C00B3ACXX:7:1206:8335:28258 1:N:0:ATCACG
# It has the following parts (with meaning according to Wikipedia):
# @HWI-ST558: the unique instrument name (same for all runs): 0
# 60: the run id (same within each run): 1
# C00B3ACXX: the flowcell id (same within each run): 2
# 7: flowcell lane (same within each sample): 3
# 1206: tile number within the flowcell lane (changes within a sample) **need this**: 4
# 8335: 'x'-coordinate of the cluster within the tile (changes within a sample) **need this**: 5
# 28258: 'y'-coordinate of the cluster within the tile (changes within a sample) **need this**: 6:0
# 1: the member of a pair, 1 or 2 (paired-end or mate-pair reads only): 6:1
# N: Y if the read fails filter (read is bad), N otherwise (I assume we do not have any Y samples.): 7
# 0: 0 when none of the control bits are on, otherwise it is an even number: 8
# ATCACG: index sequence: 9

import re #We will need the regular expressions module
import sys #We want to be able to output error messages to the screen.
import gzip #We want to be able to open zipped files.
import os #We need to be able to talk to the operating system
#We need these next two to be able to reverse-complement sequences.
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

#First, get the files names from the command line or print some instructions if not enough names are found..

Usage = '''
fastq_comb.py - version 4.0 3 April 2014
Cleans and concatenates sequence data to prepare it for fastq_parsebarcodes3.py.
Usage:
	fastq_comb_4.py (name(s) of file(s) with forward reads (containing barcodes)) > (name of log file you want to 
The file names can be either the full path, e.g., /home/moorea/... 
or the file name only if you are in that directory.

The sequence files can be zipped (with a ".gz" ending).
'''

if len(sys.argv) < 2:
	print Usage
else:
	InFileList = sys.argv[1:]

ReInFile = r"(\w+_L\d{3})_R\d(\w+)(\..+)"
SubInFileR = r"\1_R2\2\3"
SubOutFile = r"comb_\1\2\3"
SubTempFile = r"temp\1\2.fa"

print("fastq_comb_4.py\n")



# This first loop reads the forward sequences, the ones with the barcodes.
# It trims the sequences so they have 88 bp after the restriction site
# (so they will all be the same length once the barcode is removed)
# It rejects sequences without a restriction site and rejects the sequences with too many # quality symbols.

for InFileName1 in InFileList:
	InFileName2 = re.sub(ReInFile, SubInFileR, InFileName1)
	TempFileName = re.sub(ReInFile, SubTempFile, InFileName1)
	if InFileName1.endswith('.gz'):
		InFile1 = gzip.open(InFileName1)
	else:
		InFile1 = open(InFileName1,'rU')
	MyRe = r"(\w{4,8})(GATCC\w{83})\w+" #regular expression for saving 83 bases after the restriction site
	#***This is for BamHI.  If a different restriction enzyme is used, the GATCC needs to be changed accordingly.***
	MySub = r"\1\2" #regular expression for trimmed sequence
	MyRe2 = r"(\w+GGATC)(\w+)" #regular expression to see if there is a second restriction site
	#in sequence (i.e., insert is very short and sequence includes the adapters)
	MySub2 = r"\1"
	MySub3 = r"\2"
	Template = "AGATCGGAAGAGC"
	LineNum = 0
	ErrorNum = 0
	FailedSeqs= 0
	SeqNum = 0
	EndsTrimmed = 0
	#This is the data matrix we will fill.  It will be in fastq format with forward and
	#reverse ends combined, after reverse-complementation of reverse ends.
	CombinedSeqs = [ ]
	for Line in InFile1:
		Line = Line.strip('\n').strip('\r')
		SeqLine = (LineNum + 4) % 4
		if SeqLine != 1 and SeqLine != 3: #Code and "+" are copied directly
			CombinedSeqs.append(Line)
		elif SeqLine == 1: #Sequence is trimmed
			NewLine = re.sub(MyRe, MySub, Line)
			SeqNum += 1
			if len(NewLine) <= 96:
				#Checking to see if there are adapter sequences:
				NewLine2 = re.sub(MyRe2, MySub2, NewLine)
				if (NewLine == NewLine2) == False:
					#If there are, removing them:
					NewLine3 = re.sub(MyRe2, MySub3, NewLine)
					if len(NewLine3) > 13:
						if (NewLine3[0:13] == Template) == True:
							NewLine = NewLine2
							EndsTrimmed += 1
					else:
						if (NewLine3 == Template[0:len(NewLine3)]) == True:
							NewLine = NewLine2
							EndsTrimmed += 1
				CombinedSeqs.append(NewLine)
			else:
				CombinedSeqs.append("error")
				ErrorNum += 1
		else: #Quality line
			#First quality line is trimmed
			Line = Line[0:len(NewLine)]
			# Then it is examined for the presence of # symbols (very low quality sequences)
			CountFailed = Line.count('#')
			# If there are # symbols, then the sequence is rejected
			if CountFailed > 1:
				CombinedSeqs[LineNum-2] = 'error'
				FailedSeqs += 1
			CombinedSeqs.append(Line)
		LineNum += 1
		if LineNum%100000 == 0:
			sys.stderr.write("Processed the first %d lines.\n" % (LineNum-1))
	InFile1.close()
	sys.stderr.write("Forward file %s processed.\n" % InFileName1)
	sys.stderr.write("%d sequences did not have restriction sites.\n" % ErrorNum)
	sys.stderr.write("%d forward sequences were rejected due to low quality.\n" % FailedSeqs)
	sys.stderr.write("Adapter sequences were trimmed from the ends of %d sequences.\n" % EndsTrimmed)
	print("Forward file %s processed.\n" % InFileName1)
	print("%d sequences did not have restriction sites.\n" % ErrorNum)
	print("%d forward sequences were rejected due to low quality.\n" % FailedSeqs)
	print("Adapter sequences were trimmed from the ends of %d sequences.\n" % EndsTrimmed)

# This second loop reads the reverse sequences and appends the sequences and the
# quality information to the file with the forward sequences after reverse-complementing
# the reverse sequences and reversing the quality information.
# The last five bases are trimmed and sequences that still have too many # quality symbols (very low quality bases) are rejected
	LineNum = 0
	FPosition = [ ]
	RPosition = [ ]
	MyRe = r"@HWI-ST558:\d+:\w+:\d:(\d+):(\d+):(\d+) \d:\w:(\d+):\w+"
	MySub = r"\1, \2, \3"
	MyRe2 = r"(\w+GATCC)(\w+)"
	MySub2 = r"\2"
	MySub3 = r"\1"
	MyRe3 = r"(\w+GATCT)\w{4,8}GATCC\w+"
	Template2 = "CGACGCTCTTCCGATCT"
	FailedSeqs = 0
	TrimSeqs = 0
	UnTrimSeqs = 0
	EndsTrimmed2 = 0
	ToCheck = "no"
	if InFileName2.endswith('.gz'):
		InFile2 = gzip.open(InFileName2)
	else:
		InFile2 = open(InFileName2,'rU')
	for Line in InFile2:
		Line = Line.strip('\n').strip('\r')
		SeqLine = (LineNum + 4) % 4
		if SeqLine == 0: # These should be the position lines
			#We just need to make sure the coordinates are the same.  This can be taken out at the end
			# Although it is a good check to make sure that there is not some kind of a problem causing the lines to line up improperly.
			FPosition = re.sub(MyRe, MySub, CombinedSeqs[LineNum])
			RPosition = re.sub(MyRe, MySub, Line)
			if FPosition != RPosition:
				sys.stderr.write("Alas! There is a problem at line %d" % (LineNum))
				print("Alas! There is a problem at line %d" % (LineNum))
		elif SeqLine == 1: #Trim and append the sequence if there is a forward sequence
			if CombinedSeqs[LineNum] != 'error':
				Line = Line[0:95]
				ToDiscard = 0
				TrimFor = 0
				#taking the reverse complement
				Seq2F = Seq(Line, IUPAC.unambiguous_dna)
				Seq2RC = str(Seq2F.reverse_complement())
				#Checking to see if there are adapter sequences:
				Seq2RC2 = re.sub(MyRe2, MySub2, Seq2RC)
				if (Seq2RC == Seq2RC2) == False:
					#If so, removing them:
					Seq2RC3 = re.sub(MyRe2, MySub3, Seq2RC)
					if len(Seq2RC3) <= 18:
						Seq2RC = Seq2RC2
						EndsTrimmed2 += 1
					else:
						Seq2RC4 = re.sub(MyRe3, MySub3, Seq2RC)
						if (len(Seq2RC4) < len(Seq2RC3)):
							Seq2RC = Seq2RC2
							EndsTrimmed2 += 1
				#This next section looks to see if the sequences overlap.
				SubSeq = Seq2RC[0:6]
				Seq1F = CombinedSeqs[LineNum]
				if (Seq1F.count(SubSeq) >= 1):
					StrNum = 7
					while((Seq1F.count(Seq2RC[0:StrNum]) >= 1) and (StrNum <= len(Seq2RC))):
						StrNum += 1
					ToTrim = StrNum - 1
					#If they do overlap, it removes one copy of the
					#overlapping region.
					if (Seq1F[-ToTrim:] == Seq2RC[0:ToTrim]):
						CombinedSeqs[LineNum] += Seq2RC[ToTrim:]
						TrimSeqs += 1
						ToCheck = "no"
						ToDiscard = 95 - len(Seq2RC[ToTrim:])
					#If the sequences do not appear to overlap, it needs to check them.
					else:
						ToCheck = "yes"
				else:
					ToCheck = "yes"
				if ToCheck == "yes":
					#If the sequences seem like they don't overlap, it is probably because of a mismatch.
					#So this needs to be checked by actually aligning them.
					#First, writing the file to be aligned and calling muscle:
					SeqFile = ">For\n"+Seq1F+"\n>Rev\n"+Seq2RC
					TempFile = open(TempFileName, 'w')
					TempFile.write(SeqFile)
					TempFile.close()
					SeqDictTemp = { }
					#This line is for running on my laptop
					#OutLine = "muscle3.8.31_i86linux64 -in "+TempFileName+" -out "+TempFileName+".out -quiet"
					#This line is for running on oscar
					OutLine = "muscle -in "+TempFileName+" -out "+TempFileName+".out -quiet"
					os.popen(OutLine,'r')
					#Second, reading the alignment
					InFileTemp = open((TempFileName+".out"), 'rU')
					SeqName = ""
					SeqTemp = ""
					for Line in InFileTemp:
						Line = Line.strip('\r').strip('\n')
						if Line != "":
							if Line[0] ==">":
								if SeqTemp != "":
									SeqFor = SeqTemp
									SeqTemp = ""
							else:
								SeqTemp += Line
					SeqRev = SeqTemp
					InFileTemp.close()
					#Third, examining the sequences:
					#We need to make a consensus sequence and count the number of mismatches.
					CombSeq = ""
					ToTrim = 0
					NumMismatch = 0
					for bp in (range(len(SeqFor))):
						if SeqFor[bp] == "-":
							CombSeq += SeqRev[bp]
						elif SeqRev[bp] == "-":
							CombSeq += SeqFor[bp]
						elif SeqFor[bp] == SeqRev[bp]:
							CombSeq += SeqFor[bp]
						elif (SeqFor[bp] == "A" and SeqRev[bp] == "C") or (SeqRev[bp] == "A" and SeqFor[bp] == "C"):
							CombSeq += "M"
							NumMismatch += 1
						elif (SeqFor[bp] == "A" and SeqRev[bp] == "G") or (SeqRev[bp] == "A" and SeqFor[bp] == "G"):
							CombSeq += "R"
							NumMismatch += 1
						elif (SeqFor[bp] == "A" and SeqRev[bp] == "T") or (SeqRev[bp] == "A" and SeqFor[bp] == "T"):
							CombSeq += "W"
							NumMismatch += 1
						elif (SeqFor[bp] == "C" and SeqRev[bp] == "G") or (SeqRev[bp] == "C" and SeqFor[bp] == "G"):
							CombSeq += "S"
							NumMismatch += 1
						elif (SeqFor[bp] == "C" and SeqRev[bp] == "T") or (SeqRev[bp] == "C" and SeqFor[bp] == "T"):
							CombSeq += "Y"
							NumMismatch += 1
						elif (SeqFor[bp] == "G" and SeqRev[bp] == "T") or (SeqRev[bp] == "G" and SeqFor[bp] == "T"):
							CombSeq += "K"
							NumMismatch += 1
						elif SeqFor[bp] == "N":
							CombSeq += SeqRev[bp]
						elif SeqRev[bp] == "N":
							CombSeq += SeqFor[bp]
						else:
							print("Error at position %d, forward sequence has %s and reverse sequence has %s, but this cannot be resolved.\n" \
								% (bp, SeqFor[bp], SeqRev[bp]))
					#If there were not too many mismatches, we assume that the sequences really do overlap,
					#and combine them.
					if NumMismatch < 5:
						#But first, we need to make sure that the forward and reverse sequences
						#were properly trimmed.
						#Making sure that the reverse sequence does not extent past the start of the 
						#forward sequence
						if SeqFor[0] == "-":
							ForRe = r"\-+(\w+[A-Z\-]+)"
							ForSub = r"\1"
							SeqFor = re.sub(ForRe, ForSub, SeqFor)
							CombSeq = CombSeq[-len(SeqFor):]
						#And that the forward sequence does not extend past the end of the reverse sequence
						if SeqRev[-1] == "-":
							RevRe = r"([A-Z\-]+?)\-+$"
							RevSub = r"\1"
							SeqRev = re.sub(RevRe, RevSub, SeqRev)
							CombSeq = CombSeq[0:len(SeqRev)]
							CombinedSeqs[LineNum+2] = CombinedSeqs[LineNum+2][0:len(SeqRev)]
						#Once that has been done, the combined sequence can be added to the file.
						CombinedSeqs[LineNum] = CombSeq
						TrimSeqs += 1
						#And the proper number of bases of the quality score to discard can be determined.
						ToDiscard = 95 - (len(CombSeq) - len(CombinedSeqs[LineNum+2]))
					#But, if we have a lot of mismatches, we assume that that is because the sequences
					#actually don't overlap and they are joined end to end.
					else:
						CombinedSeqs[LineNum] += "NN"
						CombinedSeqs[LineNum] += Seq2RC
						CombinedSeqs[LineNum+2] += "##"
						#It is still possible that a few bases were trimmed off the start of the reverse
						#sequence (and either the trimming was a mistake or the combining the sequences
						#end to end is a mistake), but we need to take that into account so the quality
						#score and the sequence will be the same length.
						ToDiscard = 95-len(Seq2RC)
						UnTrimSeqs += 1
			else:
				ToTrim = 0 #This is to prevent it from crashing on SeqLine3 if
				#CombinedSeqs[LineNum] is "error"
		elif SeqLine == 3: #Check the sequence and trim and append the quality score
			#If the sequence needs to be examined:
			if CombinedSeqs[LineNum-2] != 'error':
				Line = Line[0:95]#Takes the start of the line (corresponding to the
				#bases in the sequence)
				Line = Line[::-1]# Reverses the line (so the quality scores are in the
				#same order as the base in the reverse-complemented sequence).
				CountFailed = Line.count('#')
				# If there is more than one # symbol (very low quality base), then the sequence is rejected
				if CountFailed > 1:
					CombinedSeqs[LineNum-2] = 'error'
					FailedSeqs += 1
					#If the sequence quality is acceptable, then the appropriate number of positions of the
					#quality score are added to the quality score from the forward sequence.
				else:			
					CombinedSeqs[LineNum] += Line[ToDiscard:]
				if len(CombinedSeqs[LineNum]) != len(CombinedSeqs[LineNum-2]):
					if CombinedSeqs[LineNum-2] != "error":
						print("Error: Sequence: %d, Qual Score: %d" % (len(CombinedSeqs[LineNum-2]), len(CombinedSeqs[LineNum])))
						print(SeqFor)
						print(SeqRev)
						print(CombSeq)
						print(CombinedSeqs[LineNum-2])
						print(CombinedSeqs[LineNum])
						print(ToDiscard)
						print(TrimFor)
		LineNum += 1
		if LineNum%100000 == 0:
			sys.stderr.write("Processed the first %d lines.\n" % (LineNum-1))
	InFile2.close()
	sys.stderr.write("Reverse file %s processed.\n" % InFileName2)
	sys.stderr.write("%d reverse sequences were rejected due to low quality.\n" % FailedSeqs)
	sys.stderr.write("The ends were trimmed from %d sequences.\n" % EndsTrimmed2)
	print("Reverse file %s processed.\n" % InFileName2)
	print("%d reverse sequences were rejected due to low quality.\n" % FailedSeqs)
	print("The ends were trimmed from %d sequences.\n" % EndsTrimmed2)

	FailedSeqs = CombinedSeqs.count('error')
	sys.stderr.write("In total, %d sequences were processed.\n" % SeqNum)
	sys.stderr.write("%d of these were trimmed because they overlapped and %d of them were not.\n" % (TrimSeqs, UnTrimSeqs))
	sys.stderr.write("In total, %d sequences were rejected due to low quality.\n" % FailedSeqs)
	print("In total, %d sequences were processed.\n" % SeqNum)
	print("%d of these were trimmed because they overlapped and %d of them were not.\n" % (TrimSeqs, UnTrimSeqs))
	print("In total, %d sequences were rejected due to low quality.\n" % FailedSeqs)

#Now make a new list (CombinedSeqsChecked), which includes only the sequences that are not rejected.

	OldLine = 0
	CombinedSeqsChecked = [ ]
	NumList = range(SeqNum)
	SeqLines = range(4)
	FailedSeqs = 0
	GoodSeqs = 0
	for SeqN in NumList:
		if CombinedSeqs[SeqN*4+1] == "error": # Skip all four items of the sequence if it is rejected
			OldLine += 4
			FailedSeqs += 1
		else: #Copy all four lines of the sequence if it is not rejected.
			for Line in SeqLines:
				CombinedSeqsChecked.append(CombinedSeqs[OldLine])
				OldLine += 1
			GoodSeqs += 1

# Now write the list of good sequences to a file.
	OutFileName = re.sub(ReInFile, SubOutFile, InFileName1)
	OutFile = gzip.open(OutFileName, 'w')
	for Line in CombinedSeqsChecked:
		OutFile.write(Line)
		OutFile.write('\n')
	OutFile.close()

	if (SeqNum-FailedSeqs) != GoodSeqs:
		sys.stderr.write("Error: %d good sequences written to a file, but %d good sequences were found.\n" % \
		(GoodSeqs, SeqNum-FailedSeqs))
		print("Error: %d good sequences written to a file, but %d good sequences were found.\n" % (GoodSeqs, SeqNum-FailedSeqs))

	sys.stderr.write("Output file %s written.\n" % OutFileName)
	sys.stderr.write("%d good sequences were written to this file; %d sequences were rejected.\n\n\n" % (GoodSeqs, FailedSeqs))
	print("Output file %s written.\n" % OutFileName)
	print("%d good sequences were written to this file; %d sequences were rejected.\n\n\n" % (GoodSeqs, FailedSeqs))

