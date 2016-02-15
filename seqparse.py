#! /usr/bin/env python

#seqparse.py version 1.0 1 Aug. 2014 Abby Moore
#This script is supposed to go through alignments of sequences for a single individual
#and determine which sequences are present in that individual.
#At this point, it makes its own decision how many sequences the individual has, although
#this may change.
#This is apparently the current working version.

from collections import defaultdict #We want to be able to make dictionaries with multiple levels.
import sys #We want to be able to get input from the command line.
import re #We want to be able to use regular expressions. (maybe)

'''Sequence names should be of the form 353632_CA53_1_13, in which:
353632 [0]: sequence number--don't need
CA53 [1]: population name (taxon followed by collection number)
1 [2]: individual number
13 [3]: number of times that sequence is present'''

Usage = '''
seqparse.py version 1.0 parses fasta files made by blast_parse_output.py to reconstruct
the actual sequence(s) present in each individual
seqparse.py [minimum number of sequences that must have an alternative nucleotide for
the position to be considered polymorphic] [minimum proportion (<1) of sequences that must have
a given allele for that allele to be considered to be present] [ploidy] [individual name] 
[folder with the final sequence files] [list of sequence files to parse] '''

if len(sys.argv) < 5:
	print Usage
else:
	NumMinO = int(sys.argv[1])
	PercMin = float(sys.argv[2])
	Ploidy = int(sys.argv[3])
	IndName = sys.argv[4]
	OutFolder = sys.argv[5]
	SeqFileList = sys.argv[6:]

SeqFileRe = r".*/*a(\w+)(\.fa)"
SeqFileSub = r"a\1\2"
LocusNameSub = r"\1"
PloidyErrorList = [ ] #List for ploidy errors (to print later).
GenErrorList = [ ] #List of other errors (to print later).
SeqInfoList = [ ] #List of information about sequences (not errors).

Test = "no" #When this is yes, it prints out a bunch more stuff.

if OutFolder[-1:] != "/":
	OutFolder += "/"

for SeqInFileName in SeqFileList:
	SeqOutFileName = OutFolder+re.sub(SeqFileRe, SeqFileSub, SeqInFileName)
	LocusName = re.sub(SeqFileRe, LocusNameSub, SeqInFileName)
	if Test == "yes":
		print LocusName
	#Get the dictionaries and things ready that we need to use:
	SeqList = [ ] #The list of sequences in the alignment
	SeqDict = { } #The dictionary of sequences in which the sequence names are keys and the
	#sequences are values
	SeqNumDict = { } #The dictionary of sequences in which the sequence names are keys and the
	#number of times a sequence is present are values
	TotalSeqs = 0
	CmpdSeq = [ ] #The final sequence in list form
	AmbigDict = defaultdict(dict) #Dictionary of variable sites in which the higher level keys
	#are the positions in the alignment (starting with zero), the lower level keys are the nucleotides,
	#and the values are lists of sequences that have that nucleotide.
	NewSeqDict = { } #Dictionary of the new sequence(s) from this locus.
	NNDict = { } #Dictionary of positions where sequences are joined.
	JoinedSeqRange = [ ] #The area of overlap between two joined sequences.
	NumMin = NumMinO #Resetting NumMin to the original value.
	#First, read the fasta file for that locus.
	InFile = open(SeqInFileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\n').strip('\r')
		if Line[0] == '>':
			SeqName = Line[1:]
			SeqList.append(SeqName)
			SeqDict[SeqName] = ""
			SeqNameSplit = SeqName.split("_")
			NumSeqs = int(SeqNameSplit[3])
			SeqNumDict[SeqName] = NumSeqs
			TotalSeqs += NumSeqs
		else:
			SeqDict[SeqName] += Line
	InFile.close()
	#Second, look for differences among the sequences:
	SeqRange = range(len(SeqDict[SeqName]))
	#If we have too few sequences, print an error message.
	if TotalSeqs < (NumMin+2):
		print("Error!  We need %d sequences to make a determination of the nucleotide, but only %d sequences are present!!\n" % (NumMin+2, TotalSeqs))
		#sys.stderr.write("Error!  We need %d sequences to make a determination of the nucleotide, but only \
		#	%d sequences are present!!\n" % (NumMin+2, TotalSeqs))
		GenError = "Locus "+LocusName+" only has "+str(NumMin)+" sequences, so no final sequences were recovered.\n"
		GenErrorList.append(GenError)
	#If not, analyze the sequences.
	else:
		#First, increasing the minimum number of times a SNP has to be present to accept it, 
		#if the number of sequences at that locus is very high.
		NumMinNew = int(TotalSeqs*PercMin)
		if TotalSeqs > 10000:
			NumMin = 1000
		elif NumMinNew > 100:
			NumMin = 100
		elif NumMinNew > NumMin:
			NumMin = NumMinNew
		if Test == "yes":
			print NumMin
		#then, check to see if the sequences are concatenated from two sequences that do not overlap:
		for SeqName in SeqDict.keys():
			try:
				NNPos = SeqDict[SeqName].index('NN')
				try:
					NNDict[NNPos] += SeqNumDict[SeqName]
				except KeyError:
					NNDict[NNPos] = SeqNumDict[SeqName]
			except ValueError:
				"do nothing"
		#Now to look through the list of positions of 'NN', if they exist
		if NNDict != { }:
			NNList = [ ]
			for NNPos in NNDict:
				if NNDict[NNPos] >= NumMin:
					NNList.append(NNPos)
			NNList.sort()
			if Test == "yes":
				print NNList
			#We can only deal with them this way if there are two positions
			if len(NNList) == 2:
				NNSeqStart = { } #Dictionary of the starts of the reverse sequences
				NNSeqEnd = { } #Dictionary of the ends of the forward sequences
				for SeqName in SeqDict.keys():
					#Look through the sequences to find the overlapping bits
					NNPartTemp = SeqDict[SeqName][NNList[0]:(NNList[1]+2)]
					NNPart = ""
					for Base in NNPartTemp:
						if Base != "-":
							NNPart += Base
					#And classify them as the start of the reverse sequence
					if NNPart[:2] == "NN":
						try:
							NNSeqStart[NNPart] += SeqNumDict[SeqName]
						except KeyError:
							NNSeqStart[NNPart] = SeqNumDict[SeqName]
					#or the end of the forward sequence
					elif NNPart[-2:] == "NN":
						try:
							NNSeqEnd[NNPart] += SeqNumDict[SeqName]
						except KeyError:
							NNSeqEnd[NNPart] = SeqNumDict[SeqName]
				#If there is only one start of the reverse sequence, take it
				if Test == "yes":
					print NNSeqStart
					print NNSeqEnd
				if len(NNSeqStart) == 1:
					SeqStart = NNSeqStart.keys()[0]
				#If not, see if the alternatives are sequencing errors
				else:
					SeqStartList = [ ]
					for NNPart in NNSeqStart:
						if NNSeqStart[NNPart] >= NumMin:
							SeqStartList.append(NNPart)
					if len(SeqStartList) == 1:
						SeqStart = SeqStartList[0]
					#If not, it becomes more complicated.
					else:
						SeqStart = 'NN'
						print("There are too many potential beginnings for the reverse sequence, so that part of the sequence was trimmed.\n")
						#sys.stderr.write("There are too many potential beginnings for the reverse sequence, so that part of the sequence was trimmed.\n")
						GenError = "Locus "+LocusName+" has too many potential beginnings for the reverse sequence, so that part of the sequence was trimmed.\n"
						GenErrorList.append(GenError)
				if Test == "yes":
					print SeqStart
				#And the same for the ends of the forward sequences
				if len(NNSeqEnd) == 1:
					SeqEnd = NNSeqEnd.keys()[0]
				else:
					SeqEndList = [ ]
					for NNPart in NNSeqEnd:
						if NNSeqEnd[NNPart] >= NumMin:
							SeqEndList.append(NNPart)
					if len(SeqEndList) == 1:
						SeqEnd = SeqEndList[0]
					else:
						SeqEnd = 'NN'
						print("There are too many potential endings for the forward sequence, so that part of the sequence was trimmed.\n")
						#sys.stderr.write("There are too many potential endings for the forward sequence, so that part of the sequence was trimmed.\n")
						GenError = "Locus "+LocusName+" has too many potential endings for the forward sequence, so that part of the sequence was trimmed.\n"
						GenErrorList.append(GenError)
						#Although this is a problem I should be able to deal with***
				if Test == "yes":
					print SeqEnd
				#Figure out what the new middle of the sequence is
				JoinedSeq = SeqEnd + SeqStart[2:]
				#and which positions of the alignment are affected.
				JoinedSeqRange = range(NNList[0],(NNList[1]+2))
				if Test == "yes":
					print JoinedSeq
					print JoinedSeqRange
			if len(NNList) == 1:
				print ("The sequences were joined, but only one or the other is present, so they can be treated like normal sequences.\n")
				#sys.stderr.write("The sequences were joined, but only one or the other is present, so they can be treated like normal sequences.\n")
				GenError = "Locus "+LocusName+" was joined, but all sequences are in one direction, so it should be checked to make sure it is correct.\n"
				GenErrorList.append(GenError)		
			elif len(NNList) > 2:
				JoinedSeq = "NN"
				JoinedSeqRange = range(NNList[0],(NNList[1]+2))
				print ("There are too more than two places where the sequences have been artificially joined, so %d bases in the middle of the sequence are ambiguous.\n" % (len(JoinedSeqRange)))
				#sys.stderr.write("There are too more than two places where the sequences have been artificially joined, so %d bases in the middle of the sequence are ambiguous.\n" % (len(JoinedSeqRange)))
				GenError = "Locus "+LocusName+" has more than two places where the sequences have been artificially joined, so "+str(len(JoinedSeqRange))+" bases in the middle of the sequence are ambiguous.\n"
				GenErrorList.append(GenError)
				if Test == "yes":
					print JoinedSeq
					print JoinedSeqRange
		#Then look through the sequences, skipping the non-overlapping bits, if they exist
		for SeqPos in SeqRange:
			if (SeqPos in JoinedSeqRange) == False:
				#resetting the counters (This could be done more elegantly using some kind of list function thingie***)
				A = 0
				C = 0
				G = 0
				T = 0
				Gap = 0
				NucTemp = [ ]
				for SeqName in SeqList:
					#determining which nucleotide is present at that position
					if SeqDict[SeqName][SeqPos] == 'A':
						A += SeqNumDict[SeqName]
					elif SeqDict[SeqName][SeqPos] == 'C':
						C += SeqNumDict[SeqName]
					elif SeqDict[SeqName][SeqPos] == 'G':
						G += SeqNumDict[SeqName]
					elif SeqDict[SeqName][SeqPos] == 'T':
						T += SeqNumDict[SeqName]
					elif SeqDict[SeqName][SeqPos] == '-':
						Gap += SeqNumDict[SeqName]
				#Deciding which nucleotides are present often enough not to ignore
				NonGap = A + C + G + T
				#If the position is mostly a gap, we don't want to include it in the final sequence.
				if 10*NonGap < Gap:
					NucTemp.append('Gap')
				else:
					if A >= NumMin:
						NucTemp.append('A')
					if C >= NumMin:
						NucTemp.append('C')
					if G >= NumMin:
						NucTemp.append('G')
					if T >= NumMin:
						NucTemp.append('T')
					if Gap >= NumMin:
						NucTemp.append('Gap')
				#adding the correct nucleotide to the sequence
				#This is trivial if there is only one nucleotide present
				if len(NucTemp) == 1:
					CmpdSeq.append(NucTemp[0])
				#And can't be done if there are no nucleotides present
				elif len(NucTemp) < 1:
					print("ERROR!  No nucleotides were found for position %d at locus %s!!\n" % (SeqPos, LocusName))
					sys.stderr.write("ERROR!  No nucleotides were found for position %d at locus %s!!\n" % (SeqPos, LocusName))
					CmpdSeq.append('N')
					GenError = "Locus "+LocusName+" has no nucleotides that were found at position "+str(SeqPos)+".\n"
					GenErrorList.append(GenError)
				#But needs to be a list if more than one nucleotide is present
				else:
					CmpdSeq.append(NucTemp)
					#And there needs to be an entry in AmbigDict for that position as well.
					for SeqName in SeqList:
						try:
							AmbigDict[SeqPos][SeqDict[SeqName][SeqPos]].append(SeqName)
						except KeyError:
							AmbigDict[SeqPos][SeqDict[SeqName][SeqPos]] = [ ]
							AmbigDict[SeqPos][SeqDict[SeqName][SeqPos]].append(SeqName)
			elif SeqPos == JoinedSeqRange[0]:
				#add the joined bit back
				CmpdSeq.append(JoinedSeq)
			else: #This should only be for the remaining bases in JoinedSeqRange:
				CmpdSeq.append('Gap')
		if Test == "yes":
			print CmpdSeq
	#Third, determining the sequence(s) for that locus.
	if len(CmpdSeq) == 0:
		#Just an error message if there were too few sequences to determine the final sequence
		print ("There were not enough data for locus %s for the sequence(s) to be determined.\n" % (LocusName))
		#sys.stderr.write("There were not enough data for locus %s for the sequence(s) to be determined.\n" % (LocusName))
		GenError = "Locus "+LocusName+" did not have enough sequence data for its sequence to be determined.\n"
		GenErrorList.append(GenError)
	elif len(AmbigDict) == 0:
		#Then there is only one sequence, and it is easy.
		NewSeq = ''
		for base in CmpdSeq:
			if base != 'Gap':
				NewSeq += base
		NewSeqDict[IndName] = NewSeq
		print ("There are %d sequences and one allele present at locus %s.\n" % (TotalSeqs, LocusName))
		#sys.stderr.write("There are %d sequences and one allele present at locus %s.\n" % (TotalSeqs, LocusName))
		SeqInfo = LocusName+"\t1\t"+str(TotalSeqs)+"\n"
		SeqInfoList.append(SeqInfo)
	elif len(AmbigDict) > 0:
		SumAlleleSeqs = 0
		AlleleList = [ ]
		#This is more difficult if there are multiple potential sequences.
		#In this case, it must be determined which SNPs go together.
		AmbigPosList = sorted(AmbigDict.keys())
		AmbigTempDict = { }
		for SeqName in SeqDict.keys():
			SeqTemp = ""
			for Position in AmbigPosList:
				SeqTemp += SeqDict[SeqName][Position]
			try:
				AmbigTempDict[SeqTemp] += SeqNumDict[SeqName]
			except KeyError:
				AmbigTempDict[SeqTemp] = SeqNumDict[SeqName]
		NumAlleles = 0
		if Test == "yes":
			print AmbigTempDict
		for Allele in AmbigTempDict:
			if AmbigTempDict[Allele] >= NumMin:
				NumAlleles += 1
				NewSeq = ""
				AmbigPos = 0
				for SeqPos in SeqRange:
					if (SeqPos in AmbigPosList) == False:
						if CmpdSeq[SeqPos] != "Gap":
							NewSeq += CmpdSeq[SeqPos]
					else:
						if Allele[AmbigPos] != "-":
							NewSeq += Allele[AmbigPos]
						AmbigPos += 1
				AlleleName = IndName+"_"+str(NumAlleles)
				NewSeqDict[AlleleName] = NewSeq
				SumAlleleSeqs += AmbigTempDict[Allele]
				AlleleList.append(str(AmbigTempDict[Allele]))
				print ("Allele %s is present %d times.\n" % (AlleleName, AmbigTempDict[Allele]))
				#sys.stderr.write("Allele %s is present %d times.\n" % (AlleleName, AmbigTempDict[Allele]))
		print ("There are %d sequences and %d alleles present at locus %s.\n" % (TotalSeqs, NumAlleles, LocusName))
		#sys.stderr.write("There are %d sequences and %d alleles present at locus %s.\n" % (TotalSeqs, NumAlleles, LocusName))
		SeqInfo = LocusName+"\t"+str(NumAlleles)+"\t"+','.join(AlleleList)+"\n"
		SeqInfoList.append(SeqInfo)
		if NumAlleles > Ploidy:
			print ("ATTENTION: there were %d alleles, but the plant is only %d-ploid!!\n" % (NumAlleles, Ploidy))
			#sys.stderr.write("ATTENTION: there were %d alleles, but the plant is only %d-ploid!!\n" % (NumAlleles, Ploidy))
			PloidyError = "Locus "+LocusName+" has "+str(NumAlleles)+" alleles, but its ploidy should be "+str(Ploidy)+".\n"
			PloidyErrorList.append(PloidyError)
		if SumAlleleSeqs < (TotalSeqs/2):
			print ("ATTENTION: there were a total of %d sequences at locus %s, but only %d were included in the final output.\n" % (TotalSeqs, LocusName, SumAlleleSeqs))
			sys.stderr.write("ATTENTION: there were a total of %d sequences at locus %s, but only %d were included in the final output.\n" % (TotalSeqs, LocusName, SumAlleleSeqs))
			GenError = "Locus "+LocusName+" had "+str(TotalSeqs)+" sequences, but only "+str(SumAlleleSeqs)+" were taken into account in the alleles output.\n"
			GenErrorList.append(GenError)
	if Test == "yes":
		print NewSeqDict
	#print the sequences for this locus to the file
	AlleleNameList = NewSeqDict.keys()
	AlleleNameList.sort()
	OutSeqList = [ ]
	for AlleleName in AlleleNameList:
		Line = ">"+AlleleName+"\n"+NewSeqDict[AlleleName]+"\n"
		OutSeqList.append(Line)
	OutFileName = SeqOutFileName
	if OutSeqList != [ ]:
		OutFile = open(OutFileName, 'a')
		for Line in OutSeqList:
			OutFile.write(Line)
		OutFile.close()
	

#print the other files with information about the loci
OutFileName = OutFolder+IndName+"_GenErrors.txt"
OutFile = open(OutFileName, 'w')
for Line in GenErrorList:
	OutFile.write(Line)
OutFile.close()
print ("General errors involving the sequences were written to the file %s.\n" % (OutFileName))
sys.stderr.write("General errors involving the sequences were written to the file %s.\n" % (OutFileName))

OutFileName = OutFolder+IndName+"_PloidyErrors.txt"
OutFile = open(OutFileName, 'w')
for Line in PloidyErrorList:
	OutFile.write(Line)
OutFile.close()
print ("Errors involving too many sequences given the ploidy of the plant were written to the file %s.\n" % (OutFileName))
sys.stderr.write("Errors involving too many sequences given the ploidy of the plant were written to the file %s.\n" % (OutFileName))

OutFileName = OutFolder+IndName+"_SeqInfo.txt"
OutFile = open(OutFileName, 'w')
for Line in SeqInfoList:
	OutFile.write(Line)
OutFile.close()
print ("The number of alleles recovered for each locus was written to the file %s.\n" % (OutFileName))
sys.stderr.write("The number of alleles recovered for each locus was written to the file %s.\n" % (OutFileName))
