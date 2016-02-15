#! /usr/bin/env python

#This script reads the alignment files made by seqparse.py and makes lists of the sequences
#and individuals in the alignments after filtering them to exclude uncommon sequences or
#uncommon individuals/populations.  It also makes a list of the sequences in the best 
#individual in each population and the best population in each taxon.
#Version 1.3 5 July 2015 Abby Moore
#Version 1.1 has been modified to remove the taxa and populations that did not have any
#"good" individuals, instead of leaving blank taxa/populations in the output.
#Version 1.2 has been modified to be able to give the option of removing the most variable x% of the loci
#Version 1.3 has been modified to only consider sequences from taxa/populations/individuals of interest
#The sequence names are of the form: >CA217_1_1
#> shows it is a sequence
#CA taxon code
#217 population number
#1 individual number
#1 sequence number (not always present)
'''
Example command:
alstats.py ~/Documents/finalseqs1/ ~/Documents/finalseqs1/filelisthead.txt 3 0.3 0.9 0.9 0.9 3 taxa ~/Documents/finalseqs1/SSTaxa_Calcifuge.txt ~/Documents/finalseqs1output_temp/
alstats.py ~/Documents/finalseqs1/ ~/Documents/finalseqs1/filelist.txt 3 0.1 0.7 0.7 0.7 3 none none ~/Documents/finalseqs1output_temp/
alstats.py ~/Documents/finalseqs1/ ~/Documents/finalseqs1/filelist.txt 1 0.1 0.7 0.7 0.7 3 taxa ~/Documents/finalseqs1/SSTaxa_LLC.txt ~/Documents/finalseqs1output_3LLC/
alstats.py InFolder SeqFileListFile MinNumTaxa MinPercSeqs PercTaxa PercPops PercInds AmtRemove SSType[none, taxa, pops, inds] SSListName OutFolder 
'''

import sys #We want to be able to send error messages to the screen.
from collections import defaultdict #We want to be able to make dictionaries with multiple levels.
from Bio import AlignIO #to read and parse sequence alignments
import math #We want to be able to take logs of things
import numpy #calculating means and percentiles
import matplotlib.pyplot #plotting data
import pandas #making data frames
import statsmodels.formula.api #doing regression

Usage = '''
alstats.py version 1.3 reads alignment files, makes lists of the plants in each 
alignment, and makes lists of the sequences in the alignments after filtering 
them according to various criteria.
Input should be in the following format:
alstats.py
[folder containing alignment files]
[file containing list of alignment files]
[minimum number of taxa that must have a sequence for it to be considered 
further]
[minimum proportion of the sequences that an individual must have for it to be 
analyzed]
[final proportion of taxa in which the sequences must occur]
[final proportion of populations in which the sequences must occur]
[final proportion of individuals in which the sequences must occur]
[distance above the median (in terms of number of interquartile ranges) at which
sequences will be considered outliers and removed]
[category to be used in making the subset of the data; options are taxa, pops, 
inds, and none]
[file containing the list of taxa/pops/inds to include, or none, if none]
[folder in which output files should be written] 
'''

SSTypeList = ['none', 'taxa', 'pops', 'inds']

print ("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 12:
	sys.exit("ERROR!! alstats.py requires 11 additional arguments and you supplied %d.\n %s" % (len(sys.argv)-1, Usage))
InFolder = sys.argv[1]
SeqFileListFile = sys.argv[2]
MinNumTaxa = int(sys.argv[3])
MinPercSeqs = float(sys.argv[4])
PercTaxa = float(sys.argv[5])
if (PercTaxa < 0 or PercTaxa > 1):
	sys.exit("ERROR!!  The proportion of taxa must be between 0 and 1 and you wrote %.3f.\n %s" % (PercTaxa, Usage))
PercPops = float(sys.argv[6])
if (PercPops < 0 or PercPops > 1):
	sys.exit("ERROR!!  The proportion of populations must be between 0 and 1 and you wrote %.3f.\n %s" % (PercPops, Usage))
PercInds = float(sys.argv[7])
if (PercInds < 0 or PercInds > 1):
	sys.exit("ERROR!!  The proportion of individuals must be between 0 and 1 and you wrote %.3f.\n %s" % (PercInds, Usage))
AmtRemove = float(sys.argv[8])
if (AmtRemove < 1):
	sys.exit("ERROR!!  The distance in terms of the interquartile range a locus has to be above the mean in order for it to be considered an outlier must be greater than 1 and you wrote %.3f.\n %s" % (AmtRemove, Usage))
SSType = sys.argv[9]
if (SSType in SSTypeList) == False:
	sys.exit("ERROR!!  You wanted the data to be subsetted by %s, but the only options are %s.\n%s" % (SSType, ",".join(SSTypeList), Usage))
SSListName = sys.argv[10]
if ((SSListName == "none") and (SSType != "none")):
	sys.exit("ERROR!!  You did not enter a list to show how the sequences should be divided, but you wanted them to be divided by %s.\n%s" % (SSType, Usage))
if ((SSType == "none") and (SSListName != "none")):
	print("WARNING!!  You did not specify the level at which the sequences should be divided, but you specified the file that should be used to divide them.  This file will be ignored.\n\n")
	sys.stderr.write("WARNING!!  You did not specify the level at which the sequences should be divided, but you specified the file that should be used to divide them.  This file will be ignored.\n\n")
OutFolder = sys.argv[11]

#making sure that the folder names end with a slash
if InFolder[-1] != "/":
	InFolder += "/"
print ("Sequence files will be read from the folder %s.\n" % (InFolder))
sys.stderr.write("Sequence files will be read from the folder %s.\n" % (InFolder))
if OutFolder[-1] != "/":
	OutFolder += "/"

#The main lists and dictionaries that will be filled out:
SeqFileList = []#The list of files of sequences (read from SeqFileListFile)
SSList = [ ]#The list of desired taxa/populations/individuals to be included
SeqsDict = defaultdict(dict) #defaultdict where sequences are the keys and individuals/taxa
#are the values
IndsDict = defaultdict(dict) #defaultdict where taxa/pops/inds are the keys and sequences are the values
SeqsDictN = defaultdict(dict) #same as SeqsDict but only with common seqs/inds
IndsDictN = defaultdict(dict) #same as IndsDict but only with commonseqs/inds
SeqVarDict = { } #SeqVarDict[SeqName] = CorVar = NumSeqs*PercVar/AlLength
VarList = [ ]#list of the variances
NumIS = defaultdict(int) #The dictionary with the number of individuals that have each locus.
NumPS = defaultdict(int) #The dictionary with the number of populations that have each locus.
NumTS = defaultdict(int) #The dictionary with the number of taxa that have each locus.
NumSI = defaultdict(int) #The dictionary with the number of loci each individual has.
NumSP = defaultdict(int) #The dictionary with the number of loci each population has.
#and the same for common sequences and individuals only:
NumISN = defaultdict(int)
NumPSN = defaultdict(int)
NumTSN = defaultdict(int)
NumSIN = defaultdict(int)
NumSPN = defaultdict(int)

#reading the list of sequence files
InFile = open(SeqFileListFile, 'rU')
for Line in InFile:
	Line = Line.strip('\r').strip('\n')
	SeqFileList.append(Line)
InFile.close()
print("%d files of sequences will be read.\n" % (len(SeqFileList)))
sys.stderr.write("%d files of sequences will be read.\n" % (len(SeqFileList)))

#reading the list of taxa/populations/individuals to be included, if necessary
if SSListName != "none":
	InFile = open(SSListName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n')
		SSList.append(Line)
	InFile.close()
	print("Only sequences from the %s included in the file %s will be output.\n" % (SSType, SSListName))
	sys.stderr.write("Only sequences from the %s included in the file %s will be output.\n" % (SSType, SSListName))
else:
	print("All sequences will be included in the output, as no file was given to tell which subset of the data to use.\n")
	sys.stderr.write("All sequences will be included in the output, as no file was given to tell which subset of the data to use.\n")

#Now to read each file of sequences and figure out which individuals are there:
OutList = [ ]
OutList2 = [ ]
#lambda function (!!) to calculate pairwise differences between sequences
diff2 = lambda seq1, seq2: sum(((seq1[i]!=seq2[i])and(seq1[i]in['A','C','G','T'])and(seq2[i]in['A','C','G','T'])) for i in range(len(seq1)))/float(len(seq1))
LocusList = [ ]
NumSeqsList = [ ]
PairwiseDiffList = [ ]
for SeqFile in SeqFileList:
	SeqFileName = InFolder+SeqFile
	SeqName = SeqFile[2:-3]
	try:
		Alignment = AlignIO.read(SeqFileName, 'fasta')
		NumSeqs = 0
		UsedSeqs = 0
		UnUsedSeqs = 0
		VarSites = 0.0
		#first, determining which taxa, populations, and individuals are present
		for Record in Alignment:
			NumSeqs += 1 #There must be a better way to calculate this!!!
			IndSeq = Record.id
			Line = IndSeq.split("_")
			TaxonName = Line[0][0:2]
			PopNum = Line[0]
			IndName = Line[0]+"_"+Line[1]
			if (SSType == "taxa") and ((TaxonName in SSList) == False):
				UnUsedSeqs += 1
			elif (SSType == "pops") and ((PopNum in SSList) == False):
				UnUsedSeqs += 1
			elif (SSType == "inds") and ((IndName in SSList) == False):
				UnUsedSeqs += 1
			#Putting all of this information into a defaultdict by sequence name.
			else:
				UsedSeqs += 1				
				try:
					SeqsDict[SeqName][TaxonName][PopNum][IndName].append(IndSeq)
				except KeyError:
					try:
						SeqsDict[SeqName][TaxonName][PopNum] = defaultdict(list)
						SeqsDict[SeqName][TaxonName][PopNum][IndName].append(IndSeq)
					except KeyError:
						SeqsDict[SeqName] = defaultdict(dict)
						SeqsDict[SeqName][TaxonName][PopNum] = defaultdict(list)
						SeqsDict[SeqName][TaxonName][PopNum][IndName].append(IndSeq)
				#Putting all of this information into a defaultdict by taxon, etc.
				try:
					IndsDict[TaxonName][PopNum][IndName].append(SeqName)
				except KeyError:
					try:
						IndsDict[TaxonName][PopNum] = defaultdict(list)
						IndsDict[TaxonName][PopNum][IndName] = [SeqName]
					except KeyError:
						IndsDict[TaxonName] = defaultdict(dict)
						IndsDict[TaxonName][PopNum] = defaultdict(list)
						IndsDict[TaxonName][PopNum][IndName] = [SeqName]
		print("Locus %s, out of %d sequences, %d were included and %d were not.\n" % (SeqName, NumSeqs, UsedSeqs, UnUsedSeqs))
		#*********It may be better to make a new alignment with the desired sequences??******
		#then, calculating sequence variability
		AlLength = Alignment.get_alignment_length()
		#Looking through each position in the alignment to see if it contains a SNP
		#only non-ambiguous nucleotides will be accepted.
		for SeqPos in range(AlLength):
			SeqPosList = list(Alignment[:,SeqPos])
			SeqPosUList = [ ]
			if 'A' in SeqPosList:
				SeqPosUList.append('A')
			if 'C' in SeqPosList:
				SeqPosUList.append('C')
			if 'G' in SeqPosList:
				SeqPosUList.append('G')
			if 'T' in SeqPosList:
				SeqPosUList.append('T')
			#including the site in the list of variable sites, if appropriate
			if len(SeqPosUList) > 1:
				VarSites += 1
		#Looking for the pairwise differences between the sequences
		DiffList = [ ]
		if NumSeqs > 1:
			LocusList.append(SeqName)
			NumSeqsList.append(NumSeqs)
			for Record1 in Alignment:
				DiffList2 = [ ]
				for Record2 in Alignment:
					if Record1.id != Record2.id:
						SeqDiff = diff2(str(Record1.seq),str(Record2.seq))
						DiffList.append(SeqDiff)
						DiffList2.append(SeqDiff)
				MeanDiff2 = numpy.mean(DiffList2)
			MeanDiff = numpy.mean(DiffList)
			PairwiseDiffList.append(MeanDiff)
		else:
			MeanDiff = 0.0
		matplotlib.pyplot.plot(NumSeqs,MeanDiff, '.', c='black')
		#print ("%s\t%d\t%.4f\n" % (SeqName, NumSeqs, MeanDiff))
	except IOError:
		"alas, no sequences in file"
		print("Sequence file %s not found.\n" % (SeqFileName))
print("%d files of sequences were successfully read.\n" % (len(SeqsDict)))
sys.stderr.write("%d files of sequences were successfully read.\n" % (len(SeqsDict)))

#looking through this information to find outliers
PairDiff_df = pandas.DataFrame({ 'NumSeqs' : NumSeqsList, 'PairwiseDiff' : PairwiseDiffList}, index=LocusList)
PairDiff_Model = statsmodels.formula.api.ols(formula='PairwiseDiff~NumSeqs', data = PairDiff_df).fit()
PairDiff_df['residuals'] = PairDiff_Model.resid
Slope = PairDiff_Model.params[1]
Intercept = PairDiff_Model.params[0]
PairDiff_df['pred_PD'] = Intercept+Slope*(PairDiff_df['NumSeqs'])
#PairDiff_df['pred_PD'] = PairDiff_Model.predict(PairDiff_df['NumSeqs'])#This gives me a ValueError and I have **no idea** why.
perc50 = numpy.percentile(PairDiff_df['residuals'], 50)
perc25 = numpy.percentile(PairDiff_df['residuals'], 25)
perc75 = numpy.percentile(PairDiff_df['residuals'], 75)
UpperBound = (perc75-perc25)*AmtRemove+perc50
print("The median of the data is %.4f, while the lower and upper quartiles are %.4f and %.4f, respectively.\n" % (perc50, perc25, perc75))
sys.stderr.write("The median of the data is %.4f, while the lower and upper quartiles are %.4f and %.4f, respectively.\n" % (perc50, perc25, perc75))
print("Any value with a residual greater than %.4f will be considered an outlier.\n" % (UpperBound))
sys.stderr.write("Any value with a residual greater than %.4f will be considered an outlier.\n" % (UpperBound))

OutFileName = OutFolder+"PairwiseDiff_AllSeqs.txt"
PairDiff_df.to_csv(OutFileName, sep="\t", float_format="%.5f")
print("The full data matrix was written to the file %s.\n" % (OutFileName))
sys.stderr.write("The full data matrix was written to the file %s.\n" % (OutFileName))

#adding a line going through these predicted values to the plot using matplotlib
matplotlib.pyplot.plot(PairDiff_df['NumSeqs'], PairDiff_df['pred_PD'], c='red', linewidth=2)

Outliers_df = PairDiff_df[PairDiff_df['residuals'] >= UpperBound]
matplotlib.pyplot.plot(Outliers_df['NumSeqs'], Outliers_df['PairwiseDiff'], '.', c='blue')
OutFileName = ("%sPairwiseDiff_Outliers_%.1f.txt" % (OutFolder, AmtRemove))
Outliers_df.to_csv(OutFileName, sep="\t", float_format="%.5f")
print("A total of %d outliers was found.\n" % (len(Outliers_df)))
sys.stderr.write("A total of %d outliers was found.\n" % (len(Outliers_df)))
print("The outliers were written to the file %s.\n" % (OutFileName))
sys.stderr.write("The outliers were written to the file %s.\n" % (OutFileName))
GoodSeqs = list(PairDiff_df[PairDiff_df['residuals'] < UpperBound].index)
print("%d sequences were not outliers.\n" % (len(GoodSeqs)))

#matplotlib.pyplot.show()
OutFileName = ("%sNumSeqs_PairDiff_%.1f.jpg" % (OutFolder, AmtRemove))
matplotlib.pyplot.savefig(OutFileName)
print("The plot of the number of sequences versus the average pairwise difference per locus was written to %s.\n" % (OutFileName))

#Making lists of sequences, populations, and taxa for the subsequent output
SeqsList = sorted(SeqsDict.keys())
TaxonList = sorted(IndsDict.keys())
PopsList = [ ]
for TaxonName in TaxonList:
	PopsList += IndsDict[TaxonName].keys()
PopsList = sorted(PopsList)
NumTaxa = len(TaxonList)
NumPops = len(PopsList)

print("These files contained sequences from %d populations of %d taxa.\n" % (len(PopsList), len(TaxonList)))
sys.stderr.write("These files contained sequences from %d populations of %d taxa.\n" % (len(PopsList), len(TaxonList)))

#Figuring out which sequences and individuals were present often enough to use:
for SeqName in SeqsList:
	for TaxonName in SeqsDict[SeqName].keys():
		NumTS[SeqName] += 1
		for PopNum in SeqsDict[SeqName][TaxonName].keys():
			NumPS[SeqName] += 1
			NumSP[PopNum] += 1
			for IndName in SeqsDict[SeqName][TaxonName][PopNum].keys():
				NumIS[SeqName] += 1
				NumSI[IndName] += 1

#Making a list of sequences that are present in more than three taxa.
CommonSeqs = [ ]
UnCommonSeqs = [ ]
for SeqName in SeqsList:
	if NumTS[SeqName] > MinNumTaxa:
		CommonSeqs.append(SeqName)
	else:
		UnCommonSeqs.append(SeqName)

print("Of these %d sequences, %d were present in more than %d taxa and will be considered further.\n" %\
	(len(SeqsList), len(CommonSeqs), MinNumTaxa))
sys.stderr.write("Of these %d sequences, %d were present in more than %d taxa and will be considered further.\n" %\
	(len(SeqsList), len(CommonSeqs), MinNumTaxa))

GoodInds = [ ]
BadInds = [ ]
for IndName in NumSI.keys():
	TestValue = float(NumSI[IndName])/len(SeqsList)
	if TestValue > MinPercSeqs:
		GoodInds.append(IndName)
	else:
		BadInds.append(IndName)
NumGInds = len(GoodInds)

print("Of these %d individuals, %d had at least %.2f proportion of the sequences and will be considered further.\n" %\
	 (len(NumSI.keys()), len(GoodInds), MinPercSeqs))
sys.stderr.write("Of these %d individuals, %d had at least %.2f proportion of the sequences and will be considered further.\n" %\
	 (len(NumSI.keys()), len(GoodInds), MinPercSeqs))

#Now to go back through the dictionaries and prune out the sequences and individuals that did not do well.
#First, make a new SeqsDict (SeqsDictN) that only includes the common sequences.
#I also filter this list for the sequences that are too variable.
for SeqName in CommonSeqs:
	if SeqName in GoodSeqs:
		for TaxonName in SeqsDict[SeqName].keys():
			for PopNum in SeqsDict[SeqName][TaxonName].keys():
				for IndName in SeqsDict[SeqName][TaxonName][PopNum].keys():
					if (IndName in GoodInds) == True:
						try:
							SeqsDictN[SeqName][TaxonName][PopNum][IndName] = SeqsDict[SeqName][TaxonName][PopNum][IndName]
						except KeyError:
							try:
								SeqsDictN[SeqName][TaxonName][PopNum] = defaultdict(list)
								SeqsDictN[SeqName][TaxonName][PopNum][IndName] = SeqsDict[SeqName][TaxonName][PopNum][IndName]
							except KeyError:
								SeqsDictN[SeqName] = defaultdict(dict)
								SeqsDictN[SeqName][TaxonName][PopNum] = defaultdict(list)
								SeqsDictN[SeqName][TaxonName][PopNum][IndName] = SeqsDict[SeqName][TaxonName][PopNum][IndName]
						#Putting all of this information into a defaultdict by taxon, etc.
						try:
							IndsDictN[TaxonName][PopNum][IndName].append(SeqName)
						except KeyError:
							try:
								IndsDictN[TaxonName][PopNum] = defaultdict(list)
								IndsDictN[TaxonName][PopNum][IndName] = [SeqName]
							except KeyError:
								IndsDictN[TaxonName] = defaultdict(dict)
								IndsDictN[TaxonName][PopNum] = defaultdict(list)
								IndsDictN[TaxonName][PopNum][IndName] = [SeqName]

#Make the new lists of taxa and populations that have at least one good individual. 
TaxonListN = sorted(IndsDictN.keys())
PopsListN = [ ]
for TaxonName in TaxonListN:
	PopsListN += IndsDictN[TaxonName].keys()
PopsListN = sorted(PopsListN)
NumTaxaN = len(TaxonListN)
NumPopsN = len(PopsListN)

print("Of the %d original taxa, %d have enough sequences to be considered further.\n" % (NumTaxa, NumTaxaN))
print("Of the %d original populations, %d have enough sequences to be considered further.\n" % (NumPops, NumPopsN))
sys.stderr.write("Of the %d original taxa, %d have enough sequences to be considered further.\n" % (NumTaxa, NumTaxaN))
sys.stderr.write("Of the %d original populations, %d have enough sequences to be considered further.\n" % (NumPops, NumPopsN))

#Then count how many taxa, populations, and individuals have each sequence.
SeqsListN = sorted(SeqsDictN.keys())
for SeqName in SeqsListN:
	for TaxonName in SeqsDictN[SeqName].keys():
		NumTSN[SeqName] += 1
		for PopNum in SeqsDictN[SeqName][TaxonName].keys():
			NumPSN[SeqName] += 1
			NumSPN[PopNum] += 1
			for IndName in SeqsDictN[SeqName][TaxonName][PopNum].keys():
				NumISN[SeqName] += 1
				NumSIN[IndName] += 1

#Then figure out what the minimum number of taxa, populations, and individuals are that have to have
#each sequence for it to be used further...
MinTaxa = NumTaxa*PercTaxa
SeqsTaxaF = [ ]
MinPops = NumPops*PercPops
SeqsPopsF = [ ]
MinInds = NumGInds*PercInds
SeqsIndsF = [ ]
#and make a list of the sequences that are present in enough taxa, populations, and individuals.
for SeqName in SeqsListN:
	if NumTSN[SeqName] > MinTaxa:
		SeqsTaxaF.append(SeqName)
	if NumPSN[SeqName] > MinPops:
		SeqsPopsF.append(SeqName)
	if NumISN[SeqName] > MinInds:
		SeqsIndsF.append(SeqName)

print("A sequence must be present in at least %.1f taxa to be included.  %d sequences meet this criterion.\n" \
	% (MinTaxa, len(SeqsTaxaF)))
print("A sequence must be present in at least %.1f populations to be included.  %d sequences meet this criterion.\n" \
	% (MinPops, len(SeqsPopsF)))
print("A sequence must be present in at least %.1f individuals to be included.  %d sequences meet this criterion.\n" \
	% (MinInds, len(SeqsIndsF)))

sys.stderr.write("A sequence must be present in at least %.1f taxa to be included.  %d sequences meet this criterion.\n" \
	% (MinTaxa, len(SeqsTaxaF)))
sys.stderr.write("A sequence must be present in at least %.1f populations to be included.  %d sequences meet this criterion.\n" \
	% (MinPops, len(SeqsPopsF)))
sys.stderr.write("A sequence must be present in at least %.1f individuals to be included.  %d sequences meet this criterion.\n" \
	% (MinInds, len(SeqsIndsF)))

#Making the output files.  First, these are the output files for all of the sequences, taxa, and populations:
#First, the number of populations of each taxon that have a given sequence.
TaxOut = [ ]
Head = "Taxon\t"+"\t".join(SeqsList)
TaxOut.append(Head)
for TaxonName in TaxonList:
	Line = [TaxonName]
	for SeqName in SeqsList:
		try:
			PopList = SeqsDict[SeqName][TaxonName].keys()
			NumPops = str(len(PopList))
			Line.append(NumPops)
		except KeyError:
			Line.append('0')
	TaxOut.append('\n')
	TaxOut.append('\t'.join(Line))

#Second, the number of individuals in each population that have a given sequence.
PopsOut = [ ]
Head = "Population\t"+"\t".join(SeqsList)
PopsOut.append(Head)
for PopNum in PopsList:
	TaxonName = PopNum[0:2]
	Line = [PopNum]
	for SeqName in SeqsList:
		try:
			IndList = SeqsDict[SeqName][TaxonName][PopNum].keys()
			NumInds = str(len(IndList))
			Line.append(NumInds)
		except KeyError:
			Line.append('0')
	PopsOut.append('\n')
	PopsOut.append('\t'.join(Line))

#Writing the output files:
#First, the file by taxon:
OutFileName = OutFolder+"tseq_counts.txt"
OutFile = open(OutFileName, 'w')
for Line in TaxOut:
	OutFile.write(Line)
OutFile.close()
print("The full table of the sequences present in each taxon was written to %s.\n" % (OutFileName))
sys.stderr.write("The full table of the sequences present in each taxon was written to %s.\n" % (OutFileName))

#Second, the file by population:
OutFileName = OutFolder+"pseq_counts.txt"
OutFile = open(OutFileName, 'w')
for Line in PopsOut:
	OutFile.write(Line)
OutFile.close()
print("The full table of the sequences present in each population was written to %s.\n" % (OutFileName))
sys.stderr.write("The full table of the sequences present in each population was written to %s.\n" % (OutFileName))

#Now, these are the same output files, but for the pruned set of sequences and populations.
#First, the number of populations of each taxon that have a given sequence.
TaxOut = [ ]
CommonSeqs = sorted(CommonSeqs)
Head = "Taxon\t"+"\t".join(CommonSeqs)
TaxOut.append(Head)
for TaxonName in TaxonListN:
	Line = [TaxonName]
	for SeqName in CommonSeqs:
		try:
			PopList = SeqsDictN[SeqName][TaxonName].keys()
			NumPops = str(len(PopList))
			Line.append(NumPops)
		except KeyError:
			Line.append('0')
	TaxOut.append('\n')
	TaxOut.append('\t'.join(Line))

#Second, the number of individuals in each population that have a given sequence.
PopsOut = [ ]
Head = "Population\t"+"\t".join(CommonSeqs)
PopsOut.append(Head)
for PopNum in PopsListN:
	TaxonName = PopNum[0:2]
	Line = [PopNum]
	for SeqName in CommonSeqs:
		try:
			IndList = SeqsDictN[SeqName][TaxonName][PopNum].keys()
			NumInds = str(len(IndList))
			Line.append(NumInds)
		except KeyError:
			Line.append('0')
	PopsOut.append('\n')
	PopsOut.append('\t'.join(Line))

#Writing the output files:
#First, the file by taxon:
OutFileName = OutFolder+"common_tseq_counts.txt"
OutFile = open(OutFileName, 'w')
for Line in TaxOut:
	OutFile.write(Line)
OutFile.close()
print("The table of the sequences present in each taxon that only includes sequences present in at least three taxa was written to %s.\n" % (OutFileName))
sys.stderr.write("The table of the sequences present in each taxon that only includes sequences present in at least three taxa was written to %s.\n" % (OutFileName))

#Second, the file by population:
OutFileName = OutFolder+"common_pseq_counts.txt"
OutFile = open(OutFileName, 'w')
for Line in PopsOut:
	OutFile.write(Line)
OutFile.close()
print("The table of the sequences present in each population that only includes sequences present in at least three taxa was written to %s.\n" % (OutFileName))
sys.stderr.write("The table of the sequences present in each population that only includes sequences present in at least three taxa was written to %s.\n" % (OutFileName))

#Now we need to make output files for the actual lists of approved taxa, populations, and individuals, so that those can be used further.
#I think everything will take the format of a tab-delimited file, with the rows as loci.
#The first column will be the locus name.
#The remaining columns will be the sequences in each taxon/pop/ind, with multiple sequences separated by commas.
#Taxa first:
TaxaOut = [ ]
Head = "Locus\t"+"\t".join(TaxonListN)
TaxaOut.append(Head)
SeqsTaxaF = sorted(SeqsTaxaF)
for SeqName in SeqsTaxaF:
	Line = "\n"+SeqName
	for TaxonName in TaxonListN:
		TSeqList = [ ]
		try:
			for PopNum in SeqsDictN[SeqName][TaxonName].keys():
				for IndName in SeqsDictN[SeqName][TaxonName][PopNum].keys():
					TSeqList += SeqsDictN[SeqName][TaxonName][PopNum][IndName]
		except KeyError:
			"do nothing"
		Line += "\t"+",".join(TSeqList)
	TaxaOut.append(Line)

#Now populations:
PopsOut = [ ]
Head = "Locus\t"+"\t".join(PopsListN)
PopsOut.append(Head)
SeqsPopsF = sorted(SeqsPopsF)
for SeqName in SeqsPopsF:
	Line = "\n"+SeqName
	for PopNum in PopsListN:
		PSeqList = [ ]
		TaxonName = PopNum[0:2]
		try:
			for IndName in SeqsDictN[SeqName][TaxonName][PopNum].keys():
				PSeqList += SeqsDictN[SeqName][TaxonName][PopNum][IndName]
		except KeyError:
			"do nothing"
		Line += "\t"+",".join(PSeqList)
	PopsOut.append(Line)

#Now individuals:
IndsOut = [ ]
GoodInds = sorted(GoodInds)
SeqsIndsF = sorted(SeqsIndsF)
Head = "Locus\t"+"\t".join(GoodInds)
IndsOut.append(Head)
for SeqName in SeqsIndsF:
	Line = "\n"+SeqName
	for IndName in GoodInds:
		PopNum = IndName.split("_")[0]
		TaxonName = PopNum[0:2]
		try:
			ISeqList = SeqsDictN[SeqName][TaxonName][PopNum][IndName]
		except KeyError:
			ISeqList = [ ]
		Line += "\t"+",".join(ISeqList)
	IndsOut.append(Line)

#Now to print everything:
OutFileName = OutFolder+"seqs_taxa.txt"
OutFile = open(OutFileName, 'w')
for Line in TaxaOut:
	OutFile.write(Line)
OutFile.close()
print("The names of the sequences of each locus that are present in each taxon were written to %s.\n" % (OutFileName))
sys.stderr.write("The names of the sequences of each locus that are present in each taxon were written to %s.\n" % (OutFileName))

OutFileName = OutFolder+"seqs_pops.txt"
OutFile = open(OutFileName, 'w')
for Line in PopsOut:
	OutFile.write(Line)
OutFile.close()
print("The names of the sequences of each locus that are present in each population were written to %s.\n" % (OutFileName))
sys.stderr.write("The names of the sequences of each locus that are present in each population were written to %s.\n" % (OutFileName))

OutFileName = OutFolder+"seqs_inds.txt"
OutFile = open(OutFileName, 'w')
for Line in IndsOut:
	OutFile.write(Line)
OutFile.close()
print("The names of the sequences of each locus that are present in each individual were written to %s.\n" % (OutFileName))
sys.stderr.write("The names of the sequences of each locus that are present in each individual were written to %s.\n" % (OutFileName))

#Now making lists of the population with the most sequences within each taxon and the individual with the most sequences within each
#population, in case we want to run subsets of the data.

TaxonPopDict = { }
for TaxonName in TaxonList:
	MostSeqs = 0
	for PopNum in IndsDict[TaxonName].keys():
		try:
			if NumSPN[PopNum] > MostSeqs:
				BestPop = PopNum
				MostSeqs = NumSPN[PopNum]
		except KeyError:
			"no pop"
	TaxonPopDict[TaxonName] = [BestPop, str(MostSeqs)]

PopIndDict = { }
for PopNum in PopsList:
	TaxonName = PopNum[0:2]
	MostSeqs = 0
	for IndName in IndsDict[TaxonName][PopNum].keys():
		try:
			if NumSIN[IndName] > MostSeqs:
				BestInd = IndName
				MostSeqs = NumSIN[IndName]
		except KeyError:
			"no ind"
	PopIndDict[PopNum] = [BestInd, str(MostSeqs)]

#Now writing those lists to files:

TPList = [ ]
Head = "Taxon\tPopulation\tNumber_of_Sequences\n"
TPList.append(Head)
for TaxonName in TaxonList:
	Line = TaxonName+"\t"+"\t".join(TaxonPopDict[TaxonName])+"\n"
	TPList.append(Line)

PIList = [ ]
Head = "Population\tIndividual\tNumber_of_Sequences\n"
PIList.append(Head)
for PopNum in PopsList:
	Line = PopNum+"\t"+"\t".join(PopIndDict[PopNum])+"\n"
	PIList.append(Line)

OutFileName = OutFolder+"pops_taxa.txt"
OutFile = open(OutFileName, 'w')
for Line in TPList:
	OutFile.write(Line)
OutFile.close()
print("The list of the populations that had the most sequences in each taxon was written to %s.\n" % (OutFileName))
sys.stderr.write("The list of the populations that had the most sequences in each taxon was written to %s.\n" % (OutFileName))

OutFileName = OutFolder+"inds_pops.txt"
OutFile = open(OutFileName, 'w')
for Line in PIList:
	OutFile.write(Line)
OutFile.close()
print("The list of the individuals that had the most sequences in each population was written to %s.\n" % (OutFileName))
sys.stderr.write("The list of the individuals that had the most sequences in each population was written to %s.\n" % (OutFileName))

#We also need lists of the sequences of those individuals and populations, so they can be used
#to make pruned alignments for downstream analysis.
#First, making the list of sequences from the best population in each taxon.
#Since we are using a single population as a representative of the taxon, I think it is
#alright to combine individuals, so that if a sequence is not present in one individual,
#it can be taken from another individual.
BTPList = [ ]
BestPops =[ ]
Head = "Locus"
for TaxonName in TaxonListN:
	Head += "\t"+TaxonPopDict[TaxonName][0]
	BestPops.append(TaxonPopDict[TaxonName][0])
BTPList.append(Head)
for SeqName in SeqsTaxaF:
	Line = "\n"+SeqName
	for PopNum in BestPops:
		TSeqList = [ ]
		TaxonName = PopNum[0:2]
		try:
			for IndName in SeqsDictN[SeqName][TaxonName][PopNum].keys():
				TSeqList = SeqsDictN[SeqName][TaxonName][PopNum][IndName]
		except KeyError:
			"do nothing"
		Line += "\t"+",".join(TSeqList)
	BTPList.append(Line)

#Now the same thing for the best individual from each population
BPIList = [ ]
BestInds =[ ]
Head = "Locus"
for PopNum in PopsListN:
	Head += "\t"+PopIndDict[PopNum][0]
	BestInds.append(PopIndDict[PopNum][0])
BPIList.append(Head)
for SeqName in SeqsPopsF:
	Line = "\n"+SeqName
	for IndName in BestInds:
		PSeqList = [ ]
		PopNum = IndName.split("_")[0]
		TaxonName = PopNum[0:2]
		try:
			PSeqList = SeqsDictN[SeqName][TaxonName][PopNum][IndName]
		except KeyError:
			"do nothing"
		Line += "\t"+",".join(PSeqList)
	BPIList.append(Line)

#Now to write this to files:

OutFileName = OutFolder + "seqs_taxa_best.txt"
OutFile = open(OutFileName, 'w')
for Line in BTPList:
	OutFile.write(Line)
OutFile.close()
print("The names of the sequences of each locus from the best population from each taxon were written to %s.\n" % (OutFileName))
sys.stderr.write("The names of the sequences of each locus from the best population from each taxon were written to %s.\n" % (OutFileName))

OutFileName = OutFolder + "seqs_pops_best.txt"
OutFile = open(OutFileName, 'w')
for Line in BPIList:
	OutFile.write(Line)
OutFile.close()
print("The names of the sequences of each locus from the best individual from each population were written to %s.\n" % (OutFileName))
sys.stderr.write("The names of the sequences of each locus from the best individual from each population were written to %s.\n" % (OutFileName))

