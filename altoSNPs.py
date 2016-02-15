#! /usr/bin/env python

#altoSNPs version 2.1 6 Feb. 2016 Abby Moore
#This script reads alignments in fasta format and figures out which sites are SNPs.
#It can then write files of the SNPs in various different formats for different
#types of analyses.
#24 Nov. 2015: adding adegenet output format
#version 2.1, 6 Feb. 2016, added fasta output format, allowed SNPs to be output in the same order as the KeyList


#structure output format:
'''
first row: list of locus names
second row (if desired): list of recessive alleles
third row (if desired): inter-marker distances (maybe useful if using multiple SNPs per locus?)
Then information for the individuals:
first column: individual name
second column (if desired): population/taxon
third column (if desired): sampling location/population
remaining columns: genotype data (I think this has to be integers, with -9 as missing)
'''
#fastStructure output format:
'''
no header row
first column: individual name
second column: taxon name (or random number)
third column: population name (or random number)
fourth column: random number
fifth column: random number
sixth column: random number
remaining column: genotype data (integers, with -9 as missing)
individuals in two rows
'''
#GenePop output format:
'''
first row: random title
second row: list of loci, separated by commas
third row (and all rows before another population): POP
fourth row (and all rows for each individual): indname, 0101 0101 0101 (with each allele having a two- or three-digit code and alleles per locus in pairs, separated by commas)
'''
#adegenet output format:
'''
as for structure, but either taxon or population and these written out, not numbers
also .r file for analysis
'''
#fasta output format:
'''
Each individual has two lines IndName_1 and IndName_2.
There is a separate file with the name of each SNP on a separate line.
'''

'''
example command:
altoSNPs.py ~/Documents/finalseqs2_pops_best/ 2 ~/Documents/finalseqs1output/seqs_pops_best.txt ~/Documents/seqs_formatted/ trial structure all 2
altoSNPs.py ~/Documents/finalseqs2/ 2 ~/Documents/finalseqs1output/seqs_taxa_best.txt ~/Documents/seqs_formatted/ trial structure all 1
altoSNPs.py ~/Documents/finalseqs2_3LLC/inds3LLC/ 3LLC ~/Documents/finalseqs1output_3LLC/seqs_inds.txt ~/Documents/seqs_formatted_outliers/adegenet20151222/ allperpop_3LLC adegenet all 3
altoSNPs.py InFolder InFilePre AlListIn OutFolder OutFilePre OutFormat[structure, fastStructure, GenePop, adegenet] SNPsperLocus DataLevels
'''

import sys #We need to talk to the command line.
from collections import defaultdict #We want to be able to make dictionaries with multiple levels
from Bio import AlignIO #We want to be able to read alignments
import random #We want to be able to generate (pseudo-)random numbers.

Usage = '''
altoSNPs version 2.1
This script reads alignments created by newal.py and finds SNPs.
Then it can output the SNPs in various formats for further analysis.
The available formats are: structure, fastStructure, GenePop, adegenet, fasta.
When calling the program, the output format must be written exactly as it is
in the previous sentence (including capitalization).
altoSNPs
[folder containing the alignments]
[prefix for the alignments]
[alstats.py file that was used to make the alignments]
[output folder] 
[prefix for output file]
[output format]
[# SNPs per locus: either "all" or an integer]
[number of levels in output data: 1 if there is one individual per taxon, 2 if there is one
individual per population, 3 if there are multiple individuals per population]
'''

print ("%s\n" % (" ".join(sys.argv)))

OutFormatList = [ 'structure', 'fastStructure', 'GenePop', 'adegenet', 'fasta'] 

if len(sys.argv) != 9:
	sys.exit("ERROR! There were too few input parameters.  Eight additional arguments are required, and you supplied %d.\n  %s" % (len(sys.argv)-1, Usage))
InFolder = sys.argv[1]
InFilePre = sys.argv[2]
AlListIn = sys.argv[3]
OutFolder = sys.argv[4]
OutFilePre = sys.argv[5]
OutFormat = sys.argv[6]
if (OutFormat in OutFormatList) == False:
	sys.exit("ERROR!  The output format can be one of the following options:\n%s\n.  But you wrote %s." % (" ".join(OutFormatList), OutFormat))
SNPsperLocus = sys.argv[7]
DataLevels = int(sys.argv[8])

if (DataLevels != 1) and (DataLevels != 2) and (DataLevels != 3):
	sys.exit("ERROR!  The number of levels can only be 1, 2, or 3!!!")

#making sure that the folder names end with a slash
if InFolder[-1] != "/":
	InFolder += "/"
print ("Sequence files will be read from the folder %s.\n" % (InFolder))
sys.stderr.write("Sequence files will be read from the folder %s.\n" % (InFolder))
if OutFolder[-1] != "/":
	OutFolder += "/"

SNPDict = defaultdict(dict) #The dictionary of SNPs that is ready to be formatted
#It will be SNPDict[SeqName][SeqPos][OTU] = list of SNPs
NumSeqs = 0 #Counter for the number of sequences.
NumSNPs = 0 #Counter for the number of SNPs
InFileDict = { } #Dictionary of sequence files to be read and their corresponding locus names.
OTUDict = { } #Dictionary of the sequence names and their corresponding OTUs
OTUList = [ ] #The list of the various OTUs
OTUSeqDict = defaultdict(dict) #The dictionary with lists of sequences for each OTU
MaxSNPs = 0 #The maximum number of SNPs per site per individual

#First, we need to read the names of the loci from the file that was used to create
#the alignments in the first place.
InFile = open(AlListIn, 'rU')
#KeyLine = [ ] #This will be the header line, so we can keep straight what the OTU names are.
for Line in InFile:
	Line = Line.strip('\n').strip('\r').split('\t')
	#Making a dictionary of locus names and files containing those loci
	if Line[0] == "Locus":
		#making a key to the OTUs so we can associate the sequence names with the OTU as we read the file
		KeyLine = Line
		OTUList = Line[1:]
		NumOTUs = len(OTUList)
	else:
		#adding the locus name to the dictionary of loci and file names
		InFileDict[InFolder+InFilePre+Line[0]+".fa"] = Line[0]
		#adding the sequence names (OTUSeq) to the dictionary that tells which sequence name goes to which OTU
		for OTUNum in range(1,NumOTUs+1):
			SeqListTemp = Line[OTUNum].split(",")
			for OTUSeq in SeqListTemp:
				if OTUSeq != '':
					OTUDict[OTUSeq] = KeyLine[OTUNum]
InFile.close()
print("The list of %d alignments was read from the file %s.\n" % (len(InFileDict), AlListIn))
sys.stderr.write("The list of %d alignments was read from the file %s.\n" % (len(InFileDict), AlListIn))
print("The alignments contained sequences from %d OTUs.\n" % (OTUNum))
sys.stderr.write("The alignments contained sequences from %d OTUs.\n" % (OTUNum))

#Making dictionaries of which individuals belong to which taxa and populations
TaxIndDict = defaultdict(list)
PopIndDict = defaultdict(list)
for OTU in OTUList:
	Taxon = OTU[0:2]
	Pop = OTU.split("_")[0]
	TaxIndDict[Taxon].append(OTU)
	PopIndDict[Pop].append(OTU)

#Then, we need to read each alignment and find the SNPs.
GoodNucs = ['A','C','G','T'] #The list of nucleotides we want in the end.
for InFileName in InFileDict:
	Alignment = AlignIO.read(InFileName, "fasta")
	AlLength = Alignment.get_alignment_length()
	SeqName = InFileDict[InFileName]
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
		#If so, adding that SNP to the dictionary
		#It will be SNPDict[SeqName][SeqPos][OTU] = list of SNPs
		if len(SeqPosUList) > 1:
			for Record in Alignment:
				OTU = OTUDict[Record.id]
				PosNuc = Record.seq[SeqPos]
				#Replacing ambiguities (N, W, etc.) with -9
				if (PosNuc in GoodNucs) == False:
					PosNuc = '-9'
			#Adding the appropriate nucleotides to the dictionary.
				try:
					SNPDict[SeqName][SeqPos][OTU].append(PosNuc)
				except KeyError:
					SNPDict[SeqName][SeqPos] = defaultdict(list)
					SNPDict[SeqName][SeqPos][OTU] = [PosNuc]
		#getting rid of nucleotides that are listed twice at the same position (because that OTU had two or more
		#sequences for that locus)
			for OTU in SNPDict[SeqName][SeqPos].keys():
				SNPDict[SeqName][SeqPos][OTU] = list(set(SNPDict[SeqName][SeqPos][OTU]))
				if '-9' in SNPDict[SeqName][SeqPos][OTU]:
					if len(SNPDict[SeqName][SeqPos][OTU]) > 1:
						SNPDict[SeqName][SeqPos][OTU].remove('-9')
				TempMaxSNPs = len(SNPDict[SeqName][SeqPos][OTU])
				if TempMaxSNPs > MaxSNPs:
					MaxSNPs = TempMaxSNPs

print("SNPs were read from the alignment files, with names of the form %s.\n" % (InFileName))
print("A maximum of %d nucleotides was found at a single site in a single individual.\n" % (MaxSNPs))
sys.stderr.write("SNPs were read from the alignment files, with names of the form %s.\n" % (InFileName))
sys.stderr.write("A maximum of %d nucleotides was found at a single site in a single individual.\n" % (MaxSNPs))


#Finally, the SNPs need to be formatted properly in the desired output format.
#first, preparing the things that will be filled out:
KeyList = [ ]
#Then parsing the SNPDict:
SeqNameList = sorted(SNPDict.keys())
for SeqName in SeqNameList:
	SeqPosList = sorted(SNPDict[SeqName].keys())
	#Adding all of the SNPs per locus if all should be added.
	if SNPsperLocus == "all":
		NewSeqPosList = SeqPosList
	#If not, randomly choosing the desired number.
	else:
		SNPsperLocus = int(SNPsperLocus)
		if len(SeqPosList) > SNPsperLocus:
			NewSeqPosList = random.sample(SeqPosList, SNPsperLocus)
			NewSeqPosList = sorted(NewSeqPosList)
		else:
			NewSeqPosList = SeqPosList
	#Now, going through each SNP in that locus:
	for SeqPos in NewSeqPosList:
		SNPName = SeqName+"_"+str(SeqPos)
		KeyList.append(SNPName)
		for OTU in OTUList:
			NumPoly = len(SNPDict[SeqName][SeqPos][OTU])
			OTUSeqDict[OTU][SNPName] = defaultdict(dict)
			if NumPoly == 1:
				for Num in range(MaxSNPs):
					OTUSeqDict[OTU][SNPName][Num] = SNPDict[SeqName][SeqPos][OTU][0]
			elif NumPoly == MaxSNPs:
				for Num in range(MaxSNPs):
					OTUSeqDict[OTU][SNPName][Num] = SNPDict[SeqName][SeqPos][OTU][Num]
			elif NumPoly == 0:
				#insert the value for missing data if there are no SNPs for that OTU
				for Num in range(MaxSNPs):
					OTUSeqDict[OTU][SNPName][Num] = '-9'
			elif (NumPoly > 1) and (NumPoly < MaxSNPs):
				#if the number of different nucleotides for that locus is more than one but less
				#than the ploidy, add the nucleotides that are present and randomly fill out the
				#rest.
				for Num in range(MaxSNPs):
					try:
						OTUSeqDict[OTU][SNPName][Num] = SNPDict[SeqName][SeqPos][OTU][Num]
					except IndexError:
						OTUSeqDict[OTU][SNPName][Num] = random.sample(SNPDict[SeqName][SeqPos][OTU], 1)[0]
KeyList = sorted(KeyList)

print("%d SNPs were found in these alignments.\n" % (len(OTUSeqDict[OTU].keys())))
sys.stderr.write("%d SNPs were found in these alignments.\n" % (len(OTUSeqDict[OTU].keys())))

if OutFormat == "structure":
	SNPListOut = [ ]
	MainParamsList = [ ]
	OutFileName = OutFolder+OutFilePre+"structure.txt"
	OutFileNameMP = OutFolder + OutFilePre + "mainparams"
	OutFileNameStr = OutFolder+OutFilePre+"out"
	#If necessary, convert taxon and population to numbers
	if DataLevels > 1:
		TaxonList = [ ]
		TaxDict = { }#Dictionary saying which # belongs to which taxon
		PopList = [ ]
		PopDict = { }#Dictionary saying which # belongs to which population
		for OTU in OTUList:
			Taxon = OTU[0:2]
			TaxonList.append(Taxon)
			Pop = OTU.split("_")[0]
			PopList.append(Pop)
		TaxonList = list(set(TaxonList))
		TaxonList = sorted(TaxonList)
		PopList = list(set(PopList))
		PopList = sorted(PopList)
		TaxNum = 1
		PopNum = 1
		for Taxon in TaxonList:
			TaxDict[Taxon] = TaxNum
			TaxNum += 1
		for Pop in PopList:
			PopDict[Pop] = PopNum
			PopNum += 1
	#making the header for the main output file
	Line = ("\t").join(KeyList)+"\n"
	SNPListOut.append(Line)
	#making the body of the main output file
	for OTU in OTUList:
		if DataLevels == 1:
			IndName = OTU+"\t"
		elif DataLevels == 2:
			IndName = OTU+"\t"+str(TaxDict[OTU[0:2]])+"\t"
		elif DataLevels == 3:
			IndName = OTU+"\t"+str(TaxDict[OTU[0:2]])+"\t"+str(PopDict[OTU.split("_")[0]])+"\t"
		for Num in range(MaxSNPs):
			TempString = ""
			for SNPName in KeyList:
				TempString += OTUSeqDict[OTU][SNPName][Num]+"\t"
			#remove the last tab character
			TempString = TempString[:-1]
			TempString = TempString.replace('A','0')
			TempString = TempString.replace('C','1')
			TempString = TempString.replace('G','2')
			TempString = TempString.replace('T','3')
			Line = IndName+TempString+"\n"
			SNPListOut.append(Line)
	#making the parameter file
	Line = "#define MAXPOPS 3\n#define BURNIN 50000\n#define NUMREPS 100000\n"
	MainParamsList.append(Line)
	Line = "#define INFILE "+OutFileName+"\n#define OUTFILE "+OutFileNameStr+"\n"
	MainParamsList.append(Line)
	Line = "#define NUMINDS\t"+str(len(OTUList))+"\n#define NUMLOCI "+str(len(OTUSeqDict[OTU].keys()))+"\n"
	MainParamsList.append(Line)
	Line = "#define PLOIDY "+str(MaxSNPs)+"\n#define MISSING -9\n#define ONEROWPERIND 0\n#define LABEL 1\n"
	MainParamsList.append(Line)
	if DataLevels == 1:
		Line = "#define POPDATA 0\n#define POPFLAG 0\n#define LOCDATA 0\n"
	elif DataLevels == 2:
		Line = "#define POPDATA 1\n#define POPFLAG 0\n#define LOCDATA 0\n"
	elif DataLevels == 3:
		Line = "#define POPDATA 1\n#define POPFLAG 0\n#define LOCDATA 1\n"
	MainParamsList.append(Line)
	Line = "#define PHENOTYPE 0\n#define EXTRACOLS 0\n#define MARKERNAMES 1\n#define RECESSIVEALLELES 0\n"
	MainParamsList.append(Line)
	Line = "#define MAPDISTANCES 0\n#define PHASED 0\n#define PHASEINFO 0\n"
	MainParamsList.append(Line)
	OutFile = open(OutFileNameMP, 'w')
	for Line in MainParamsList:
		OutFile.write(Line)
	OutFile.close()
	print("In order to run the Structure file, you will need to go to the folder where Structure is found, \
	and type the following:./structure -m %s -K [number of populations you want to test] \n" % (OutFileNameMP))
	sys.stderr.write("In order to run the Structure file, you will need to go to the folder where Structure is found, \
	and type the following:./structure -m %s -K [number of populations you want to test] \n" % (OutFileNameMP))
	#writing keys to the taxon and population numbers, if necessary
	if DataLevels > 1:
		OutTaxList = [ ]
		for Taxon in TaxDict:
			Line = Taxon+"\t"+str(TaxDict[Taxon])+"\n"
			OutTaxList.append(Line)
		OutFileNameTL = OutFolder+OutFilePre+"Taxon_Key.txt"
		OutFile = open(OutFileNameTL, 'w')
		for Line in OutTaxList:
			OutFile.write(Line)
		OutFile.close()
		print("The list of numbers corresponding to the various taxa was written to the file %s.\n" % (OutFileNameTL))
		sys.stderr.write("The list of numbers corresponding to the various taxa was written to the file %s.\n" % (OutFileNameTL))
	if DataLevels > 2:
		OutPopList = [ ]
		for Pop in PopDict:
			Line = Pop+"\t"+str(PopDict[Pop])+"\n"
			OutPopList.append(Line)
		OutFileNamePL = OutFolder+OutFilePre+"Pop_Key.txt"
		OutFile = open(OutFileNamePL, 'w')
		for Line in OutPopList:
			OutFile.write(Line)
		OutFile.close()
		print("The list of numbers corresponding to the various populations was written to the file %s.\n" % (OutFileNamePL))
		sys.stderr.write("The list of numbers corresponding to the various populations was written to the file %s.\n" % (OutFileNamePL))
elif OutFormat == "fastStructure":
	#lists of taxa and populations for the labels file
	TaxList = [ ]
	PopList = [ ]
	#what we will actually print
	SNPListOut = [ ]
	OutFileName = OutFolder+OutFilePre+"fstructure.str"
	#no header for fastStructure
	#making the body of the main output file
	for OTU in OTUList:
		if DataLevels == 1:
			IndName = OTU+"\t1\t2\t3\t4\t5\t"
		elif DataLevels == 2:
			#IndName = OTU+"\t"+str(TaxDict[OTU[0:2]])+"\t1\t2\t3\t4\t"
			IndName = OTU+"\t"+OTU[0:2]+"\t1\t2\t3\t4\t"
		elif DataLevels == 3:
			#IndName = OTU+"\t"+str(TaxDict[OTU[0:2]])+"\t"+str(PopDict[OTU.split("_")[0]])+"\t1\t2\t3\t"
			IndName = OTU+"\t"+OTU[0:2]+"\t"+OTU.split("_")[0]+"\t1\t2\t3\t"
		for Num in range(MaxSNPs):
			TempString = ""
			for SNPName in KeyList:
				TempString += OTUSeqDict[OTU][SNPName][Num]+"\t"
			#remove the last tab character
			TempString = TempString[:-1]
			TempString = TempString.replace('A','0')
			TempString = TempString.replace('C','1')
			TempString = TempString.replace('G','2')
			TempString = TempString.replace('T','3')
			Line = IndName+TempString+"\n"
			SNPListOut.append(Line)
		TaxList.append(OTU[0:2])
		PopList.append(OTU.split("_")[0])
	#There is no parameter file, either.
	print("In order to run fastStructure, you will need to copy the output file to a subdirectory of the directory in which fastStructure is found.\
	The commands will be of the following form: python structure.py -K 2 --input=dir/%s --output=dir/outfilehead --format=str\
	[if desired] --prior=logistic\n" % (OutFileName))
	sys.stderr.write("In order to run fastStructure, you will need to copy the output file to a subdirectory of the directory in which fastStructure is found.\
	The commands will be of the following form: python structure.py -K 2 --input=dir/%s --output=dir/outfilehead --format=str\
	[if desired] --prior=logistic\n" % (OutFileName))
	#writing keys to the taxon and population numbers, if necessary
	if DataLevels > 1:
		OutFileNameTL = OutFolder+OutFilePre+"Taxon_Key.txt"
		OutFile = open(OutFileNameTL, 'w')
		OutFile.write("\n".join(TaxList))
		OutFile.close()
		print("The list of numbers corresponding to the various taxa was written to the file %s.\n" % (OutFileNameTL))
		sys.stderr.write("The list of numbers corresponding to the various taxa was written to the file %s.\n" % (OutFileNameTL))
	if DataLevels > 2:
		OutFileNamePL = OutFolder+OutFilePre+"Pop_Key.txt"
		OutFile = open(OutFileNamePL, 'w')
		OutFile.write("\n".join(PopList))
		OutFile.close()
		print("The list of numbers corresponding to the various populations was written to the file %s.\n" % (OutFileNamePL))
		sys.stderr.write("The list of numbers corresponding to the various populations was written to the file %s.\n" % (OutFileNamePL))
elif OutFormat == "GenePop":
	OutFileName = OutFolder+OutFilePre+"GenePop.txt"
	SNPListOut = ["GenePop file "+OutFilePre+" from the alignments listed in "+AlListIn+"\n"]
	Line = ", ".join(KeyList)+"\n"
	SNPListOut.append(Line)
	GPOTUDict = { }
	for OTU in OTUList:
		TempString = ""
		for SNPName in KeyList:
			TempString += " " 
			for Num in OTUSeqDict[OTU][SNPName]:
				TempString+=OTUSeqDict[OTU][SNPName][Num]
		TempString = TempString.replace('A','01')
		TempString = TempString.replace('C','02')
		TempString = TempString.replace('G','03')
		TempString = TempString.replace('T','04')
		TempString = TempString.replace('-9','00')
		GPOTUDict[OTU] = OTU+","+TempString+"\n"
	if DataLevels == 1:
		Line = "POP\n"
		SNPListOut.append(Line)
		for OTU in OTUList:
			SNPListOut.append(GPOTUDict[OTU])
	elif DataLevels == 2:
		OutListKey = [ ]
		for Taxon in TaxIndDict:
			Line = "POP\n"
			SNPListOut.append(Line)
			for OTU in TaxIndDict[Taxon]:
				SNPListOut.append(GPOTUDict[OTU])
			KeyLine = Taxon+"\t"+str(len(TaxIndDict[Taxon]))+"\n"
			OutListKey.append(KeyLine)
		OutFileTName = OutFolder+OutFilePre+"Taxon_Key.txt"
		OutFileT = open(OutFileTName, 'w')
		for KeyLine in OutListKey:
			OutFileT.write(KeyLine)
		OutFileT.close()
		print("The individuals are divided according to taxa and a key to the taxa was written to the file %s.\n" % (OutFileTName))
		sys.stderr.write("The individuals are divided according to taxa and a key to the taxa was written to the file %s.\n" % (OutFileTName))
	elif DataLevels == 3:
		OutListKey = [ ]
		for Pop in PopIndDict:
			Line = "POP\n"
			SNPListOut.append(Line)
			for OTU in PopIndDict[Pop]:
				SNPListOut.append(GPOTUDict[OTU])
			KeyLine = Pop+"\t"+str(len(PopIndDict[Pop]))+"\n"
			OutListKey.append(KeyLine)
		OutFilePName = OutFolder+OutFilePre+"Pop_Key.txt"
		OutFileP = open(OutFilePName, 'w')
		for KeyLine in OutListKey:
			OutFileP.write(KeyLine)
		OutFileP.close()
		print("The individuals are divided according to population and a key to the populations was written to the file %s.\n" % (OutFilePName))
		sys.stderr.write("The individuals are divided according to population and a key to the populations was written to the file %s.\n" % (OutFilePName))
elif OutFormat == "adegenet":
	SNPListOut = [ ]
	RScriptList = [ ]
	OutFileName = OutFolder+OutFilePre+"adegenet.str"
	OutFileNameR = OutFolder+OutFilePre+"adegenet.r"
	#making the body of the main output file
	for OTU in OTUList:
		IndName = OTU+"\t"+OTU[0:2]+"\t"+OTU.split("_")[0]+"\t"
		for Num in range(MaxSNPs):
			TempString = ""
			for SNPName in KeyList:
				TempString += OTUSeqDict[OTU][SNPName][Num]+"\t"
			#remove the last tab character
			TempString = TempString[:-1]
			TempString = TempString.replace('A','0')
			TempString = TempString.replace('C','1')
			TempString = TempString.replace('G','2')
			TempString = TempString.replace('T','3')
			Line = IndName+TempString+"\n"
			SNPListOut.append(Line)
	#making the parameter file
	Line = 'library("ape")\nlibrary("pegas")\nlibrary("seqinr")\nlibrary("ggplot2")\nlibrary("adegenet")\nlibrary("ade4")\n\nrm(list=ls())\n'
	RScriptList.append(Line)
	LineStart = 'GBSData <- read.structure("'+OutFileName+'", n.ind = '+str(len(OTUList))+', n.loc = '+str(len(OTUSeqDict[OTU].keys()))+', onerowperind = FALSE, col.lab = 1, col.pop = '
	if (DataLevels == 1) or (DataLevels == 2):
		Line = '#These data can only be analyzed by taxon:\n'
		Line += LineStart+'2, col.others = 3, row.marknames = FALSE, pop = 0)\n\n'
	elif DataLevels == 3:
		Line = '#These data can be analyzed by population or by taxon.\n#Comment out the one you do not want to use.\n'
		Line += '#for analysis by populations:\n'+LineStart+'3, col.others = 2, row.marknames = FALSE, pop = 0)\n'
		Line += '#for analysis by taxa:\n'+LineStart+'2, col.others = 3, row.marknames = FALSE, pop = 0)\n\n'
	RScriptList.append(Line)
	#The rest of this stuff is all the same, regardless of 
	#doing the original PCA
	Line = '#missing values are treated as the mean of the other values\nX <- scaleGen(GBSData, NA.method = "mean")\n#pca, will ask you to choose the PCs to use\n'
	Line += 'pca1 <- dudi.pca(X, cent=FALSE, scale=FALSE)\n#add number of PCs to this to repeat the same analysis\n#pca1 <- dudi.pca(X, cent=FALSE, scale=FALSE, scannf=FALSE, nf=xx)\n#barplot of the PCs\nbarplot(pca1$eig[1:50], main="PCA eigenvalues", col = heat.colors(50))\n\n\n'
	RScriptList.append(Line)
	#plotting the PCA
	Line = '#simple plot with all individulas labeled\ns.label(pca1$li)\ntitle("PCA of GBS Data")\nadd.scatter.eig(pca1$eig[1:20], 3, 1, 2)\n\n'
	Line += '#plot with populations labeled, PCs 1 and 2\ns.class(pca1$li, pop(GBSData))\ntitle("PCA of GBS Data\naxes 1-2")\nadd.scatter.eig(pca1$eig[1:20], 3, 1, 2)\n\n'
	Line += '#plot with populations labeled, PCs 1 and 3\ns.class(pca1$li, pop(GBSData), xax = 1, yax = 3, sub="PCA 1-3", csub=2)\ntitle("PCA of GBS Data\naxes 1-3")\nadd.scatter.eig(pca1$eig[1:20], nf=3, xax=1, yax = 3)\n\n'
	Line += '#colorful plot, PCs 1 and 2\ncol <- funky(nPop(GBSData))\ns.class(pca1$li, pop(GBSData), xax=1, yax=2, col = transp(col,.6), axesell = FALSE,\n        cstar = 0, cpoint = 3, grid = FALSE)\n\n\n'
	RScriptList.append(Line)
	#finding groups for DAPC
	if (DataLevels == 1) or (DataLevels == 2):
		Line = '###DAPC\n#choose fewer than NPops/3 PCs\ngrp <- find.clusters(GBSData, max.n.clust=nPop(GBSData))\n'
	elif DataLevels == 3:
		Line = '###DAPC\n#choose fewer than NPops/3 PCs\ngrp <- find.clusters(GBSData, max.n.clust=nPop(GBSData)*2)\n'
	Line += '#to fill out to replicate analysis\n#grp <- find.clusters(GBSData, max.n.clust=nPop(GBSData)*2, n.pca=xx, n.clust=xx)\n\n'
	Line += '#comparing the recovered groups with the populations\ntable.value(table(pop(GBSData), grp$grp), col.lab = paste ("inf", levels(grp$grp)), row.lab = levels(pop(GBSData)))\n'
	Line += '#Like with Structure, it is not clear what the best number of groups is.\n#You want to try several different numbers.\n#And you need to run the find.clusters several times (at least 10) for each number of groups, because you will likely get different group compositions in the different runs.\n\n' 
	RScriptList.append(Line)
	#figuring out the best number of PCs to retain with DAPC
	Line = '#running DAPC:\n#to find the optimal number of pcs to retain:\n#first run the analysis with a very large number of pcs\n'
	Line += 'dapc2 <- dapc(GBSData, n.da=100, n.pca=nInd(GBSData))\n#then get a rough idea of the a.score (this function evaluates multiple values and interpolates for the rest)\ntemp <- optim.a.score(dapc2)\n\n'
	Line += '#Then iteratively run dapc until you figure out the number of PCs to optimize the a.score\n#Then you can investigate the various options within that window more closely\ndapc2 <- dapc(GBSData, n.da=100, n.pca=xx)\ntemp <- optim.a.score(dapc2)\n\n'
	RScriptList.append(Line)
	#running the DAPC with the chosen number of PCs
	Line = '#now rerunning it with the optimal number of PCs\ndapc1 <- dapc(GBSData, grp$grp)\n#fill out xx to replicate original analysis\n'
	Line += '#n.pca should be the optimal number of pcs from the above analysis\n#dapc1 <- dapc(GBSData, grp$grp, n.pca=xx, n.da=xx)\n\n'
	RScriptList.append(Line)
	#plotting the results
	Line = '#complicated plot\nMyCol <- c ("darkblue", "purple", "green", "orange", "red", "blue")\n'
	Line += 'scatter(dapc1, ratio.pca = 0.3, bg = "white", pch = 20, cell = 0, cstar = 0,\n        col = MyCol, solid = 0.4, cex = 3, clab = 0, scree.da =FALSE,\n        leg = TRUE, txt.leg = paste("Cluster", 1:summary(dapc1)$n.pop))\n'
	Line += 'par(xpd = TRUE)\npoints(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4, cex=3, lwd=8, col="black")\npoints(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4, cex=3, lwd=2, col=MyCol)\n'
	Line += 'myInset <- function(){\n  temp <- dapc1$pca.eig\n  temp <- 100*cumsum(temp)/sum(temp)\n  plot(temp, col=rep(c("black", "lightgrey"), c(dapc1$n.pca,1000)), ylim=c(0,100),\n       xlab = "PCA axis", ylab="Cumulated variance (%)", cex=1, pch=20, type="h", lwd=2)\n}\n'
	Line += 'add.scatter(myInset(), posi="bottomright", inset=c(-0.03, -0.01), ratio=.28, bg=transp("white"))\n\n'
	RScriptList.append(Line)
	#plots for how DAPC populations fell out on the different PC axes
	Line = '#or looking at how the individuals fall out with a single discriminant function:\n#(This does not work if some of the groups only have one individual.)\n'
	Line += 'scatter(dapc1, 1, 1, col=MyCol, scree.da=FALSE, legend=TRUE, solid=0.4)\ntitle(paste("Discriminant Function 1: eigenvalue ", dapc1$eig[1]))\n'
	Line += 'scatter(dapc1, 2, 2, col=MyCol, scree.da=FALSE, legend=TRUE, solid=0.4)\ntitle(paste("Discriminant Function 2: eigenvalue ", dapc1$eig[2]))\n'
	Line += 'scatter(dapc1, 3, 3, col=MyCol, scree.da=FALSE, legend=TRUE, solid=0.4)\ntitle(paste("Discriminant Function 3: eigenvalue ", dapc1$eig[3]))\n\n'
	RScriptList.append(Line)
	#and how well the individuals are assigned to their previous populations
	Line = '#assignment plot to show how well individuals are assigned to their populations\nassignplot(dapc1)\n'
	Line += '#assignplot(dapc1, subset = 1:30) #to only look at a subset of the data\n\n'
	RScriptList.append(Line)
	#structure-like plots
	Line = '#or similar to structure\ncompoplot(dapc1, posi="bottomright", txt.leg=paste("Cluster", 1:summary(dapc1)$n.pop),\n          lab="", ncol=1, xlab="individuals")\n'
	Line += '#looking for most admixed individuals\nadmixedInds <- which(apply(dapc1$posterior, 1, function(e) all(e<0.9)))\nadmixedInds\n'
	Line += '#plotting them\ncompoplot(dapc1, posi="bottomright", txt.leg=paste("Cluster", 1:summary(dapc1)$n.pop),\n          ncol=1, subset=admixedInds)\n\n'
	RScriptList.append(Line)
	#comparing the DAPC groups to those found with random assignment
	Line = '#the a score compares the groups found with dapc to those found with random assignment\ngroupTesting <- a.score(dapc1)\n'
	Line += 'groupTesting$pop.score\ngroupTesting$mean\n'
	RScriptList.append(Line)
	OutFile = open(OutFileNameR, 'w')
	for Line in RScriptList:
		OutFile.write(Line)
	OutFile.close()
	print("The R script to analyze this file is %s.  It is best executed in a gui program such as RStudio.\n" % (OutFileNameR))
	sys.stderr.write("The R script to analyze this file is %s.  It is best executed in a gui program such as RStudio.\n" % (OutFileNameR))
	print("If you make changes to the R script (to add the number of pcs you actually use, etc.), then you need to save it with a different name, such as %sadegenet_date.r, to avoid writing over it when this script is rerun.\n" % (OutFilePre))
	sys.stderr.write("If you make changes to the R script (to add the number of pcs you actually use, etc.), then you need to save it with a different name, such as %sadegenet_date.r, to avoid writing over it when this script is rerun.\n" % (OutFilePre))
elif OutFormat == "fasta":
	#This format is for finding SNPs that are different between different populations.
	SNPListOut = [ ]
	OutFileName = OutFolder+OutFilePre+"SNPs.fa"
	OutFileNameList = OutFolder+OutFilePre+"SNP_list.txt"
	#making the body of the main output file
	for OTU in OTUList:
		for Num in range(MaxSNPs):
			Line = ">"+OTU+"_"+str(Num+1)+"\n"
			for SNPName in KeyList:
				SNP = OTUSeqDict[OTU][SNPName][Num]
				if SNP == "-9":
					SNP = "-"
				Line += SNP
			Line += "\n"
			SNPListOut.append(Line)
	OutFile = open(OutFileNameList, 'w')
	for KeyName in KeyList:
		OutFile.write(KeyName+"\n")
	OutFile.close()
	print("The list of SNP names was written to the file %s.\n" % (OutFileNameList))
	sys.stderr.write("The list of SNP names was written to the file %s.\n" % (OutFileNameList))


#writing the main output file (regardless of format)
OutFile = open(OutFileName, 'w')
for Line in SNPListOut:
	OutFile.write(Line)
OutFile.close()
print("The output in %s format was written to the file %s.\n" % (OutFormat, OutFileName))
print("A maximum of %s SNPs were included for each locus.\n" % (str(SNPsperLocus)))
sys.stderr.write("The output in %s format was written to the file %s.\n" % (OutFormat, OutFileName))
sys.stderr.write("A maximum of %s SNPs were included for each locus.\n" % (str(SNPsperLocus)))
if (DataLevels == 1) or (OutFormat == "fasta"):
	print("The samples were not divided into subgroups.\n")
	sys.stderr.write("The samples were not divided into subgroups.\n")
elif DataLevels == 2:
	print("The samples were divided into taxa, but not into populations.\n")
	sys.stderr.write("The samples were divided into taxa, but not into populations.\n")
elif DataLevels == 3:
	print("The samples were divided into taxa and populations.\n")
	sys.stderr.write("The samples were divided into taxa and populations.\n")
