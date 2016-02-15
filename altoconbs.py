#! /usr/bin/env python

#altoconbs.py version 1.0 29 March 2015 Abby Moore
#This script reads alignments in fasta format and concatenates them.  Depending on
#how inspired I get, it may also output the data in various formats.  But I bet AlignIO can
#do that well....  It also makes bootstrapped alignments with different genes
#present or absent.
#At this point, if an OTU has multiple sequences for a given locus, it just writes them
#over each other and so be default takes the last one.

'''
example command:
altoconbs.py ~/Documents/finalseqs2_3/inds3/ 3 ~/Documents/finalseqs1output_3/seqs_inds.txt ~/Documents/seqs_formatted_outliers/SVDQuartets20151222/ 3_inds_SVDQ_ nexus 0
altoconbs.py InFolder InFilePre AlListIn OutFolder OutFilePre OutFormat [fasta, phylip, nexus] BSNum
'''

import sys #We need to talk to the command line.
from collections import defaultdict #We want to be able to make dictionaries with multiple levels
from Bio import AlignIO #We want to be able to read alignments
import numpy #We want to be able to generate new lists by sampling with replacement from old lists

Usage = '''
altoconbs.py version 1.0
This script reads alignments created by newal.py, concatenates them, and makes
bootstrapped datasets.  It can output the alignments in various formats
The available formats are: fasta, phylip, and nexus.
When calling the program, the output format must be written exactly as it is
in the previous sentence (including capitalization).
altoconbs.py
[folder containing the alignments]
[prefix for the alignments]
[alstats.py file used to make the alignments]
[output folder] 
[prefix for output file]
[output format]
[# of bootstrap datasets]
'''

FormatList = ['fasta', 'phylip', 'nexus']

print ("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 8:
	sys.exit("ERROR! There were too %d additional arguments, but seven are requried.\n %s" % (len(sys.argv)-1, Usage))
InFolder = sys.argv[1]
InFilePre = sys.argv[2]
AlListIn = sys.argv[3]
OutFolder = sys.argv[4]
OutFilePre = sys.argv[5]
OutFormat = sys.argv[6]
if (OutFormat in FormatList) == False:
	sys.exit("ERROR!  The output format can only be one of the following: %s, and you wrote %s.\n" % (" ".join(FormatList), OutFormat))
BSNum = int(sys.argv[7])

#making sure that the folder names end with a slash
if InFolder[-1] != "/":
	InFolder += "/"
print ("Sequence files will be read from the folder %s.\n" % (InFolder))
sys.stderr.write("Sequence files will be read from the folder %s.\n" % (InFolder))
if OutFolder[-1] != "/":
	OutFolder += "/"

#setting file extension:
if OutFormat == "fasta":
	OutFilePost = ".fa"
elif OutFormat == "phylip":
	OutFilePost = ".phy"
elif OutFormat == "nexus":
	OutFilePost = ".nex"

#Making the main dictionaries and lists that will be filled out:
SeqDict = defaultdict(dict) #The dictionary of the sequences by locus and individual
#It will be SeqDict[Locus][OTU] = Seq
AlDict = defaultdict(dict) #The dictionary of the actual alignments
#It will be AlDict[AlName][OTU] = List of Seqs Making up the ConSeq where AlName can be Orig or BS##
LocusList = [ ] #List of the loci
LLDict = { } #dictionary of how long the loci are (locus length dict)
ALenDict = { } #dictionary of how long the alignments are (alignment length dict)
#LLDict[Locus] = AlLength
NumSeqs = 0 #Counter for the number of sequences.
BSRange = range(1, BSNum+1) #The list of numbers of the bootstrap replicates
InFileDict = { } #Dictionary of sequence files to be read and their corresponding locus names.
OTUDict = { } #Dictionary of the sequence names and their corresponding OTUs
OTUSeqDict = defaultdict(dict) #The dictionary with lists of sequences for each OTU

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
		#adding the locus name to the list of loci
		LocusList.append(Line[0])
		#adding the sequence names (OTUSeq) to the dictionary that tells which sequence name goes to which OTU
		for OTUNum in range(1,NumOTUs+1):
			SeqListTemp = Line[OTUNum].split(",")
			for OTUSeq in SeqListTemp:
				if OTUSeq != '':
					OTUDict[OTUSeq] = KeyLine[OTUNum]
InFile.close()
NumLoci = len(LocusList)
LocusRange = range(0,NumLoci)
print("The list of %d alignments was read from the file %s.\n" % (NumLoci, AlListIn))
sys.stderr.write("The list of %d alignments was read from the file %s.\n" % (NumLoci, AlListIn))
print("The alignments contained sequences from %d OTUs.\n" % (OTUNum))
sys.stderr.write("The alignments contained sequences from %d OTUs.\n" % (OTUNum))


#Then, we need to read each alignment and add it to the dictionary.
for Locus in LocusList:
	InFileName = InFolder+InFilePre+Locus+".fa"
	Alignment = AlignIO.read(InFileName, "fasta")
	LLDict[Locus] = Alignment.get_alignment_length()
	#Reading the alignment into the SeqDict
	for Record in Alignment:
		OTU = OTUDict[Record.id]
		SeqDict[Locus][OTU] = str(Record.seq)

print("Sequences were read from the alignment files, with names of the form %s.\n" % (InFileName))
sys.stderr.write("Sequences were read from the alignment files, with names of the form %s.\n" % (InFileName))

#Now to make the concatenated sequences
AlLength = 0
AlDict['Orig'] = defaultdict(list)
for Locus in LocusList:
	BlankSeq = LLDict[Locus]*'N'
	AlLength += LLDict[Locus]
	for OTU in OTUList:
		try:
			AlDict['Orig'][OTU].append(SeqDict[Locus][OTU])
		except KeyError:
			SeqDict[Locus][OTU] = BlankSeq
			AlDict['Orig'][OTU].append(SeqDict[Locus][OTU])
ALenDict['Orig'] = AlLength

#Now to make the replicate alignments for bootstrapping:
RepsDone = 0
for BSRep in BSRange:
	LocusListTemp = numpy.random.choice(LocusList, size=NumLoci, replace=True)
	AlName = "BS"+str(BSRep)
	AlDict[AlName] = defaultdict(list)
	AlLength = 0
	for Locus in LocusListTemp:
		AlLength += LLDict[Locus]
		for OTU in OTUList:
			AlDict[AlName][OTU].append(SeqDict[Locus][OTU])
	ALenDict[AlName] = AlLength
	RepsDone += 1

print("%d replicate resampled datasets were created.\n" % (RepsDone))
sys.stderr.write("%d replicate resampled datasets were created.\n" % (RepsDone))

NumAlFiles = 0
for AlName in AlDict.keys():
	OutList = [ ]
	if OutFormat == 'fasta':
		for OTU in AlDict[AlName]:
			Line = ">"+OTU+"\n"+"".join(AlDict[AlName][OTU])+"\n"
			OutList.append(Line)
	elif OutFormat == 'phylip':
		Line = "\t"+str(NumOTUs)+"\t"+str(ALenDict[AlName])+"\n"
		OutList.append(Line)
		for OTU in AlDict[AlName]:
			Line = OTU+"\t"+"".join(AlDict[AlName][OTU])+"\n"
			OutList.append(Line)
	elif OutFormat == 'nexus':
		#We need dictionaries of the species and populations so we can make sets for paup.
		SpDict = defaultdict(list)
		PopDict = defaultdict(list)
		Line = "#NEXUS\nBEGIN DATA;\n  DIMENSIONS NTAX="+str(NumOTUs)+" NCHAR="+str(ALenDict[AlName])+";\n  FORMAT DATATYPE=DNA\n"
		Line += "  GAP=-\n  ;\nMATRIX\n"
		OutList.append(Line)
		OTUNum = 1
		for OTU in AlDict[AlName]:
			#first, add the sequence to the alignment file
			Line = "["+str(OTUNum)+"] "+OTU+"\t"+"".join(AlDict[AlName][OTU])+"\n"
			OutList.append(Line)
			#Then increase the counter by one, so the sequences can be properly numbered
			OTUNum += 1
			#Then add the individual to the two dictionaries.
			SpName = OTU[:2]
			PopName = OTU.split("_")[0]
			SpDict[SpName].append(OTU)
			PopDict[PopName].append(OTU)
		#formatting the species dict for output:
		SpDictString = ""
		for SpName in SpDict:
			SpDictString += SpName+": "+" ".join(SpDict[SpName])+", "
		SpDictString=SpDictString[:-2]
		#formatting the PopDict for output:
		PopDictString = ""
		for PopName in PopDict:
			PopDictString += PopName+": "+" ".join(PopDict[PopName])+", "
		PopDictString=PopDictString[:-2]
		#adding them to the nexus file:
		Line = ";\nEND;\nBEGIN SETS;\nTAXPARTITION species = "+SpDictString+";\nTAXPARTITION pops = "+PopDictString+";\nEND;\n"
		OutList.append(Line)
	OutFileName = OutFolder+OutFilePre+AlName+OutFilePost
	OutFile = open(OutFileName, 'w')
	for Line in OutList:
		OutFile.write(Line)
	OutFile.close()
	NumAlFiles += 1

print("%d alignments were written, which have names such as %s and are in %s format.\n" % (NumAlFiles, OutFileName, OutFormat))
sys.stderr.write("%d alignments were written, which have names such as %s and are in %s format.\n" % (NumAlFiles, OutFileName, OutFormat))

#scripts to analyze the data:
if OutFormat == 'phylip':
	OutList = ['#! /bin/bash\n\n']
	Line = "rm "+OutFolder+"RAxML*"+OutFilePre+"*\n"
	Line += "raxmlHPC -s "+OutFolder+OutFilePre+"Orig.phy -n "+OutFilePre+"Orig"+" -m GTRCAT -p 1234 -f d -w "+OutFolder+"\n"
	OutList.append(Line)
	BSTreeList = [ ]
	for BSRep in BSRange:
		Line = "raxmlHPC -s "+OutFolder+OutFilePre+"BS"+str(BSRep)+".phy -n "+OutFilePre+"BS"+str(BSRep)+" -m GTRCAT -p 1234 -f d -w "+OutFolder+"\n"
		BSTreeList.append(OutFolder+"RAxML_bestTree."+OutFilePre+"BS"+str(BSRep))
		OutList.append(Line)
	Line = "cat "+" ".join(BSTreeList)+" > "+OutFolder+OutFilePre+"bsreps.tre\n"
	Line += "raxmlHPC -t "+OutFolder+"RAxML_bestTree."+OutFilePre+"Orig -z "+OutFolder+OutFilePre+"bsreps.tre -n "+OutFilePre+"final -m GTRCAT -s "+OutFolder+OutFilePre+"Orig.phy -f b -w "+OutFolder+"\n"
	#Line += "raxmlHPC -t "+OutFolder+"RAxML_bestTree."+OutFilePre+"Orig -z "+OutFolder+OutFilePre+"bsreps.tre -n "+OutFilePre+"final -f b -w "+OutFolder+"\n"
	OutList.append(Line)
	OutFileName = OutFolder+OutFilePre+"RAxML_script.sh"
	OutFile = open(OutFileName, 'w')
	for Line in OutList:
		OutFile.write(Line)
	OutFile.close()
	print("The script for analyzing the alignments using RAxML was written to %s.\n" % (OutFileName))
	sys.stderr.write("The script for analyzing the alignments using RAxML was written to %s.\n" % (OutFileName))
elif OutFormat == 'nexus':
	OutList = ['BEGIN PAUP;\n']
	Line = "CD "+OutFolder+";\n"
	Line += "Execute "+OutFilePre+"Orig.nex;\n"
	OutList.append(Line)
	if len(SpDict) > 3:
		Line = "SVDQuartets partition=species speciesTree=yes bootstrap = yes nreps = 500;\n"
		Line += "Savetree file="+OutFilePre+"spp_bs500.tre;\n"
		OutList.append(Line)
	if len(PopDict) > 3:
		Line = "SVDQuartets partition=pops speciesTree=yes bootstrap = yes nreps = 500;\n"
		Line += "Savetree file="+OutFilePre+"pops_bs500.tre;\n"
		OutList.append(Line)
	Line = "SVDQuartets speciesTree=no bootstrap=yes nreps=500;\n"
	Line += "Savetree file="+OutFilePre+"inds_bs500.tre;\n"
	Line += "End;\n"
	OutList.append(Line)
	OutFileName = OutFolder+OutFilePre+"paup_SVDQ_script.nex"
	OutFile = open(OutFileName, 'w')
	for Line in OutList:
		OutFile.write(Line)
	OutFile.close()
	print("The paup script for analyzing the alignment using SVDQuartets was written to %s.\n" % (OutFileName))
	sys.stderr.write("The paup script for analyzing the alignment using SVDQuartets was written to %s.\n" % (OutFileName))
	print("In order to run the analysis, open paup and type 'Execute %s;'\n" % (OutFileName))
	sys.stderr.write("In order to run the analysis, open paup and type 'Execute %s;'\n" % (OutFileName))
