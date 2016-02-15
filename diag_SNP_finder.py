#! /usr/bin/env python

#diag_SNP_finder.py version 1.1 10 Feb. 2016 Abby Moore
#This script finds SNPs that have a minimum frequency difference between two groups.
#The groups can either be two subgroups of the individuals, or they can be
#one subgroup and the rest of the individuals.
#It expects the input sequences to be in fasta format, with two lines per individual.
#version 1.1 adds the option of including a third group of individuals in the output,
#even though they will not be used to choose the SNPs.  This is the putative hybrid
#population.

import sys #We need to talk to the command line.
from collections import defaultdict #We want to be able to make dictionaries with multiple levels
from Bio import AlignIO #We want to be able to read alignments
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Align import AlignInfo

'''
diag_SNP_finder.py AlFileName SNPFileName Group1FileName Group2FileName Group3FileName FreqDiff OutFolder OutFilePre
diag_SNP_finder.py allperpop_3LLCSNPs.fa allperpop_3LLCSNP_list.txt CA_inds_LLC.txt all none 0.8 ~/Documents/seqs_formatted_outliers/diagSNPs20160206/ 3LLC_
pwd: /home/abby/Documents/seqs_formatted_outliers/diagSNPs20160206/
'''

Usage = '''
diag_SNP_finder.py version 1.0
This script reads fasta files of SNPs created by altoSNPs.py and finds SNPs
that are diagnostic for certain groups.
diag_SNP_finder.py
[alignment name]
[list of SNPs, one per line, in the same order as in the alignment, or "none", 
if you do not have this file and do not want to know the SNP names]
[file containing the list of individuals in the first group]
[file containing the list of individuals in the second group or "all", if the 
first group of individuals is to be compared to all of the remaining individuals]
[file containing the list of individuals in the third group--these are 
individuals that are not included in the comparison, but will be included in the
output, or "none", if no third group--if the second group is "all", the thrid
group must be "none"]
[difference in the frequency of a given nucleotide between the two groups for it
to be considered diagnostic, e.g., 0.8, or 1.0 if you are only interested in
diagnostic SNPs that are fixed in one of the groups]
[output folder] 
[prefix for output file]
'''

print ("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 9:
	sys.exit("ERROR! There were too few input parameters.  Eight additional arguments are required, and you supplied %d.\n  %s" % (len(sys.argv)-1, Usage))
AlFileName = sys.argv[1]
SNPFileName = sys.argv[2]
Group1FileName = sys.argv[3]
Group2FileName = sys.argv[4]
Group3FileName = sys.argv[5]
if (Group2FileName == "all") and (Group3FileName != "none"):
	sys.exit("ERROR! If the comparison group is all remaining sequences, there cannot be a third group included in the output.\n %s" % (Usage))
FreqDiff = float(sys.argv[6])
OutFolder = sys.argv[7]
if OutFolder[-1] != "/":
	OutFolder += "/"
OutFilePre = sys.argv[8]
if OutFilePre == "none":
	OutFilePre == ""

#############################################################################################################

#CaptureColumn makes a list from a specified column in a file.  This is useful
#for reusing various text files that previous programs needed.
#The column number has to follow python numbering, so 0 for the first, 1 for the
#second, etc.
#from tseq_placer_dup.py
def CaptureColumn(FileName, ColNum):
	TempList = [ ]
	InFile = open(FileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\n').strip('\r').split('\t')
		TempList.append(Line[ColNum])
	InFile.close()
	print("Column number %d was read from the file %s.\n" % (ColNum+1, FileName))
	#sys.stderr.write("Column number %d was read from the file %s.\n" % (ColNum+1, FileName))
	return TempList
	#This is SNPList, Grp1List, possibly Grp2List


#############################################################################################################

#reading the alignment
Alignment = AlignIO.read(AlFileName, "fasta")
print("The alignment of %d SNPs in %d individuals was read from file %s.\n" % (Alignment.get_alignment_length(), len(Alignment)/2, AlFileName))
sys.stderr.write("The alignment of %d SNPs in %d individuals was read from file %s.\n" % (Alignment.get_alignment_length(), len(Alignment)/2, AlFileName))
#sys.stderr.write("The alignment of %d SNPs in %d individuals was read from file %s.\n" % (Alignment.get_alignment_length(), len(Alignment)/2, AlFileName))

#figuring out which sequences go with which individuals
IndDict = { }
AllIndsList = [ ]
for Sequence in Alignment:
	SeqName = Sequence.id
	IndName = ("_").join(SeqName.split("_")[:2])
	IndDict[SeqName] = IndName
	AllIndsList.append(IndName)
AllIndsList = sorted(list(set(AllIndsList)))

#reading the names of the SNPs
if SNPFileName == "none":
	print("No file of names for the SNPs was given.\n")
	sys.stderr.write("No file of names for the SNPs was given.\n")
else:
	SNPList = CaptureColumn(SNPFileName, 0)
	print("The names of %d SNPs were read.\n" % (len(SNPList)))
	sys.stderr.write("The names of %d SNPs were read.\n" % (len(SNPList)))

#reading the names of the individuals in the first group
Grp1List = CaptureColumn(Group1FileName, 0)
print("%d individuals will be in the first group.\n" % (len(Grp1List)))
sys.stderr.write("%d individuals will be in the first group.\n" % (len(Grp1List)))

#making a list of individuals in the second group
#either from the second list file
if Group2FileName != "all":
	Grp2List = CaptureColumn(Group2FileName, 0)
	OVSet = set(Grp1List).intersection(Grp2List)
	if len(OVSet) != 0:
		sys.exit("ERROR!!!  The following individuals are present in groups one and two: %s\n" % (", ".join(sorted(list(OVSet)))))
	print("%d individuals will be in the second group.\n" % (len(Grp2List)))
	sys.stderr.write("%d individuals will be in the second group.\n" % (len(Grp2List)))
	print("And %d individuals will not be examined further.\n" % (len(AllIndsList)-len(Grp1List)-len(Grp2List)))
	sys.stderr.write("And %d individuals will not be examined further.\n" % (len(AllIndsList)-len(Grp1List)-len(Grp2List)))
#or by taking the rest of the sequences that were not part of the first group
elif Group2FileName == "all":
	Grp2List = list(set(AllIndsList)-set(Grp1List))
	print("The remaining %d individuals will be in the second group.\n" % (len(Grp2List)))
	sys.stderr.write("The remaining %d individuals will be in the second group.\n" % (len(Grp2List)))

#make the list for the third group, if applicable
if Group3FileName != "none":
	Grp3List = CaptureColumn(Group3FileName, 0)
	OVSet = list(set(Grp3List).intersection(Grp2List))
	if len(OVSet) != 0:
		sys.exit("ERROR!! The following individuals are present in groups two and three: %s\n" % (", ".join(sorted(list(OVSet)))))
	OVSet = list(set(Grp3List).intersection(Grp1List))
	if len(OVSet) != 0:
		sys.exit("ERROR!! The following individuals are present in groups one and three: %s\n" % (", ".join(sorted(list(OVSet)))))

#Make two separate alignments for the two groups
for Sequence in Alignment:
	IndName = IndDict[Sequence.id]
	if IndName in Grp1List:
		try:
			Alignment1.append(Sequence)
		except NameError:
			Alignment1 = MultipleSeqAlignment([Sequence])
	elif IndName in Grp2List:
		try:
			Alignment2.append(Sequence)
		except NameError:
			Alignment2 = MultipleSeqAlignment([Sequence])
	elif (Group3FileName != "none") and (IndName in Grp3List):
		try:
			Alignment3.append(Sequence)
		except NameError:
			Alignment3 = MultipleSeqAlignment([Sequence])
if Group3FileName ==  "none":
	print("Two new alignments have been made.  There are %d sequences and %d individuals in the first alignment.\n" % (len(Alignment1), len(Alignment1)/2))
	print("And there are %d sequences and %d individuals in the second alignment.\n" % (len(Alignment2), len(Alignment2)/2))
else:
	print("Three new alignments have been made.  There are %d sequences and %d individuals in the first alignment.\n" % (len(Alignment1), len(Alignment1)/2))
	print("There are %d sequences and %d individuals in the second alignment.\n" % (len(Alignment2), len(Alignment2)/2))
	print("And there are %d sequences and %d individuals in the third alignment.\n" % (len(Alignment3), len(Alignment3)/2))


#Go through each position in the alignment and find the nucleotide frequencies at that position
#for SeqPos in range(0,Alignment1.get_alignment_length()):
AlSum1 = AlignInfo.SummaryInfo(Alignment1)
AlSum2 = AlignInfo.SummaryInfo(Alignment2)
NucList = ['A', 'C', 'G', 'T']
NucListGap = ['A', 'C', 'G', 'T', '-']
GapIgnList = [ ]
GapInclList = [ ]
for SeqPos in range(0,Alignment1.get_alignment_length()):
	#print(SeqPos)
	#I do not know what AlSum1 and AlSum2 are doing here, but it seems to require them.
	Site1Freqs = AlSum1._get_letter_freqs(SeqPos, Alignment1, NucList, ['-'])
	Gap1Freqs = AlSum1._get_letter_freqs(SeqPos, Alignment1, NucListGap, [ ])
	#print(Site1Freqs)
	#print(Gap1Freqs)
	Site2Freqs = AlSum2._get_letter_freqs(SeqPos, Alignment2, NucList, ['-'])
	Gap2Freqs = AlSum2._get_letter_freqs(SeqPos, Alignment2, NucListGap, [ ])
	#print(Site2Freqs)
	#print(Gap2Freqs)
	for Nuc in NucList:
		if (Gap1Freqs['-'] < 0.5) and (Gap2Freqs['-'] < 0.5):
			NucDiff = Site1Freqs[Nuc]-Site2Freqs[Nuc]
			if abs(NucDiff) >= FreqDiff:
				#print("At site %d, there are significantly different nucleotide frequencies when gaps are ignored: %s in alignment 1 is %0.2f while in alignment 2 it is %0.2f.\n" % (SeqPos, Nuc, Site1Freqs[Nuc], Site2Freqs[Nuc]))
				GapIgnList.append(SeqPos)
	for Nuc in NucListGap: 
		NucDiff = Gap1Freqs[Nuc]-Gap2Freqs[Nuc]
		if abs(NucDiff) >= FreqDiff:
			#print("At site %d, there are significantly different nucleotide frequencies when gaps are included: %s in alignment 1 is %0.2f while in alignment 2 it is %0.2f.\n" % (SeqPos, Nuc, Gap1Freqs[Nuc], Gap2Freqs[Nuc]))
			GapInclList.append(SeqPos)
GapIgnList = sorted(list(set(GapIgnList)))
GapInclList = sorted(list(set(GapInclList)))
print("There was a total of %d sites with nucleotide frequencies that differed by at least %0.2f when gaps are ignored.\n" % (len(GapIgnList), FreqDiff))
print("There was a total of %d sites with nucleotide frequencies that differed by at least %0.2f when gaps are included.\n" % (len(GapInclList), FreqDiff))

#Making sub-alignments:
#First for when the gaps are ignored
for SeqPos in GapIgnList:
	try:
		#SeqPos:SeqPos+1 is obviously the same thing as SeqPos, but writing it this way gives you small alignment and writing it the other way gives you a string.
		NewGapIgnAl += Alignment1[:, SeqPos:SeqPos+1]
		GapIgnAl2 += Alignment2[:, SeqPos:SeqPos+1]
		if Group3FileName != "none": 
			GapIgnAl3 += Alignment3[:, SeqPos:SeqPos+1]
	except NameError:
		NewGapIgnAl = Alignment1[:, SeqPos:SeqPos+1]
		GapIgnAl2 = Alignment2[:, SeqPos:SeqPos+1]
		if Group3FileName != "none":
			GapIgnAl3 = Alignment3[:, SeqPos:SeqPos+1]
print("For the alignments of the sites the significant differences when gaps are ignored, the first alignment had %d sequences and was %d nucleotides long.\n" % (len(NewGapIgnAl), NewGapIgnAl.get_alignment_length()))
print("The second alignment had %d sequences and was %d nucleotides long.\n" % (len(GapIgnAl2), GapIgnAl2.get_alignment_length()))
if Group3FileName != "none":
	print("The third alignment had %d sequences and was %d nucleotides long.\n" % (len(GapIgnAl3), GapIgnAl3.get_alignment_length()))
if Group3FileName != "none":
	NewGapIgnAl.extend(GapIgnAl3)
NewGapIgnAl.extend(GapIgnAl2)
print("The combined alignment had %d sequences and was %d nucleotides long.\n" % (len(NewGapIgnAl), NewGapIgnAl.get_alignment_length()))
AlFileNameOut = OutFolder+OutFilePre+"gaps_ignored.fa"
OutFile = open(AlFileNameOut, "w")
AlignIO.write(NewGapIgnAl, OutFile, "fasta")
OutFile.close()
print("This alignment was written to the file %s.\n" % (AlFileNameOut))
sys.stderr.write("The new alignment when gaps are ignored was written to the file %s.\n" % (AlFileNameOut))
SNPFileNameOut = OutFolder+OutFilePre+"SNPnames_gaps_ignored.txt"
OutFile = open(SNPFileNameOut, 'w')
for SeqPos in GapIgnList:
	OutFile.write(SNPList[SeqPos]+"\n")
OutFile.close()
print("The list of different SNPs was written to the file %s.\n" % (SNPFileNameOut))

#Then for when the gaps are treated as characters
for SeqPos in GapInclList:
	try:
		NewGapInclAl += Alignment1[:, SeqPos:SeqPos+1]
		GapInclAl2 += Alignment2[:, SeqPos:SeqPos+1]
		if Group3FileName != "none":
			GapInclAl3 += Alignment3[:, SeqPos:SeqPos+1]
	except NameError:
		NewGapInclAl = Alignment1[:, SeqPos:SeqPos+1]
		GapInclAl2 = Alignment2[:, SeqPos:SeqPos+1]
		if Group3FileName != "none":
			GapInclAl3 = Alignment3[:, SeqPos:SeqPos+1]
print("For the alignments of the sites the significant differences when gaps are treated as characters, the first alignment had %d sequences and was %d nucleotides long.\n" % (len(NewGapInclAl), NewGapInclAl.get_alignment_length()))
print("The second alignment had %d sequences and was %d nucleotides long.\n" % (len(GapInclAl2), GapInclAl2.get_alignment_length()))
if Group3FileName != "none":
	print("The third alignment had %d sequences and was %d nucleotides long.\n" % (len(GapInclAl3), GapInclAl3.get_alignment_length()))
if Group3FileName != "none":
	NewGapInclAl.extend(GapInclAl3)
NewGapInclAl.extend(GapInclAl2)
print("The combined alignment had %d sequences and was %d nucleotides long.\n" % (len(NewGapInclAl), NewGapInclAl.get_alignment_length()))
AlFileNameOut = OutFolder+OutFilePre+"gaps_coded.fa"
OutFile = open(AlFileNameOut, "w")
AlignIO.write(NewGapInclAl, OutFile, "fasta")
OutFile.close()
print("This alignment was written to the file %s.\n" % (AlFileNameOut))
sys.stderr.write("The new alignment when gaps were treated as characters was written to the file %s.\n" % (AlFileNameOut))
SNPFileNameOut = OutFolder+OutFilePre+"SNPnames_gaps_coded.txt"
OutFile = open(SNPFileNameOut, 'w')
for SeqPos in GapInclList:
	OutFile.write(SNPList[SeqPos]+"\n")
OutFile.close()
print("The list of different SNPs was written to the file %s.\n" % (SNPFileNameOut))
