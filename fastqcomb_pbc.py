#! /usr/bin/env python

#This program is supposed to take the output from fastq_comb_3.py, which are the
#concatenated forward and reverse reads, with failed sequences removed, and sort them
#according to barcode.
#It assumes that there is only one individual per barcode (it does not sort according to
#lane first).
#First, it reads the barcodes from an input file (not compressed) and creates a dictionary
#for that lane.
#Then it opens the files of sequences, which can be compressed or not.
#For each sequence file, it does the following:
#Labels individuals according to barcode, trims off the barcode from the front of
#the sequence, and trims those bases from Phred score
#Makes a defaultdict with one (top-level) entry per individual
#Writes separate files for each individual with only its sequences
#(named after the original sequence file).
#These are fastq files but the first line of each file is the individual name.
#Then it opens the next sequence file and repeats the process.

#This generates many files, so at the end of this process, the program writes one file per
#individual with all of its sequences in that file.  The first line of this file is the
#individual's name/number.
#Then the program writes a shell script that can be executed to delete the intermediate files.
#The script should be opened and read before it is executed, to make sure that no other
#files are being deleted.


import re #We will need the regular expressions module
import sys #We want to be able to output error messages to the screen.
import gzip #We want to be able to open zipped files.
from collections import defaultdict #We need to make dictionaries with multiple levels.

#First, get the files names from the command line or print some instructions if not enough names are found..

Usage = '''
fastqcomb_pbc.py - version 1.0
Sorts cleaned, concatenated sequence data (output from fastq_comb_3.py) according to barcode.
Usage:
	fastqcomb_pbc.py (barcode file name) (sequence file names)
The file names can be either the full path, e.g., /home/moorea/... 
or the file name only if you are in that directory.

Sequence file names can either be a list of all sequence file names
or something of the form "combined*.fastq".
And the sequence files (but not the barcode file) can also be zipped (with a ".gz" ending).
'''

if len(sys.argv) <3:
	print Usage
	#DictName = "/home/abby/Documents/sandbox/Barcodes1.txt"
	#This is the dictionary (for testing purposes).
	#InFileList = ["/home/abby/Documents/sandbox/comb_Moore_1_L007_R1_head.fastq.gz"]
	#This is the infile (for testing purposes) 
else:
	DictName = sys.argv[1]
	InFileList = sys.argv[2:]

print("fastqcomb_pbc.py")
	
#%%%%%%%%%%%THE DICTIONARY%%%%%%%%%%%%%%%%%%

#First, make the dictionary.  It is a tab-delimited file with one line per barcode of the format:
#CCACAA	39_15, where CCACAA is the barcode (without restriction site) and 39_15 is the individual number
#Individual names may consist of letters, numbers, and _ only.

Dict = { }
NumBarcodes = 0

InFileDict = open(DictName, 'rU')
for Line in InFileDict:
	Line = Line.strip('\n').strip('\r').split("\t")
	Dict[Line[0]] = Line[1]
	NumBarcodes += 1
InFileDict.close()

sys.stderr.write("%d barcodes were read from the file %s.\n" % (NumBarcodes, DictName))
print("%d barcodes were read from the file %s.\n" % (NumBarcodes, DictName))
	
#%%%%%%%%%%%%PREPARATION FOR SORTING INDIVIDUALS%%%%%%%%%%%%%%%%%%%%%

TotalSeqs = 0
BadBarcodesSum = 0
GoodBarcodesSum = 0

#***This is for BamHI.  If a different restriction enzyme is used, the GATCC needs to be changed accordingly.***
MyRe = r"(\w{4,8})(GATCC\w+)" #regular expression that finds the barcode
MySub = r"\1" #regular expression for barcode
MySub2 = r"\2" #regular expression for the sequence without the barcode
MyRe3 = r"\d+%(\w+)" #regular expression that finds the individual's name
MySub3 = r"\1" #regular expression for individual's name
InFileRe1 = r"(.+/)(\w+).fastq.gz"
InFileRe2 = r"(.+/)(\w+).fastq"
InFileRe3 = r"(\w+).fastq.gz"
InFileRe4 = r"(\w+).fastq"

for InFileName in InFileList:
	#This is for if the files have the full path (from my computer):
	'''if InFileName.endswith('.gz'):
		InFile = gzip.open(InFileName) #Opens in the file
		# and prepares to name the output files.
		InFilePre = re.sub(InFileRe1, MySub, InFileName)
		InFilePost = re.sub(InFileRe1, MySub2, InFileName)
	else:
		InFile = open(InFileName,'rU')
		InFilePre = re.sub(InFileRe2, MySub, InFileName)
		InFilePost = re.sub(InFileRe2, MySub2, InFileName)'''
	#This is for the files that only have the file names (for Oscar):
	if InFileName.endswith('.gz'):
		InFile = gzip.open(InFileName) #Opens in the file
		# and prepares to name the output files.
		InFilePre = "" 
		InFilePost = re.sub(InFileRe3, MySub, InFileName)
	else:
		InFile = open(InFileName,'rU')
		InFilePre = ""
		InFilePost = re.sub(InFileRe4, MySub, InFileName)
	BadBarcodes = 0
	GoodBarcodes = 0
	TotalInds = 0
	LineNum = 0
	SeqNum = 0
	SeqDict = defaultdict(dict)

#%%%%%%%%%%%%MAKING DEFAULTDICT WITH SEQUENCES SORTED ACCORDING TO BARCODE%%%%%%%%%%%%%%%%
	for Line in InFile:
		Line = Line.strip('\n').strip('\r')
		SeqLine = (LineNum + 4) % 4
		if SeqLine == 0: #Code is copied directly (first line of fastq formatted file)
			SeqName = Line
			LineNum += 1
		elif SeqLine == 1: #Barcode is read and sequence is copied
			Barcode = re.sub(MyRe, MySub, Line)
			NewLine = re.sub(MyRe, MySub2, Line)
			LenLine = len(NewLine)
			TempSeq = NewLine+"\n"
			LineNum += 1
			SeqNum += 1
			try: #If a valid barcode is found, the individual name
			#is written onto third line.
				IndName = Dict[Barcode]
				TempSeq = TempSeq + "+" + IndName + "\n"
				GoodBarcodes += 1
			except KeyError: #If no valid barcode is found, the individual
			#name is changed to "error".
				IndName = 'error'
				BadBarcodes += 1
		elif SeqLine == 2:#The third line of the sequence (just '+' in these sequences) is not read.
			LineNum += 1
		elif SeqLine == 3: #The Phred score is copied without the bases for the barcode
			Line = Line[-LenLine:]
			TempSeq = TempSeq + Line
			LineNum += 1
			try:
				SeqDict[IndName][SeqName] = TempSeq
			except KeyError: 
				SeqDict[IndName] = defaultdict(dict)
				SeqDict[IndName][SeqName] = TempSeq
	InFile.close()
	GoodBarcodesSum += GoodBarcodes
	BadBarcodesSum += BadBarcodes
	sys.stderr.write("\nIn total, %d sequences were read from the file %s.\n" % (SeqNum, InFileName))
	sys.stderr.write("%d of them had good barcodes, while %d of them had bad barcodes.\n" % (GoodBarcodes, BadBarcodes))
	print("\nIn total, %d sequences were read from the file %s.\n" % (SeqNum, InFileName))
	print("%d of them had good barcodes, while %d of them had bad barcodes.\n" % (GoodBarcodes, BadBarcodes))

#%%%%%%%%%%%%%%%%%%%%%%%WRITING SEPARATE FILES FOR EACH INDIVIDUAL%%%%%%%%%%%%%%%%%%

#So the program does not run out of memory, it writes one file per individual per original
#fastq file, instead of one file per individual per lane.  These are then combined later.

	for IndName in Dict.values():
		NumIndSeqs = 0
		#The Outfile is named after the individual and is a zipped file.
		OutFileName = InFilePre + IndName + "_" + InFilePost + ".fastq.gz"
		OutFile = gzip.open(OutFileName, 'w')
		for SeqName in SeqDict[IndName]: #writing the sequences to the file.
			OutFile.write(SeqName)
			OutFile.write('\n')
			OutFile.write(SeqDict[IndName][SeqName])
			OutFile.write("\n")
			NumIndSeqs += 1
		OutFile.close()
		TotalSeqs += NumIndSeqs
		TotalInds += 1
		sys.stderr.write("%d sequences were found for individual %s and written to the file %s.\n" % (NumIndSeqs, IndName, OutFileName))
		print("%d sequences were found for individual %s and written to the file %s.\n" % (NumIndSeqs, IndName, OutFileName))

	sys.stderr.write("%d sequences from %d individuals from the file %s were written to separate files.\n" % (TotalSeqs, TotalInds, InFileName))
	print("%d sequences from %d individuals from the file %s were written to separate files.\n" % (TotalSeqs, TotalInds, InFileName))

#%%%%%%%%%%%%%%%%%%%%%%%PERFORMING SOME FINAL CHECKS%%%%%%%%%%%%%%%%%%%%

sys.stderr.write("In %d infiles, %d sequences with good barcodes and %d sequences with bad barcodes were found.\n" % (len(InFileList), GoodBarcodesSum, BadBarcodesSum))
#print("In %d infiles, %d sequences with good barcodes and %d sequences with bad barcodes were found.\n" % (len(InFileList), GoodBarcodesSum, BadBarcodesSum))

if TotalSeqs != GoodBarcodesSum:
	sys.stderr.write("ERORR: %d sequences with good barcodes were found, but only %d of them were written to a file.\n" % (GoodBarcodesSum, TotalSeqs))
	sys.stderr.write("The rest have gotten lost somewhere.")
	print("ERORR: %d sequences with good barcodes were found, but only %d of them were written to a file.\n" % (GoodBarcodesSum, TotalSeqs))
	print("The rest have gotten lost somewhere.")

#%%%%%%%%%%%%%%%%%%%%%%%%COMBINING THE FILES FOR EACH INDIVIDUAL%%%%%%%%%%%%%
IndFilesWritten = 0

RemoveScriptName = InFilePre + "Script_for_Removing_Files.sh"
RemoveScript = open(RemoveScriptName, 'w')
RemoveScript.write("#! /bin/bash\n\n")
TotalLines = 0
IndFilesWritten = 0

for IndName in Dict.values():
	IndFileNameOut = InFilePre + IndName + "_all.fastq.gz"
	IndFileOut = gzip.open(IndFileNameOut, 'w')
	IndFileOut.write(IndName)
	NumLines = 0
	#The first line of the outfile is the individual name
	IndFileOut.write('\n')
	for InFileName in InFileList:
		#This is for if the files have the full path (from my computer):
		'''if InFileName.endswith('.gz'):
			InFilePre = re.sub(InFileRe1, MySub, InFileName)
			InFilePost = re.sub(InFileRe1, MySub2, InFileName)
		else:
			InFilePre = re.sub(InFileRe2, MySub, InFileName)
			InFilePost = re.sub(InFileRe2, MySub2, InFileName)'''
		#This is for the files that only have the file names (for Oscar):
		if InFileName.endswith('.gz'):
			InFilePre = "" 
			InFilePost = re.sub(InFileRe3, MySub, InFileName)
		else:
			InFilePre = ""
			InFilePost = re.sub(InFileRe4, MySub, InFileName) 
		IndFileNameIn = InFilePre + IndName + "_" + InFilePost + ".fastq.gz"
		IndFileIn = gzip.open(IndFileNameIn, 'r')
		for Line in IndFileIn:
			Line = Line.strip('\n').strip('\r')
			IndFileOut.write(Line)
			IndFileOut.write("\n")
			NumLines += 1
		IndFileIn.close()
		RemoveScript.write("rm " + IndFileNameIn + "\n")
	IndFileOut.close()
	IndFilesWritten += 1
	TotalLines += NumLines
	sys.stderr.write("%d sequences were found for individual %s.\n" % (NumLines/4, IndName))
	print("%d sequences were found for individual %s.\n" % (NumLines/4, IndName))
RemoveScript.close()

sys.stderr.write("In total, %d sequences were found.\n" % (TotalLines/4))
print("In total, %d sequences were found.\n" % (TotalLines/4))

sys.stderr.write("Files for %d individuals were written.\n" % (IndFilesWritten))
print("Files for %d individuals were written.\n" % (IndFilesWritten))

sys.stderr.write("To get rid of the intermediate files, please execute the file %s.\n" % (RemoveScriptName))
sys.stderr.write("In order to do so, you will need to make it executable by typing chmod u+x %s \n" % (RemoveScriptName))
sys.stderr.write("Then simply type %s to execute it.\n" % (RemoveScriptName))
sys.stderr.write("Please open and read the list of files to be deleted to make sure that files you need are not on that list.\nTHE AUTHOR CAN ACCEPT NO RESPONSIBILITY IF THIS DELETES THE ENTIRE CONTENTS OF YOUR HARD DRIVE!!!!!!!!!\n")

print("To get rid of the intermediate files, please execute the file %s.\n" % (RemoveScriptName))
print("In order to do so, you will need to make it executable by typing chmod u+x %s \n" % (RemoveScriptName))
print("Then simply type %s to execute it.\n" % (RemoveScriptName))
print("Please open and read the list of files to be deleted to make sure that files you need are not on that list.\nTHE AUTHOR CAN ACCEPT NO RESPONSIBILITY IF THIS DELETES THE ENTIRE CONTENTS OF YOUR HARD DRIVE!!!!!!!!!\n")

sys.stderr.write("\n\n\n")
