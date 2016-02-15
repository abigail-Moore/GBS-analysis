#! /usr/bin/env python

InFileName = "rtd.mIDs"
OutFileName = "rtd_0001of0001_subj.fa"

SubjectSeqs = [ ]
SplitLine = [ ]
AllUniqued = 0
UniqSubjects = 0
WrittenSeqs = 0

InFile = open(InFileName, 'rU')
for Line in InFile:
	Line = Line.strip('\n').strip('\r').split()
	AllUniqued += 1
	if len(Line) > 2:
		NewSeq = Line[0]
		SubjectSeqs.append(NewSeq)
		UniqSubjects += 1
InFile.close()
print("A total of %d unique sequences were found.  %d of these were found in more than one individual.\n" % (AllUniqued, UniqSubjects))
if UniqSubjects == 0: #If there are no sequences that appear more than once, then just randomly choose the first 10 for the file.
	InFile = open(InFileName, 'rU')
	LineNum = 0
	for Line in InFile:
		Line = Line.strip('\n').strip('\r').split()
		if LineNum < 10:
			NewSeq = Line[0]
			SubjectSeqs.append(NewSeq)
			LineNum += 1
	InFile.close()


OutFile = open(OutFileName, 'w')
for Line in SubjectSeqs:
	OutFile.write(">")
	OutFile.write(Line)
	OutFile.write("\n")
	SplitLine = Line.split(".")
	OutFile.write(SplitLine[2])
	OutFile.write("\n")
	WrittenSeqs += 1
OutFile.close()

print("%d sequences were written to the file %s.\n" % (WrittenSeqs, OutFileName))
