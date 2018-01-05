This is mainly a set of wrapper scripts to get the RTD pipeline by Peterson et al. to run with a large number of individuals on a desktop.
Most of these scripts are probably not useful anymore, because improved analysis methods (and more computational power) are generally available; they are simply included for completeness.
The most useful scripts are probably those starting at alstats.py, which take fasta-formatted alignments of each locus, choose the loci to use based on various criteria, and create output files formatted for different programs. 

Filter and trim reads:
(fastq_comb_4.py)
fastq_comb_4.py [name(s) of files with forward reads] >> [name of log file]
**This step takes a long time and should definitely be run overnight or, preferably, on a cluster.**
This looks through all sequences, looks for restriction sites (GATC), trims the ends off the sequences (both the ends before the restriction sites and any adapter sequences from the end of the sequence), and rejects the sequences that lack restriction sites or have too many bases that are of low quality.  Then it reverse-complements the reverse sequences, determines whether forward and reverse sequences overlap, and concatenates forward and reverse sequences (while removing overlapping bases, if any exist).  It tries to do this on its own, but it uses Muscle to try to combine the sequences that it can't combine (for example because of errors in sequencing).  Thus, you may need to change the name it uses for Muscle, depending on what your computer calls it. This is easy to find in the file by searching for “muscle”.  Make sure you change the line that is not commented out.

Sort reads by barcode:
(fastqcomb_pbc.py)
fastqcomb_pbc.py [barcode file name] [sequence file name(s)] >> [name of log file]
This reads a text file of the individuals present and their barcodes, labels all of the sequences by barcode, and then writes files for each individual containing its barcodes.  To avoid overwhelming the computer’s memory, it writes separate files for each individual for each input file.  It combines these files at the end so there is one file per individual.  It also writes a shell script that deletes the temporary files (although it does not run this script automatically).  Naming the files only works if the sequence file names are the full path, instead of just the file name within the directory.

Make a files of the unique sequences:
(fastqpbc_uniqued_rtd.py)
fastqpbc_uniqued_rtd.py [desired length for uniqued sequences, 30 bp works well] [name of outfile, without “.uniqued.gz”, which the program adds automatically] [name of input files—the individuals that you want to be part of the master alignment] >> [name of log file]
Makes a file with sequences, quality scores, and information about which individuals have each sequence.  I think it is best to use just a subset of individuals for this process (say 25 individuals spread across the taxa of interest).  This allows it to run much faster and makes it easier to add more individuals as the study progresses.  This script will make two files of uniqued sequences, one with sequences of the length you told it to make when you called the script and one with full-length sequences.  You will need the shorter sequences for the first alignment step and the longer sequences for the second alignment step.

Alignment the first 30 bases:
(rtd_run.py)
This uses a pipeline written by Peterson for RAD sequences (Peterson, B. K., J. N. Weber, E. H. Kay, H. S. Fisher, and H. E. Hoekstra.  2012.  Double digest RADseq: an inexpensive method for de novo SNP discovery and genotyping in model and non-model species.  PLoS ONE 7(5): e37135.).  It is written in python, but calls a bunch of other programs.  It also does not work perfectly and makes very large temporary files.
rtd_run.py [outfile path] [infile name(s)]
cd [outfilepath]
new_subj_rtd.py
rtd_run.py [outfile path] [infile name(s)]

(You have to call rtd_run twice, because it somehow can't make one of the files that it needs, so it stops partway through with an error, then new_subj_rtd.py makes that file, and then it can go on.  It should be able to make the alignment on the second run, but, even then, it generally isn’t able to finish without running into some errors.  However, it should be able to get far enough to make the .cluni and the .clstats files that you need for the next step.)

Figure out which of these alignments you want to use:

You get a whole bunch of files with very long names, but you only need two of them:

The (long name).clstats file:
0	25946	154	0.416008621028

0: alignment number
25946: number of sequences in this alignment
154: number of individuals that have a sequence in this alignment
0.416...: „dirt“

And the (somewhat different long name) .cluni file:
0	1000076.uniqued_allinds30.uniqued	GATCCGTACCACAGTTAGGCTGATTCACGT	1	?@<DDBBDF?:22<?F9<<C?CH9F+A8)?	50_3	1	.	.

group number, sequence name, sequence, 1= forward read, quality score, first individual that has that sequence, # of times that sequence is present in that individual (maybe?), something about reverse reads (which are not present here)

Figure out which sequences are interesting to analyze further:
(bam30_uniq_split_5.py)
bam30_uniq_split_5.py [lower limit of group size] [upper limit of group size] [lower limit of group size for smaller clusters] [.clstats file name from rtd_run.py] [.cluni file name from rtd_run.py] [file containing uniqued (full) sequences—the second file made by fastqpbc_uniqued_rtd.py] [outfile path] [outfile prefix--without the path] >> [name of log file]

Its default is the maximum sequence number of 1000 per alignment and the maximum dirt of 0.05.  You need to input the upper and lower limit of the group size (meaning the number of individuals per group).  If you analyze sequences from 25 individuals, then a lower limit of 20 and an upper limit of 25 would be good numbers to start out with.  The “group size for smaller clusters” is the group size for the second set of rtd_run.py runs, where you are analyzing each individual cluster separately.  I do not know what the best number is to use here, but I have been trying 5.  The .cluni and .clstats files have very long names.  Since you will only have one of each, you can just write [filepath]/*.cluni and [filepath]/*.clstats.
This looks through the .clstats file to figure out which of the alignments meet the search criteria and  then looks through the .cluni file for all of the sequences in those alignments.  It makes a list with the sequences and the alignment number.  (When you look at the list, it will be in a random order, because it comes from a python dictionary.)
This file also writes a shell script that runs rtd_run.py on the individual sets of sequences and also calls clusters2template.py, template_concat.py, blast, and blast_parse_template_2.py.  These are described below, even though you will not have to call them directly:

(clusters2template.py)
clusters2template.py [lower limit of group size] [group name] [.clstats file] [.cluni file] [out file folder]
This script takes the output from the second rtd_run.py run (the one that analyzes one group of sequences, not all of the sequences) and decides which of those sequences should be used for the final analysis (for the reference-guided assembly you will run on all of your sequences, not just those of the 25 individuals used to choose which regions to analyze).

(template_concat.py)
template_concat.py [output file prefix] [list of template files]
This script takes the results of the various clusters2template.py runs and concatenates them together in a single, large file.  Then it writes a script to delete the intermediate files.

(blast)
makeblastdb -in [input file name (.fa)] -out [output database file name (no suffix)] -dbtype nucl
blastn -db [database name (no suffix)] -query [input file name (.fa)] -out [output file name] -outfmt '6 std qlen slen'
In this step, BLAST aligns the potential reference sequences against each other and against the primer sequences so that sequences that are reverse-complements of each other, other duplicate sequences, and any remaining stray primer sequences can be taken out.  So in this case, the input file for making the database and query file for the actual blast search are the same.  They are all fasta files, except the database, which is a group of files in various formats (but only needs to be called by the overall database name, without the extensions).

(blast_parse_template_4.py)
blast_parse_template_4.py [output filename for unique sequences] [output file prefix, generally just the folder, for sets of matching sequences] [input filename—BLAST against sequences] [input filename—BLAST against primers] [input filename—actual sequences] >> [name of log file]
The output files will be fasta files and one text file and the input files are the two files output by the two BLAST searches and the original sequence file (since the BLAST output files do not include the sequences).  This script looks through the BLAST output to find the duplicate sequences and sequences with primers and writes a new fasta file for further analysis by BLAST that does not include these sequences.  It looks through the duplicate sequences that each only have one match to see if the matches are simply reverse-complements of each other and, if so, writes those to the file of unique sequences.  It writes the remaining sets of matching sequences to individual files and writes a shell script for Muscle to align them and a text file listing the alignments and sequence names for annotation.

**At this point, it is important to look through the alignments and mark the sequences in the text file, either by *(tab) for wanted sequences and x(tab) for unwanted sequences, or by deleting the unwanted sequences.**

(template_making.py)
template_making.py [file name for the marked up list of grouped sequences, output by blast_parse_template_4.py and modified by hand] [fasta file containing all of the sequences—both with and without matches] [fasta file of template sequences without matches output by blast_parse_template_4.py] >> [name of log file]
This script simply adds the sequences selected from the groups of matching sequences to the file of non-matching sequences to make a fasta file from which the new database can be made.

(blastdb again to make the cleaned database file, suggested suffix: _cldb)

Making the fasta file for reference-guided assembly using BLAST:
(fastqpbc_uniqued_fasta.py)
fastqpbc_uniqued_fasta.py  [output file name—no suffix] [input file names] >> [name of log file]
This script makes a fasta file including all of the unique sequences.  The sequences are named as follows: [sequence number (sequential number from 1 to the number of sequences)]_[name of first individual that has that sequence]_[number of times that sequence is in individual #1]_[2nd individual's name]_[number of times that sequence is in the 2nd individual]_[etc.]

Running BLAST: 

blastn -db [database name, [your_name]_cldb] -query [file output from fastqpbc_uniqued_fasta.py [your_name]_uniq.fa] -out [output file name] -outfmt '6 std qlen slen'
This BLASTs a set of sequences against the database.  It is on the slow side, and may be faster if you only do one population/group of 5-10 individuals at once.  (I think it has RAM problems if you analyze too many sequences at once.)  If you do them in groups, it takes about an hour and a half to BLAST a whole set of 25 individuals on my laptop.

Analyzing the BLAST data: (Basically a test run)
(blast_parse_output.py)
blast_parse_output.py [prefix for output files for statistics] [folder in which to output the sequence alignments] [file of BLAST output] [fasta file of uniqued sequences] >> [name of log file]
This script looks through the output from BLASTing the uniqued sequences against the cleaned template sequences.  It keeps the sequences that have one match to a template sequence with more than 90% sequence similarity over more than 90% of their length (although this can be changed).  It writes these sequences to separate files for each template sequence.  (The sequences are just appended to the files, so existing sequences in those files are not overwritten.)  It also writes various files of statistics: a list of all sequences and the sequences they match (more than 90,90): _SeqStats.txt, a list of the number of sequences that had a given number of matches: _NumMatchesDict.txt, a list of each pair of matching sequences and how often that pair occurred: _SeqGroups.txt, and a list of sequences with multiple matches and the template sequences to which they matched: _AmbigSeqs.txt.

(seq_group_concat.py)
seq_group_concat.py [output filename for the template sequences and their numbers of matches in each population] [fasta file including all template sequences] [output filename for the fasta file of filtered template sequences] [list of _SeqGroups.txt files from blast_parse_output.py] >> [log filename]
This script finds the template sequences that had the most non-unique matches (currently 100 or more) over all of the populations.  Then it makes a new fasta file of template sequences that does not include these sequences.  The entire process must now be repeated with this fasta file.

Makeblastdb (same as always)

Analyzing the BLAST data (final run!):

blastn for all of the files of uniqued sequences

blast_parse_ind.py  [folder in which to output the sequence alignments] [file of BLAST output] [fasta file of uniqued sequences] >> [name of log file]
This script looks through the output from BLASTing the uniqued sequences from one individual against the final template sequences.  It keeps the sequences that have one match to a template sequence with more than 90% sequence similarity over more than 90% of their length (although this can be changed).  It writes these sequences to separate files for each template sequence. It also writes a script to align all of the sequence files for each locus using muscle.
***This script is run once for each individual, and the folder for the sequence alignments needs to be a folder for that individual.  Thus, if the individual is CA217_1 in group Moore1, then the folder could be Moore1CA217_1

Muscle to align all of the sequence files: AlignmentScript.sh written by blast_parse_ind.py

(seqparse.py)
seqparse.py [minimum number of sequencess that must have an alternative nucleotide for
the position to be considered polymorphic—I try 3] [minimum percent of sequences that must have
a given allele for that allele to be considered to be present—I try 0.05] [ploidy] [individual name] 
[folder with the final sequence files] [list of sequence files to parse—this can be something like *.fa, but it is not a file with a list, just the list itself]
Example: ./seqparse.py 3 .1 2 CA217_1 ./finalseqs1 ./Moore1CA217_1/a*.fa >> Moore1inds.log
It uses both a minimum number of sequences and a minimum percent of sequences so you have some coverage when the total number of sequences is small (like 10) and also something realistic when the total number of sequences is very large (like 10,000).  At this point, ploidy just gives you error messages when you have more alleles at a given locus than the ploidy level of the plant.  Individual name is just the prefix for the sequences.
This script looks through the aligned fasta files produced by blast_parse_ind.py and adds the (potential) true alleles (after filtering for sequence errors) to the final sequence files for each locus. 

Muscle to align the various files.

Take the alignment files and figure out which loci are worth using:
(alstats.py)
alstats.py [folder containing alignment files] [file containing list of alignment files] [minimum number of taxa that must have a sequence for it to be considered further] [minimum proportion of the sequences that an individual must have for it to be analyzed] [final proportion of taxa in which the sequences must occur] [final proportion of populations in which the sequences must occur] [final proportion of individuals in which the sequences must occur] [distance above the median variability (in terms of number of interquartile ranges) at which sequences will be considered too variable and be removed] [category to be used in making the subset of the data; options are taxa, pops, inds, and none] [file containing the list of taxa/pops/inds to include, or none, if none] [folder in which output files should be written]
This script looks through the finished alignments produced by seqparse.py and looks at which individuals have which sequences.  Then it takes out the loci that are too variable, the uncommon loci, and the failed individuals and calculates everything again.  It makes lists of the numbers of sequences from each allele that each individual, population, and taxon has and makes lists of the sequences in each alignment from each individual, the best individual from each population, and the best population of each taxon.  These will be used to make the alignments for the final analysis.
Output files (no prefix, so they have to go in a new folder so previous files aren't overwritten):
common_pseq_counts.txt: table of populations by loci, showing the number of individuals from each population that have at least one sequence from each locus (only for the loci that are present in more than a certain percent of individuals)
common_tseq_counts.txt: table of taxa by loci, showing the number of populations from taxon that have at least one sequence from each locus (sequences can be from different individuals for the different loci) (only for the loci that are present in more than a certain percent of individuals)
inds_pops.txt: table of the individual in each population that has the most sequences and the number of sequences it has
pops_taxa.txt: table of the population of each taxon that has the most sequences and the number of sequences it has
pseq_counts.txt:  table of populations by loci, showing the number of individuals from each population that have at least one sequence from each locus (for all loci)
seqs_inds.txt: table showing all sequences for each locus in each individual
seqs_pops_best.txt: table showing the sequences for each locus for the best individual of each population (the same individual for each locus)
seqs_pops.txt: table showing all sequences for each locus for each population
seqs_taxa_best.txt: table showing the sequences for each locus for the best population of each taxon (individuals can differ at different loci)
seqs_taxa.txt: table showing all sequences for each locus for each taxon
tseq_counts.txt: table of taxa by loci, showing the number of populations from taxon that have at least one sequence from each locus (sequences can be from different individuals for the different loci) (for all loci)

Make new sets of alignments for the subset of loci/individuals you want to analyze:
(newal.py)
newal.py [folder containing alignment files] [prefix for old alignment files] [suffix for old alignment files] [output file from alstats.py that should be used to make the new alignments] [number of sequences per set (taxon, population, or individual, depending on the alstats.py output file selected) to be used in the alignment--all if all are to be used] [folder into which new alignment files should be put] [prefix for new alignments--none if no prefix]
newal.py makes new alignments following one of the lists produced by alstats.py.  It also allows you to specify the number of sequences per OTU that each alignment should have, which is important for something like structure, where all individuals will need to have the same ploidy. 

Taking these alignment files and formatting them for further analyses:
Analyses of SNP data:
(altoSNPs.py)
altoSNPs [folder containing the alignments] [prefix for the alignments] [alstats.py output file that was used to make the alignments] [output folder]  [prefix for output file] [output format] [# SNPs per locus: either "all" or an integer] [number of levels in output data: 1 if there is one individual per taxon, 2 if there is one individual per population, 3 if there are multiple individuals per population]
The available output formats are structure, fastStructure, GenePop, and adegenet (a package in R).
The commands or a file to run structure, fastStructure, and adegenet are included in the output or written to the screen.  GenePop is a format that is used by various different programs, and no specific commands are given for it.
altoSNPs.py takes the new alignments and converts them into the format of choice.

Analyses of sequence data:
(altoconbs.py)
altoconbs.py [folder containing the alignments] [prefix for the alignments] [alstats.py file used to make the alignments]  [output folder] [prefix for output file] [output format: fasta, phylip, or nexus] [# of bootstrap datasets]
altoconbs.py makes a concatenated sequence file from the separate files for each locus made by newal.py.  It can also make bootstrapped datasets where the loci (not the individual alignment positions) are resampled.  These files can be output in fasta, phylip, or nexus format.
When phylip output is selected, it also makes a script for raxml.  When nexus output is selected, it makes a paup script for running an SVDQuartets analysis.  This script has 500 bootstrap replicates, at each of the taxon, population, and individual levels, so it will take a while to run.