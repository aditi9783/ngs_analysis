import sys, re, string

file = open('Tcell_HIV_CTTGTA_L001_R1_R2_bwa_mem_alnpe.bam.sorted.bam.pileup','r');
#cov_file = open('Tcell_HIV_CTTGTA_L001_R1_R2_bwa_mem_alnpe.bam.sorted.bam.pileup.cov','r');

# read reference position, bases at that site, and base qualities from pileup file.
# for every position, list the number of A/T/C/G in a list format (only consider bases that have base quality >= 30).

#baseC = []; # list of base counts for each position in the reference
QSCORE_THRES = 30;

for line in file:
	content = re.split('\s+',line);  # split each line in pileup output at spaces
	# content: [seq_id, pos, N, cov, base_string, qual_string]

#	print content;
	# if cov (# reads) < 50, don't consider that position
	if (int(content[3]) < 50):
#		baseC.append( [0,0,0,0] ); # append empty base counts
		continue;
	else:
		# if cov > 50, get base counts from bases with base qual >=30.
		seq = re.sub('\^.','',content[4]); # character following "^" has ASCII char with q-score. This ASCII char can be [atcg], thus removing it before the next step, where non [atcg] characters are removed.
		seq = re.sub('[^atcgATCG]+', '', seq); # only keep bases from the sequence string 
		seq_u = seq.upper(); 	# convert to uppercase
		seq_num = seq_u.translate(string.maketrans("ATCG", "0123")); # transliterate bases to numbers that serve as index in count list	
		bases = [int(b) for b in seq_num]; # get the numerical representation of bases in a list format
		counts = [0,0,0,0]; # counts for A,T,C,G

		# convert base qualities into q-score, and count base only if q-score >= 30
		q_scores = [ord(c)-33 for c in content[5]]; # ord gets ASCII value for character, and q-score is ASCII value minus 33.
		if (len(bases) != len(q_scores)):
			print seq;
			print len(bases), bases, '\n', len(q_scores), q_scores;
			print "number of bases from reads do not match number of q-scores at position ", content[1], "\nExiting\n\n!";
			sys.exit(1);
		for i in range(len(bases)):
			if (q_scores[i] >= QSCORE_THRES):
				counts[ bases[i] ] += 1;
			else:
				continue;
#		baseC.append( counts );	
		sumC = float(sum(counts));
		countsfrac = [c/sumC for c in counts];
#		print content[1], seq_u, bases, counts, q_scores, sumC, countsfrac;
		print content[1], content[3], sumC, countsfrac[0], countsfrac[1], countsfrac[2], countsfrac[3];
#	else:
#		print "Exiting!\n";
#		sys.exit(1);
