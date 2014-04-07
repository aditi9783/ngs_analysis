# read mpileup output and write out the positions that have a variation (i.e. more than one base at a position).
# also calculate base frequencies at these positions.
# insertions and deletions from the reference are written to a different file, along with the position that encountered that indel.
# quality control used a quality score threshold of 30, so base qualities are not looked at here. All bases are considered.
# script mpileup_snps.py has code to consider quality scores as well.

import sys, re, string

mpileup_filename = sys.argv[1]; # mpileup file
file = open(mpileup_filename,'r');

base_f = open(mpileup_filename+'.mismatches', 'w');
indel_f = open(mpileup_filename+'.indels','w');

QSCORE_THRES = 30; # qscore threshold for removing low quality bases
COV_THRES = 2; # threshold for coverage at each pos: is cov is less than thres, don't consider that position

base_f.write("#<reference position> <reference base> <mpileup coverage> <effective coverage= sum of A, T, C, G that have q scores >= threshold> <freq:A> <freq:T> <freq:C> <freq:G> <counts:A> <counts:T> <counts:C> <counts:G>\n");
for line in file:
	content = re.split('\s+',line);  # split each line in pileup output at spaces
	# content: [seq_id, pos, N, cov, base_string, qual_string]

#	print content;
	# if cov (# reads) < 5, don't consider that position
	if (int(content[3]) < COV_THRES):
		continue;
	else:
		# if the pileup sting contains bases (indicating mismatch/indel) then process that line in mpileup file
		var_present = re.search('[ATCGatcg]+', content[4]);
	#	print "pos:", content[1], " cov:", content[3], " var_present:", var_present, " pileup:", content[4];

		if var_present:
			seq = re.sub('\^.','',content[4]); # character following "^" has ASCII char with read mapping quality. This ASCII char can be [atcg], thus removing it before the next step, where non [atcg] characters are removed.
			# '^' marks the begining of a read, and '$' marks the end, so remove '$' too
			seq = re.sub('\$', '', seq); # removing '$' symbols
	#		seq = re.sub('\*', '', seq); # removing '*' symbol that indicates deletion in the read (or a gap?)
			# deleting '*' results in sum of base counts to be different from coverage reported. Thus, remove '*' count before calculating base fractions, and print '*' count.
		        seq_upper = seq.upper(); # convert to upper case	
			# split seq string into an array
			seq_list = list(seq_upper);

			if ( int(content[3]) != len(seq_list)):
				print "Coverage and seq list are different lengths: cov is ", content[3], " and seq len is ", len(seq_list);

			# convert base qualities into q-score, and count base only if q-score >= 30
			q_scores = [ord(c)-33 for c in content[5]]; # ord gets ASCII value for character, and q-score is ASCII value minus 33.

			#print line;
			#print q_scores;
	
			# check if sequence string and quality string are the same length. Exit if they aren't
			#if (len(seq_list) != len(q_scores)):
			#	print seq;
			#	print len(seq_list), seq_list, '\n', len(q_scores), q_scores;
			#	print "number of bases from reads do not match number of q-scores at position ", content[1], "\nExiting\n\n!";
			#	sys.exit(1);

			low_qual_bases = 0;
			for qs in q_scores:
				if (qs < QSCORE_THRES):
					low_qual_bases += 1;

			# since we dont' want to consider low quality bases, if the effective coverage is less than coverage threshold, skip the position
			effective_cov = int(content[3]) - low_qual_bases;
			if (effective_cov < COV_THRES):
				continue;
			else:
	
				counts = {'A':0, 'T':0, 'C':0, 'G':0, 'N':0, '*':0}; # counts for A,T,C,G
				indel_flag = 0;
				indel_len = 0;
				indel = [];
				refbase = content[2].upper();
				for i in range(len(seq_list)):
					if (q_scores[i] >= QSCORE_THRES): # only considering bases that have q-scores greather than threshold
						b = seq_list[i]; # ith base in the sequence string
						if (indel_flag == 0):
							if (b == '.' or b == ','):
								counts[ refbase ] += 1;
							elif (b == '+' or b == '-'):
								indel_flag = 1; 
								indel.append(b);
							else:
								counts[ b ] += 1;
						elif (indel_flag == 1):
							if (indel_len == 0):
								if (b.isalpha() or b == '*'): # b is a alphabet or '*'
									indel_len = int( "".join( indel[1:] ) ); # indel list: indel[0] is +/-, followed by indel len before first base
								indel.append(b); # save the current indel base, which is also the first indel base
								if (indel_len == 1):# only one indel base, i.e. current base, so write that out
									indel_f.write(str(content[1]) + " : " + "".join( indel ) +"\n" );
									# reset flags/lists
									indel = [];					
									indel_flag = 0;
									indel_len = 0;
							
							else: # indel len is > 1, so keep saving correct indel bases and adjusting remaining indel_len.
								indel.append(b);
								indel_len -= 1;
								if (indel_len == 1): # indel len == 1 means current base is last indel base, so write it out. 
									indel_f.write(str(content[1]) + " : " + "".join( indel ) +"\n" );
									# reset flags/lists
									indel = [];					
									indel_flag = 0;
									indel_len = 0;
			
				# only interested in counts for A, T, C, G (not interested in N or *)
				sumC = float( counts['A'] + counts['T'] + counts['C'] + counts['G'] ); # this is effective coverage
	#		print content[1], refbase, content[3]; 
	#		print "Printing counts:";
				freq = {'A':0.0, 'T':0.0, 'C':0.0, 'G':0.0}; # freq for A,T,C,G
				novariation_flag = 0; # after getting the high qscore bases only, it may turn out that there is no mutation/variation at the pos
				for b, c in freq.iteritems():
	#			print b, " : ", c;
					freq[b] = counts[b]/sumC;
					if (freq[b] == 1.0):
						novariation_flag = 1; # if any one base has freq=1.0, there is no variation (thus not an interesting position)
				
				# print only if there is a variation at the position
				if (novariation_flag == 0):
					base_f.write(str(content[1]) +" "+ refbase +" "+ str(content[3]) +" "+ str(sumC) +" "+ str(freq['A']) +" "+ str(freq['T']) +" "+ str(freq['C'])+" "+ str(freq['G']) +" "+ str(counts['A']) +" "+ str(counts['T']) +" "+ str(counts['C'])+" "+ str(counts['G']) +"\n" );
	#		if (sumC != float(content[3])): # sum of bases may not match coverage if some bases are removed due to low qscore
	#			print "position ", content[1], " Error!! Coverage doesn't match sum of bases.\n";
	#			print "Exiting!";
	#			sys.exit(1);
#	else:
#		print "Exiting!\n";
#		sys.exit(1);
