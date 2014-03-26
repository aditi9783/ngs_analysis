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


base_f.write("#<reference position> <reference base> <mpileup coverage> <freq:A> <freq:T> <freq:C> <freq:G> <freq:N> <counts:A> <counts:T> <counts:C> <counts:G> <counts:N> <counts:*>\n");
for line in file:
	content = re.split('\s+',line);  # split each line in pileup output at spaces
	# content: [seq_id, pos, N, cov, base_string, qual_string]

#	print content;
	# if cov (# reads) < 5, don't consider that position
	if (int(content[3]) < 5):
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
			counts = {'A':0, 'T':0, 'C':0, 'G':0, 'N':0, '*':0}; # counts for A,T,C,G
			indel_flag = 0;
			indel_len = 0;
			indel = [];
			refbase = content[2].upper();
			for b in seq_list:
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
						
			sumC = float(sum( counts.values() ));
	#		print content[1], refbase, content[3]; 
	#		print "Printing counts:";
			freq = {'A':0.0, 'T':0.0, 'C':0.0, 'G':0.0, 'N':0.0, '*':0.0}; # freq for A,T,C,G
			for b, c in counts.iteritems():
	#			print b, " : ", c;
				freq[b] = c/(sumC - counts['*']);
				
			base_f.write(str(content[1]) +" "+ refbase +" "+ str(content[3]) +" "+ str(freq['A']) +" "+ str(freq['T']) +" "+ str(freq['C'])+" "+ str(freq['G']) +" "+ str(freq['N']) +" "+ str(counts['A']) +" "+ str(counts['T']) +" "+ str(counts['C'])+" "+ str(counts['G']) +" "+ str(counts['N']) +" "+ str(counts['*'])+"\n" );
			if (sumC != float(content[3])):
				print "position ", content[1], " Error!! Coverage doesn't match sum of bases.\n";
				print "Exiting!";
				sys.exit(1);
#	else:
#		print "Exiting!\n";
#		sys.exit(1);
