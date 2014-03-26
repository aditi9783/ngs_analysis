#!/usr/bin/python

# identify probable SNP positions based on following rule:
# 1. Most abundant base have freq < 99.5%
# 2. Remaining bases have a non-uniform base distribution => at least one non-major base is more frequent than others, and thus is a variant

# This is to be followed by part 2, where snps that have linkage information with other SNPs are kept and the rest are discarded.

import sys, os, re

mismatch_file = sys.argv[1];
maindir = sys.argv[2]; # location of main directory for a dataset. For eg: /mnt/scratch/agupta/MiSeq_Mar2013/TcellHIVanalysis
out_loc = maindir + "/variants"; # create a subdirectory 'variants' to save snps output by the program

# create output directory
os.makedirs(out_loc);

# mismatch file: containing counts of all bases at each reference position
mismatch_f = open(mismatch_file, 'r');
out_f = open(out_loc+"/snp_pos_based_on_basefreq.out", 'w');
out_snps = open(out_loc+"/snps_based_on_basefreq.out", 'w');

next(mismatch_f); #discard first line that contains a comment

for line in mismatch_f:
	content = re.split('\s+', line);
	# index guide for content:
	# 0: ref pos
	# 1: ref base
	# 2: mpileup coverage
	# 3-7: freq of A, T, C, G, N
	# 8-12: counts of A, T, C, G, N
	# 13: counts of '*'
	
	# skip position if coverage is < 10
	if ( int(content[2]) < 10 ):
		continue;
	else:
		# get max freq at this position
		freq = [];
		maxf_idx = 0;
		max_freq = 0.0;
		base_f = 0.0;
		bases = ["A", "T", "C", "G"]; # order of bases for freq and counts.
		for i in range(3,7): # gets indices 3,4,5,6 => freq of bases
			base_f = float(content[i]);
			freq.append( base_f );
			if ( base_f > max_freq ):
				max_freq = base_f;
				maxf_idx = i;
			else:
				continue;

		if (max_freq < 0.995):
			# test if remaining bases form a uniform distribution
			counts = [];
			for i in range(8,12): # get counts for bases A,T,C,G
				if (i != maxf_idx+5): # index of base with max freq + 5 gives the counts for same base
					counts.append( float(content[i]) );
			sum_c = sum(counts);
			# expected count for each non-major base
			exp_c = sum_c/3; # since 4th base is the abundant base, thus divide by 3
			
			# compute chi-square statistic
			chi_sq_cal = 0.0;
			for c in counts:
				chi_sq_cal += (c - exp_c)**2 / exp_c;
			chi_sq_crit = 5.991; # chi sq value for df=2 and alpha=0.05

#			print "line is ", line;
#			print "sum: ", sum_c, " exp_c ", exp_c, " and chisq_cal ", chi_sq_cal;

			if (chi_sq_crit < chi_sq_cal): # position is a candidate snp, print it
				out_f.write(line); # print the entire line to the snp file
				
				# print every non-major base as a snp that has a base freq > 0.33% (= 0.0033). 
				# Rationale being if 1% is total error rate, per base error is 0.33%.
				# base with higher freq is a probable snp
				# ASSUMPTION: major base is the reference base
				for idx in range(0,4):
					if (idx == maxf_idx-3): # maxf_idx is index of major base, with indexing starting from 0 instead of 3.
						continue;
					if (freq[idx] > 0.0033):
						out_snps.write(content[0]+" "+bases[idx]+" "+content[idx+8]+" "+content[2]+"\n"); # print snp pos, base, # times this base was seen at this position, and coverage at this position
		else:
			continue; # not a variant position as majority base is 99.5% abundant

