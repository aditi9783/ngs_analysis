#!/usr/bin/python

# read list of probable snps as output by step 1, and do the following:
# : output a file that contains the simple integer representation of read ids for each true snp

# the list of probable snps is in a file called 'snps_based_on_basefreq.out' in 'variants' subdirectory in the main analysis directory for a dataset. 
# File 'snp_pos_based_on_basefreq.out' contains detailed information about base counts at this snp position.

# upload modules GNU/4.4.5, Python/2.7.2, and PySAM/0.6 to use pysam
import sys, re, pysam
from operator import itemgetter

main_dir = sys.argv[1]; # main dir for a dataset. eg: /mnt/scratch/agupta/MiSeq_Mar2013/TcellHIVanalysis
var_dir = main_dir+"/variants"; # location of var directory for a dataset. For eg: /mnt/scratch/agupta/MiSeq_Mar2013/TcellHIVanalysis/variants
bamfile = main_dir+"/mapping/merged.sorted.bam"; # reads mapped to reference, in bam format

# read input files
var_f = open(var_dir+"/snps_based_on_basefreq.out", 'r'); # input snp file
alnfile = pysam.Samfile(bamfile, "rb"); # read alignment file using pysam package

print "header ", alnfile.header;
print "ref ", alnfile.references;

# open output streams for writing
#snps_out = open(var_dir+"/snps_based_on_basefreq_and_linkmat.out", 'w'); # final snp file
hapread_out = open(var_dir+"/hapread.out", 'w'); # file containing read-id, snp pos, and snp base
readidmap_out = open(var_dir+"/readidmap.out", 'w'); # file mapping original read-id to its integer representation

snplist = [];
readid_map = {};
hapread = [];
counter = 0;  # to create simple int representation for each read id

for line in var_f:
	content = re.split('\s+', line);
	# index guide for content:
	# 0: snp pos
	# 1: snp base
	# 2: # times this snp base is observed at this position
	# 3: mpileup coverage at this snp pos
	snplist.append( [ int(content[0]), content[1] ] ); # save snp pos and base
#	snpindex.append( content[0]+content[1] );
#print "snplist made ... ", snplist;

#maxl = 6;
# POSSIBLE RUNTIME IMPROVEMENT: currently pysam is called for every snp, but some snps share a position. Improvement could be calling pysam once for each position.
for i in range(len(snplist)):
#	if (maxl == 0):
#		break;
#	maxl -= 1;
	readidlist = []; # list of all readids that map to this snp
	snppos = snplist[i][0];
	snpbase = snplist[i][1];
	print "Reading snp: ", snppos, " ", snpbase;
	# get reads that map to this snp pos
	# pysam position = pileup position-1. So give snppos-1, snppos to get reads.
	for pileupcolumn in alnfile.pileup( 'K03455|HIVHXB2CG', snppos-1, snppos ):
		if (pileupcolumn.pos == snppos-1): # snp position reached
			for pileupread in pileupcolumn.pileups:
				readbase = pileupread.alignment.seq[pileupread.qpos];
				if (readbase == snpbase):
					readid = pileupread.alignment.qname;
					readidlist.append( readid );
				else:
					continue;
			# exit loop since no need to look further, since snp position and base has been found
			break;
		else:
			continue;
	# remove duplicates in read id list. A snp cannot appear multiple times in same read id. Some error with pysam output or my mistake?
#	print snppos, snpbase, "read id list: ", readidlist;
	unique_rids = list(set(readidlist));
#	print "unique read ids: ", unique_rids;
	for i in range(len(unique_rids)):
#		print "read id read: ", unique_rids[i];
		if (unique_rids[i] in readid_map):
#			print "read seen before at counter ", readid_map[ unique_rids[i] ];
			hapread.append( [readid_map[ unique_rids[i] ], str(snppos), snpbase] );
		else:
			readid_map[ unique_rids[i] ] = counter;
			counter +=1;
#			print "New read. counter now is ", counter;
			hapread.append( [readid_map[ unique_rids[i] ], str(snppos), snpbase] );

# sort the hapread list by readid_map int representation
sorted_hapread = sorted(hapread, key=itemgetter(0));

# print readid map
sorted_readidmap = sorted(readid_map.items(), key=itemgetter(1)); # sort by value, i.e., the integer rep of the readid
for pair in sorted_readidmap:
	readidmap_out.write(pair[0]+"\t"+str(pair[1])+"\n"); # original read id and its int representation

#print this sorted list
for hr in sorted_hapread:
	hapread_out.write(str(hr[0])+"\t"+hr[1]+"\t"+hr[2]+"\n");	

# close file handles
alnfile.close();	
var_f.close();
hapread_out.close();
readidmap_out.close();						
		
