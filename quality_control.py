#!/usr/bin/python

# script to read fastq files and do quality control: run trimmomatic 0.3 to remove Illumina adaptors and drop bad reads, fastx for trimming low quality bases, and generate fastqc output.

# COMPLETE PATHS TO INPUT FASTQ FILES AND OUTPUT DIRECTORY ARE REQUIRED IN COMMAND LINE

import sys, os, subprocess

fastq_1 = sys.argv[1]; # argv[0] is the script name. argv[1] contains first file name supplied at command line.
fastq_2 = sys.argv[2]; # contains paired read file
out_loc = sys.argv[3] + "/quality_control"; # location where output files are to be stored. A new folder called 'quality control' will be created

# change to output directory
os.makedirs(out_loc);
os.chdir(out_loc);
print ("\nWorking Directory:");
subprocess.call("pwd");
#exit(1);

#run trimmomatic from Trimmomatic-0.3 directory
trimmomatic_loc = "/mnt/home/agupta/HIV_Sequencing/analysis/Trimmomatic-0.30/trimmomatic-0.30.jar";
trimmomatic_adaptor_file = "/mnt/home/agupta/HIV_Sequencing/analysis/Trimmomatic-0.30/adapters/TrimmomaticPalindromeSimple.fa";

trim_cmd = "java -jar "+trimmomatic_loc+" PE -phred33 "+fastq_1+" "+fastq_2+" s1_palindrome_pe s1_palindrome_se s2_palindrome_pe s2_palindrome_se ILLUMINACLIP:"+trimmomatic_adaptor_file+":2:30:10 LEADING:3 TRAILING:3";

# this is correct subprocess.call, but this doesn't return output, so use Popen instead:
#subprocess.call(['java', '-jar', trimmomatic_cmd, 'PE', fastq_1, fastq_2, 's1_palindrome_pe', 's1_palindrome_se', 's2_palindrome_pe', 's2_palindrome_se', "ILLUMINACLIP:"+trimmomatic_adaptor_file+":2:30:10"]); 

trim_out= subprocess.Popen(trim_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT);

print ('\nPrinting trimmomatic command:');
print (trim_cmd, '\n');
print ('Trimmomatic Output:');
print (trim_out.communicate(), '\n');

# load fastx
subprocess.call(['bash', '-c', 'module load GNU/4.4.5']);
subprocess.call(['bash', '-c', 'module load FASTX/0.0.13']);
# load fastQC
subprocess.call(['bash', '-c', 'module load FastQC/0.10.1']);

input_files = ["s1_palindrome_pe", "s1_palindrome_se", "s2_palindrome_pe", "s2_palindrome_se"];
for file in input_files:
	outfile = file +".filt";
	print ("\nRunning fastq_quality_filter on "+file+"... Output in "+outfile);
	#qc_cmd = "fastx_trimmer -Q33 -i "+file+" | fastq_quality_filter -Q33 -q 30 -p 50 -o"+outfile;
	#subprocess.call(['fastx_trimmer', '-Q33', '-t', '10', '-i', file, '|', 'fastq_quality_filter', '-Q33', '-q', '30', '-p', '50', '-o', outfile]);
	subprocess.call(['fastq_quality_filter', '-Q33', '-q', '30', '-p', '50', '-i', file, '-o', outfile]);
	print ("\nRunning fastQC on "+outfile+" ...");
	subprocess.call(['fastqc', outfile]);
print ("\nCompleted quality control analysis. Files are in "+out_loc);

