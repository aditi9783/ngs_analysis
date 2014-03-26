import sys, re

filename = sys.argv[1];
file = open(filename,'r');
cov_file = open(filename+'.coverage','w');

for line in file:
	content = re.split('\s+',line);
	#print content[1], "\t", content[3];
	cov_file.write(content[1]+ "\t" +content[3]+ "\n");
