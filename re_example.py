#!/usr/bin/python

import re

str = "+3ATCTATGT";
ilen = re.search('\d+',str);
istr = re.search('[ATCGN]+', str);

print "ilen ", ilen;
print ilen.group();
indellen = int( ilen.group() );

print "istr ", istr;
print istr.group();
indel = istr.group();

if (indellen != len(indel) ):
	i = list(indel);
	for a in range(indellen):
		print indel[a];	
	print "Remaining bases:";
	print indel[indellen:];
	for a in range(indellen,len(indel)):
		print indel[a].isalpha();
		print indel[a];
num = "2";
print num.isalpha();
base = 'T';
print base, base.isalpha();
