#!/usr/bin/perl 
# remove XA and MD from bowtie generated sam file for cufflinks compatability
# Lifeng Tian
# Jan 2010

while(<>){
	chomp;
	~s/XA:\w*:\w*\s*//;
	~s/MD:\w*:\w*\s*//;
	~s/AS:\w*:\w*\s*//;
	chop;
	print $_,"\n";
	#print $_,"\tXS:A:-\n" ;
}
