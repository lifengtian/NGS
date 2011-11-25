#!/usr/bin/perl 

### select unique blat mapper
### column 10 in psl file is the seqID
### Lifeng Tian
### March 2010
###

my %h=();

<>;<>;<>;<>;<>;
while(<>){
	my @lines=split/\t/;
	my $key = $lines[9];

	$h{$key}++ ;	

	print $_ if $h{$key} == 1;

}
