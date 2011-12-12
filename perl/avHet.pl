
use strict;
use warnings;
use Text::CSV;

my ( $fn ) = @ARGV;

my $csv = Text::CSV->new();

open IN, $fn or die "error open $fn";

my %h;
<IN>; #skip header

while(<IN>){

	if ( $csv->parse($_) ) {
		my @col=$csv->fields();
		$h{$col[15].':'.$col[16]}++;
		#print $col[15].':'.$col[16],"\n"; 
	}
}

close IN;

my $dbsnp="/work/1a/cagapps/gatk/annovar/humandb/hg19_snp132.txt";
open SNP, $dbsnp or die "Error open $dbsnp";

while(<SNP>){
	my ($rsid,$avHET)=(split)[4,13];
	print join("\t",$rsid,$avHET),"\n";
}

close SNP;

	

