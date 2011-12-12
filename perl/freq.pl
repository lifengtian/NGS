
use strict;
use warnings;
use Text::CSV;

my ( $fn ) = @ARGV;

my $csv = Text::CSV->new();

open IN, $fn or die "error open $fn";

while(<IN>){

	if ( $csv->parse($_) ) {
		my @cols=$csv->fields();
		(print, next ) if $cols[6] =~/^$/ ;
		print $_ if $cols[6] < 0.05 ;
	}
}

