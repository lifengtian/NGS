
use strict;
use warnings;
use Text::CSV;

my ( $fn, @col ) = @ARGV;

my $csv = Text::CSV->new();

open IN, $fn or die "error open $fn";

while(<IN>){

	if ( $csv->parse($_) ) {
		my @cols=$csv->fields();
		print join("\t",@cols[@col]),"\n";
	}
}

