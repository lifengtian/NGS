
use strict;
use warnings;
use Text::CSV;

my ( $fn ) = @ARGV;

my $csv = Text::CSV->new();

open IN, $fn or die "error open $fn";
<IN>;
while(<IN>){

	if ( $csv->parse($_) ) {
		my @cols=$csv->fields();
		print join("\t",@cols[15..17,0..14,18..25]),"\n";
	}
}

