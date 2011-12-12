
use strict;
use warnings;
use Text::CSV;
use Storable;


#my ( $fn ) = @ARGV;

my $usage="perl create_avHet.pl\n";

my $csv = Text::CSV->new();

my $hash=retrieve("dbSNP132_avHet.storable");

foreach my $fn ( `ls *.csv` ) {
open IN, $fn or die "error open $fn";
open OUT, ">avHet/$fn";


print "processing $fn\n";
<IN>;#skip header
while(<IN>){

	chomp;
	if ( $csv->parse($_) ) {
		my @col=$csv->fields();
	
		if ( $col[9] =~/^rs/ ) {
			print OUT  $_,",",$hash->{$col[9]},"\n"; 
		} else {
			print OUT $_,",\"NA\"\n";
		} 
	}
}

print "End\n";
close IN;
close OUT;
}

__END__



