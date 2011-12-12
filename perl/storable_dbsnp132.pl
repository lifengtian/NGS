use Storable qw/store retrieve/;
use strict;

my $fn="dbSNP132_avHet.txt";
my $out="dbSNP132_avHet.storable";

my %hash;
open IN, $fn;
while(<IN>){
	my @a=split;
	$hash{$a[0]}=$a[1];
}

store \%hash, $out;
close IN;


#$hashref = retrieve();

