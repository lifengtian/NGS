use Modern::Perl;
use List::MoreUtils qw/zip/;

=head1 note

perl combine_coverage.pl [coverage_file1] [coverage_file2] ...


0       SAMPLE
1       LIBRARY
2       READ_GROUP
3       BAIT_SET
4       GENOME_SIZE
5       BAIT_TERRITORY
6       TARGET_TERRITORY
7       BAIT_DESIGN_EFFICIENCY
8       TOTAL_READS
9       PF_READS
10      PF_UNIQUE_READS
11      PCT_PF_READS
12      PCT_PF_UQ_READS
13      PF_UQ_READS_ALIGNED
14      PCT_PF_UQ_READS_ALIGNED
15      PF_UQ_BASES_ALIGNED
16      ON_BAIT_BASES
17      NEAR_BAIT_BASES
18      OFF_BAIT_BASES
19      ON_TARGET_BASES
20      PCT_SELECTED_BASES
21      PCT_OFF_BAIT
22      ON_BAIT_VS_SELECTED
23      MEAN_BAIT_COVERAGE
24      MEAN_TARGET_COVERAGE
25      PCT_USABLE_BASES_ON_BAIT
26      PCT_USABLE_BASES_ON_TARGET
27      FOLD_ENRICHMENT
28      ZERO_CVG_TARGETS_PCT
29      FOLD_80_BASE_PENALTY
30      PCT_TARGET_BASES_2X
31      PCT_TARGET_BASES_10X
32      PCT_TARGET_BASES_20X
33      PCT_TARGET_BASES_30X

=cut

my @keys_to_report = qw/
SAMPLE
GENOME_SIZE
TARGET_TERRITORY
TOTAL_READS
PF_UQ_BASES_ALIGNED
ON_TARGET_BASES
PCT_ON_TARGET
MEAN_TARGET_COVERAGE
ZERO_CVG_TARGETS_PCT
PCT_TARGET_BASES_20X
/;


my $hash = {};

foreach my $fn ( @ARGV ){
	&get_coverage_from_file($fn, $hash);	

}

foreach (@keys_to_report ) {
	print $_,"\t",$hash->{$_},"\n";
}



sub get_coverage_from_file {
	my ($fn, $h) = @_;

	open my $fh, $fn or die "Error open $fn!";

	while(<$fh>){
		if ( /^#/ || /^$/ ) {
			next;
		}

	my $header = $_;
	my $l = <$fh>;
	my @headers = split/\t/, $header;
	my @ls = split/\t/, $l;
	$ls[0] = ($fn=~/(.*?)\./)?$1:$fn;
	if ( $#headers != $#ls ) {
		die " in file $fn, header cols($#headers) don't equal data cols($#ls). quit!";
	}
	
	my %th = zip ( @headers, @ls);
	$th{'ZERO_CVG_TARGETS_PCT'} = sprintf("%.0f",$th{'ZERO_CVG_TARGETS_PCT'}*100);
	$th{'PCT_TARGET_BASES_20X'} = sprintf("%.0f",$th{'PCT_TARGET_BASES_20X'}*100);
	$th{'MEAN_TARGET_COVERAGE'} = sprintf("%.0f",$th{'MEAN_TARGET_COVERAGE'});
	foreach (@headers) {
		$h->{$_} .= "\t".$th{$_};
	}

	#handle the PCT_ON_TARGET
	my $pct_on_target=sprintf("%.0f",$th{'ON_TARGET_BASES'}/$th{'PF_UQ_BASES_ALIGNED'}*100);
	$h->{'PCT_ON_TARGET'} .="\t".$pct_on_target;	
	last;
	}
}
