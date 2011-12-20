
=head1 note
## net.sf.picard.metrics.StringHeader
# net.sf.picard.analysis.directed.CalculateHsMetrics BAIT_INTERVALS=../../AMD/results/baits.picard TARGET_INTERVALS=../../AMD/results/baits.picard INPUT=DS103_CAGATC_L001.sorted.bam OUTPUT=DS103_CAGATC_L001.sorted.bam.coverage VALIDATION_STRINGENCY=SILENT    METRIC_ACCUMULATION_LEVEL=[ALL_READS] TMP_DIR=/tmp/lifeng2 VERBOSITY=INFO QUIET=false COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false
## net.sf.picard.metrics.StringHeader
# Started on: Sat Oct 29 10:47:21 EDT 2011

## METRICS CLASS	net.sf.picard.analysis.directed.HsMetrics
SAMPLE	LIBRARY	READ_GROUP	BAIT_SET	GENOME_SIZE	BAIT_TERRITORY	TARGET_TERRITORY	BAIT_DESIGN_EFFICIENCY	TOTAL_READS	PF_READS	PF_UNIQUE_READS	PCT_PF_READS	PCT_PF_UQ_READS	PF_UQ_READS_ALIGNED	PCT_PF_UQ_READS_ALIGNED	PF_UQ_BASES_ALIGNED	ON_BAIT_BASES	NEAR_BAIT_BASES	OFF_BAIT_BASES	ON_TARGET_BASES	PCT_SELECTED_BASES	PCT_OFF_BAIT	ON_BAIT_VS_SELECTEDMEAN_BAIT_COVERAGE	MEAN_TARGET_COVERAGE	PCT_USABLE_BASES_ON_BAIT	PCT_USABLE_BASES_ON_TARGET	FOLD_ENRICHMENT	ZERO_CVG_TARGETS_PCT	FOLD_80_BASE_PENALTY	PCT_TARGET_BASES_2X	PCT_TARGET_BASES_10X	PCT_TARGET_BASES_20X	PCT_TARGET_BASES_30X	HS_LIBRARY_SIZE	HS_PENALTY_10X	HS_PENALTY_20X	HS_PENALTY_30X	AT_DROPOUT	GC_DROPOUT
			baits	3137161264	1343846	1343846	1	7738802	7738802	7738802	1	1	6956673	0.898934	699877964	203242146	34186499	462449319	203242146	0.339243	0.660757	0.856014	151.239164	181.735888	0.260027	0.260027	677.920541	0.013615	3.365479	0.803834	0.779288	0.7524	0.727413		0	0	0	0	0

0	SAMPLE
1	LIBRARY
2	READ_GROUP
3	BAIT_SET
4	GENOME_SIZE
5	BAIT_TERRITORY
6	TARGET_TERRITORY
7	BAIT_DESIGN_EFFICIENCY
8	TOTAL_READS
9	PF_READS
10	PF_UNIQUE_READS
11	PCT_PF_READS
12	PCT_PF_UQ_READS
13	PF_UQ_READS_ALIGNED
14	PCT_PF_UQ_READS_ALIGNED
15	PF_UQ_BASES_ALIGNED
16	ON_BAIT_BASES
17	NEAR_BAIT_BASES
18	OFF_BAIT_BASES
19	ON_TARGET_BASES
20	PCT_SELECTED_BASES
21	PCT_OFF_BAIT
22	ON_BAIT_VS_SELECTED
23	MEAN_BAIT_COVERAGE
24	MEAN_TARGET_COVERAGE
25	PCT_USABLE_BASES_ON_BAIT
26	PCT_USABLE_BASES_ON_TARGET
27	FOLD_ENRICHMENT
28	ZERO_CVG_TARGETS_PCT
29	FOLD_80_BASE_PENALTY
30	PCT_TARGET_BASES_2X
31	PCT_TARGET_BASES_10X
32	PCT_TARGET_BASES_20X
33	PCT_TARGET_BASES_30X
34	HS_LIBRARY_SIZE
35	HS_PENALTY_10X
36	HS_PENALTY_20X
37	HS_PENALTY_30X
38	AT_DROPOUT
39	GC_DROPOUT

=cut

my $header_printed = 0;
my $paste_file_names;

# input is like /path/SID4695105858-2-TGACCA

foreach my $cov (@ARGV) {

    open F,  $cov . ".dedup.bam.target_coverage" || die "Error open $cov";
    open F2, $cov . ".dedup.bam.flagstat"        || die "Error open $cov";
    open F3, $cov . ".sorted.bam.flagstat"        || die "Error open $cov";

    open O, ">" . $cov . '.out';

    $paste_file_names .= "$cov.out ";
    my $line;

    #skip headers
    while ( $line =~ /^#/ || $line =~ /^$/ ) {
        $line = <F>;
    }

    my @header = split /\s+/, $line;
    foreach ( 0 .. $#header ) {

        #	 print $_,"\t",$header[$_],"\n";
    }

## process flagstat
# since I choose to remove all the dups, have to change this code
#    <F2>;
#    my $line2 = <F2>;
#    my $line3 = <F2>;
#    my ( $dup,    @t )  = split /\s+/, $line2;
#    my ( $mapped, @t2 ) = split /\s+/, $line3;

# to
	my ( $total , @t) = split/\s+/, <F3>;
	<F3>;
	my ( $total_mapped, @t ) = split/\s+/, <F3>;
	my ( $after_dup, @t ) = split/\s+/, <F2>;
	<F2>;
	my ( $after_dup_mapped, @t ) = split/\s+/, <F2>;

    $line = <F>;
    my @result = split /\s+/, $line;
    #$cov =~ /(.*?).dedup/;
    #$result[0] = $1;
    $result[0] = $cov;
	foreach ( 2, 4, 6, 8, 15, 17, 19, 24, 28, 31, 32 ) {
        if ( !$header_printed ) {
            print O $header[$_], "\t";
        }
        if ( $header[$_] =~ /^PCT/ || $header[$_] =~ /PCT$/ ) {
            print O sprintf( "%.2f", $result[ $_ - 2 ] * 100 ), "\n";
        }
        else {
            if ( $result[ $_ - 2 ] =~ /\d\.\d/ ) {
                print O sprintf( "%.2f", $result[ $_ - 2 ] ), "\n";
            }
            else {
                print O $result[ $_ - 2 ], "\n";
            }
        }
    }
    if ( !$header_printed ) {
	print O "Total_reads\t", $total, "\n";
        print O "Duplicate_reads\t", $total - $after_dup,    "\n";
        print O "Total_mapped_reads\t",    $total_mapped, "\n";
        print O "Duplicate_reads\/Mapped_reads\t",
          sprintf( "%.2f", ($total - $after_dup) / $total_mapped );
        print O "Mapped_reads_after_dup\t",    $after_dup_mapped, "\n";

    }
    else {
        print O $dup,    "\n";
        print O $mapped, "\n";
        print O sprintf( "%.2f", $dup / $mapped );

    }

    close F;
    close O;
    close F2;
	close F3;

    $header_printed = 1;
}

`paste $paste_file_names > target.coverage`;

