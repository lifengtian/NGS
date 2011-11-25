
=head1 NOTE
   genomecoverage.pl --- Calculate coverage for whole genome sequencing

   INPUT: samtools depth output
	chr1	9994	1
	chr1	9995	4
	chr1	9996	4
	chr1	9997	4
	chr1	9998	4
	chr1	9999	4
	chr1	10000	4
	chr1	10001	6
	chr1	10002	7
	chr1	10003	9

   OUTPUT: 
        chr	bases	1	4	6	7	9	max
	chr1	47	1	6	1	1	1	9
	sum= 47; pos=10


=cut

my $VERSION=0.1;

my $total_bases = 0;
my $covered_bases = 0;

while(<>){
        my @a=split;
	my $chr=$a[0];
	my $coverage=$a[2];
        $total_bases += $coverage;
        $covered_bases++;

        $bases_per{$chr}+= $coverage;

        $sum1{$chr}{$coverage}++;
        $cov{$coverage}++;

        if ( $coverage  >= $max{$chr} ) {
        $max{$chr} = $coverage;
    }

                if ( $coverage < 1 ) {
                #print $_;
        }
}

my @bins=sort {$a<=>$b} keys %cov;
print join("\t","chr","bases",join("\t",@bins),"max"),"\n";

foreach my $k ( keys %sum1) {
print join("\t",$k,$bases{$k});
foreach my $n ( @bins  ) {
        if ( $sum1{$k}{$n} ) {
                print "\t",$sum1{$k}{$n};
        } else {
                print "\t0";
        }
}
print "\t$max{$k}\n";

}

print "$total_bases  $total_bases covered_bases $covered_bases\n";
