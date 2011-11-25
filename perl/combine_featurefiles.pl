#!/usr/bin/perl -w 

### combine three feature files 
### extract only the gene normalized expression value
### Lifeng Tian
### April 2010
###

my $usage="perl combine_featurefiles.pl <feature file1>  <feature file2> <feature file3>";

my %geneid = ();
my $genecounter;
my @gene_loc = ();

my %gene_ave_nrm = ();
my %gene_ave_nrm2 = ();
my %gene_ave_nrm3 = ();


if ( @ARGV != 3 ){
	die $usage;
}


open(INFILE, $ARGV[0]);
<INFILE>;
<INFILE>; # this is the ----------- ... line
$genecounter = 0;
while( my $line = <INFILE>) { # this is the gene id line
    chomp($line);
    $geneid[$genecounter] = $line;

    $line = <INFILE>; # this is the header line that each gene has
    $line = <INFILE>; # this is the line for the full gene (sum of all exons)
    chomp($line);
    $line =~ s/ //g;
    my @a = split(/\t/,$line);
    $gene_loc[$genecounter] = $a[1];
    $gene_ave_nrm[$genecounter] = $a[4];
    $line = <INFILE>; # this is the first exon line
    chomp($line);
    until(($line =~ /----------/) || !($line =~ /\S/)) {
	$line = <INFILE>;
    }
    $genecounter++;
}
close(INFILE);


open(INFILE, $ARGV[1]);
<INFILE>;
<INFILE>; # this is the ----------- ... line
$genecounter = 0;
while( my $line = <INFILE>) { # this is the gene id line
    chomp($line);
    $line = <INFILE>; # this is the header line that each gene has
    $line = <INFILE>; # this is the line for the full gene (sum of all exons)
    chomp($line);
    $line =~ s/ //g;
    my @a = split(/\t/,$line);
    $gene_ave_nrm2[$genecounter] = $a[4];
    $line = <INFILE>; # this is the first exon line
    chomp($line);
    until(($line =~ /----------/) || !($line =~ /\S/)) {
	$line = <INFILE>;
    }
    $genecounter++;
}
close(INFILE);



open(INFILE, $ARGV[2]);
<INFILE>;
<INFILE>; # this is the ----------- ... line
$genecounter = 0;
while(my $line = <INFILE>) { # this is the gene id line
    chomp($line);
    $line = <INFILE>; # this is the header line that each gene has
    $line = <INFILE>; # this is the line for the full gene (sum of all exons)
    chomp($line);
    $line =~ s/ //g;
    my @a = split(/\t/,$line);
    $gene_ave_nrm3[$genecounter] = $a[4];

    $line = <INFILE>; # this is the first exon line
    chomp($line);
    until(($line =~ /----------/) || !($line =~ /\S/)) {
	$line = <INFILE>;
    }
    $genecounter++;
}
close(INFILE);



print "GENE_ID\tLocation\t$ARGV[0]_intensity\t$ARGV[1]_intensity\t$ARGV[2]_intensity\n";
for( my $i=0; $i<$genecounter; $i++) {
     print "$geneid[$i]\t$gene_loc[$i]\t$gene_ave_nrm[$i]\t$gene_ave_nrm2[$i]\t$gene_ave_nrm3[$i]\n";
}





sub make_ratio (){
    my ($a, $b) = @_;
    if ( $a == 0 || $b == 0 ) {
	return 10000;
    } elsif ( $a > $b ) {
	return $a / $b ;
    } else {
	return $b / $a ;
    }

}
