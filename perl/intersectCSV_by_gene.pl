use strict;
use warnings;
use Text::CSV;
use List::Compare;

=head1 note

=cut

my $count1;
my $count2;
my $count3;

my $in1    = $ARGV[0];
my $in2    = $ARGV[1];
my $first_fn  = $ARGV[2];
my $all_fn    = $ARGV[3];
my $second_fn = $ARGV[4];

open (my $F1, $in1);
open (my $F2, $in2);

open (my $FIRST,  ">$first_fn");
open (my $ALL,    ">$all_fn");
open (my $SECOND, ">$second_fn");

my %h1;
my %h;
my %h2;

my $csv = Text::CSV->new();

<$F1>;
<$F2>;

while (<$F1>) {
    chomp;
    if ( $csv->parse($_) ) {
        my @a = $csv->fields();
        $h1{ $a[1] }{ join( "\t", @a[ 15 .. 19 ] ) } = $_;
    }
}
close $F1;

while (<$F2>) {
    chomp;
    if ( $csv->parse($_) ) {
        my @a = $csv->fields();
        $h2{ $a[1] }{ join( "\t", @a[ 15 .. 19 ] ) } = $_;
    }
}
close $F2;

#find keys shared or unique
my @k1 = keys %h1;
my @k2 = keys %h2;

my $lc = List::Compare->new( \@k1, \@k2 );

my @l_inter = $lc->get_intersection();
my @l_only1 = $lc->get_Lonly();
my @l_only2 = $lc->get_Ronly();

foreach (@l_inter) {
    foreach my $k ( keys %{ $h1{$_} } ) {
        print $ALL $h1{$_}->{$k}, "\n";
    }
    foreach my $k ( keys %{ $h2{$_} } ) {
        print $ALL $h2{$_}->{$k}, "\n";
    }
}

foreach (@l_only2) {
    foreach my $k ( keys %{ $h2{$_} } ) {
        print $SECOND $h2{$_}->{$k}, "\n";
    }
}

foreach (@l_only1) {
    foreach my $k ( keys %{ $h1{$_} } ) {
        print $FIRST $h1{$_}->{$k}, "\n";
    }
}

close $ALL;
close $FIRST;
close $SECOND;

