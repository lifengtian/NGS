use Text::CSV;

=head1 note

gene may have multiple lines. 
=cut

my $count1;
my $count2;
my $count3;



my $in1=$ARGV[0];
my $in2=$ARGV[1];
my $first=$ARGV[2];
my $all=$ARGV[3];
my $second=$ARGV[4];

open F1,$in1;
open F2,$in2;

open FIRST,">$first";
open ALL ,">$all";
open SECOND,">$second";

my %h1;
my %h;
my %h2;

my $csv = Text::CSV->new();

<F1>;
<F2>;

while(<F1>){
	chomp;
	if ( $csv->parse($_) ) {
        my @a=$csv->fields();	
		$h1{$a[1]}{join("\t",@a[15..19])}=$_;
		$h{$a[1]}=1;
		$count1++;
	}
}
close F1;

while(<F2>){
	chomp;
	if ( $csv->parse($_) ) {
        my @a=$csv->fields(); 
		$h2{$a[1]}{join("\t",@a[15..19])}=$_;
		if ( $h{$a[1]} == 1 ) { $h{$a[1]} = 3;$count3++ } else {
			$h{$a[1]} = 2 ; $count2++
		}
	}
}
close F2;

print join("\t",$count1,$count2,$count3),"\n";

foreach (keys %h){
     if ( $h{$_} == 3 ) {
		
		foreach my $k ( keys %{$h1{$_}} ) {
			print ALL $h1{$_}->{$k},"\n"; 
		}
		foreach my $k ( keys %{$h2{$_}} ) {
                        print ALL $h2{$_}->{$k},"\n";
                }
	} elsif ( $h{$_} == 2 ) {
		foreach my $k ( keys %{$h2{$_}} ) {
                        print SECOND $h2{$_}->{$k},"\n";
                }
	} elsif ($h{$_} == 1) {
		foreach my $k ( keys %{$h1{$_}} ) {
                        print FIRST $h1{$_}->{$k},"\n";
                }
	} else {
		#print $h1{$_},"\n";
		#print $h2{$_},"\n";
		die "$h{$_} error";
	}
}

close ALL;
close FIRST;
close SECOND;


