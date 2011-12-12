use Text::CSV;


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
		$h1{join("\t",@a[15..19])}=$_;
		$h{join("\t",@a[15..19])}++;
	}
}
close F1;

while(<F2>){
	chomp;
	if ( $csv->parse($_) ) {
        my @a=$csv->fields(); 
		$h2{join("\t",@a[15..19])}=$_;
		$h{join("\t",@a[15..19])} += 2;
	}
}
close F2;

foreach (keys %h){
     if ( $h{$_} == 3 ) {
		print ALL $h1{$_},"\n";
	}                          elsif ( $h{$_} == 1 ) {
		print FIRST $h1{$_},"\n";
	} elsif ( $h{$_} == 2 ) {
		print SECOND $h2{$_},"\n";
	} else {
		#print $h1{$_},"\n";
		#print $h2{$_},"\n";
		#die "$h{$_} error";
	}
}

close ALL;
close FIRST;
close SECOND;


