my $sum = 0;
my $pos = 0;

while(<>){
	my @a=split/\s+/;
	$sum+=$a[2];
	$pos++;

	$bases{$a[0]}+=$a[2];

	if ( $a[2] >= 20 ) {
		$sum20{$a[0]}++; 
	} 
	if ( $a[2] >= 10 ) {
		$sum10{$a[0]}++;
	} 
	if ( $a[2] >= 8 ) {
		$sum8{$a[0]}++;
	} if ( $a[2] >= 4 ) {
		$sum4{$a[0]}++;
	} if ( $a[2] >=2 ) {
		$sum2{$a[0]}++;
	} if ( $a[2] >= 1 ) {
		$sum1{$a[0]}++;
	}
		if ( $a[2] < 1 ) {
		#print $_;
	}
}

print join("\t","chr","bases",20,10,8,4,2,1),"\n";

foreach my $k ( keys %sum1) {
print join("\t",$k,$bases{$k},$sum20{$k},$sum10{$k},$sum8{$k},$sum4{$k},$sum2{$k},$sum1{$k}),"\n";

}

print "sum= $sum; pos=$pos\n";

