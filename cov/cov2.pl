my $target_region = 0;
my $bases_in_target = 0;

while(<>){
	my @a=split/\s+/;
	$target_region+=($a[2] - $a[1] )+ 1;
	$bases_in_target += $a[3];
}

print $target_region,"\t",$bases_in_target,"\t",$bases_in_target / $target_region ,"\n";

