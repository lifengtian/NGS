#calculate the total length in a bed file
#Lifeng Tian

my $sum = 0;
while(<>){
     my @a=split/\s+/;
	$sum+= ($a[2]-$a[1]);  #bed is [a,b)
}

print "Length: $sum\n";

