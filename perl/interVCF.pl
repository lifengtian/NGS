
=head1 NOTE

    Find the intersection and difference between two VCF files

=cut



my ($f1, $f2, $one,$both,$two) = @ARGV;

open F1, $f1;
open F2, $f2;
open ONE, ">$one";
open BOTH,">$both";
open TWO,">$two";

while(<F1>){
	next if /^#/ ;
	my @a=split;
	$h{$a[0].':'.$a[1]} = $_;

}

while(<F2>){
	next if /^#/ ;
	my @a=split;
	 if ( $h{$a[0].':'.$a[1]} ) {
		print BOTH $h{$a[0].':'.$a[1]} ;
	} elsif ( ! $h{$a[0].':'.$a[1]} )  {
		print TWO $_;
	}

}

close F1;
close F2;



open F1, $f1;
open F2, $f2;

while(<F2>){
        next if /^#/ ;
        my @a=split;
        $g{$a[0].':'.$a[1]} = $_;

}

while(<F1>){
        next if /^#/ ;
        my @a=split;
         if ( !$g{$a[0].':'.$a[1]} ) {
        	print ONE $_;
	}

}


close ONE;
close BOTH;
close TWO;
