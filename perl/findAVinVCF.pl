

my ($av, $vcf, $out) = @ARGV;

open F1, $av;
open F2, $vcf;
open OUT, ">$out";

while(<F1>){
        my @a=split;
        $hav{$a[0].':'.$a[1]} = $_;

}

while(<F2>){
        next if /^#/ ;
        my @a=split;
         if ( $hav{$a[0].':'.$a[1]} ) {
                print OUT $_; 
        }

}

close F1;
close F2;
close OUT;


