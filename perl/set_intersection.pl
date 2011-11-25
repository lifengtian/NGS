#!/usr/bin/perl

 
my ($f1, $f2) = @ARGV;

open F1, $f1;
my @f1 = <F1>;
close F1;

open F2, $f2;
my @f2 = <F2>;
close F2;

chomp @f1; chop @f1;
chomp @f2; chop @f2;

shift @f1; shift @f1;
shift @f2; shift @f2;

foreach (@f1){
	$h{$_}++;
}

foreach (@f2){
	$h{$_}++;
}

foreach (keys %h){
	if ( $h{$_} == 2 ){
		print $_,"\n";
	}
}

