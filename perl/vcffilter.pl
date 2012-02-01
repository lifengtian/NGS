#!/bin/env perl

use strict;

my $usage = "perl vcffilter.pl \n
		[vcf filename] \n
		[chr default 'ALL'] \n
		[qual 30] \n
		[filter default '.'] \n
		[info_filter 1000g 0.05 lt,AB 0.8 lt] \n
		[Sample name, default is 'ALL'] \n
		[het, default is 'both'] \n
	";

if ( @ARGV == 0 ) {
    die $usage;
}

my ( $fn, $chr, $quality, $filter_in, $info, $sample, $het ) = @ARGV;

my $float_n = qr'[-+]?[0-9]*\.?[0-9]+';
my $int_n   = qr'\d+';
my $string  = qr'.+';

my %field_name = (
    'avHet'            => $float_n,
    '1000g'            => $float_n,
    'DP'               => $int_n,
    'AB'               => $float_n,
    'BaseQRankSum'     => $float_n,
    'AC'               => $int_n,
    'HaplotypeScore'   => $float_n,
    'FS'               => $float_n,
    'MQ0'              => $int_n,
    'QD'               => $float_n,
    'AF'               => $float_n,
    'ReadPosRankSum'   => $float_n,
    'AN'               => $int_n,
    'MQ'               => $float_n,
    'HRun'             => $float_n,
    'MQRankSum'        => $float_n,
    'MAFinPercent_EA'  => $float_n,
    'MAFinPercent_AA'  => $float_n,
    'MAFinPercent_ALL' => $float_n,
    'ExonicFunc'       => $string,
    'dbSNP132'         => $string,
    'CG69'             => $float_n,
);

my %field_op = (
    'avHet'            => [ 'lt', 1, 'gt', 1 ],
    '1000g'            => [ 'lt', 1, 'gt', 1 ],
    'DP'               => [ 'lt', 1, 'gt', 1 ],
    'AB'               => [ 'lt', 1, 'gt', 1 ],
    'BaseQRankSum'     => [ 'lt', 1, 'gt', 1 ],
    'AC'               => [ 'lt', 1, 'gt', 1 ],
    'HaplotypeScore'   => [ 'lt', 1, 'gt', 1 ],
    'FS'               => [ 'lt', 1, 'gt', 1 ],
    'MQ0'              => [ 'lt', 1, 'gt', 1 ],
    'QD'               => [ 'lt', 1, 'gt', 1 ],
    'AF'               => [ 'lt', 1, 'gt', 1 ],
    'ReadPosRankSum'   => [ 'lt', 1, 'gt', 1 ],
    'AN'               => [ 'lt', 1, 'gt', 1 ],
    'MQ'               => [ 'lt', 1, 'gt', 1 ],
    'HRun'             => [ 'lt', 1, 'gt', 1 ],
    'MQRankSum'        => [ 'lt', 1, 'gt', 1 ],
    'MAFinPercent_EA'  => [ 'lt', 1, 'gt', 1 ],
    'MAFinPercent_AA'  => [ 'lt', 1, 'gt', 1 ],
    'MAFinPercent_ALL' => [ 'lt', 1, 'gt', 1 ],
    'ExonicFunc'       => [ 'eq', 1, 'ne', 1 ],
    'dbSNP132'         => [ 'eq', 1, 'ne', 1 ],
    'CG69'             => [ 'lt', 1, 'gt', 1 ],
);

my @ops;
my $info_p = 1;    #default setting requires some filtering string.

#deal with the awful space issue, e.g., '     DP 10 lt,   AB  0.7 gt,   '
$info =~ s/^\s*//g;
$info =~ s/\s*$//g;
if ( $info =~ /,/ ) {
    $info =~ s/,$//g;

    foreach my $i ( split /,/, $info ) {
        $i =~ s/^\s*//g;
        my ( $field, $value, $op ) = split /\s+/, $i;

        my %op_allowed = @{ $field_op{$field} };

        if ( !$field_name{$field} ) {
            die "not known field. exit";
        }

        elsif ( !$op_allowed{$op} ) {
            die "op '$op' not know. exit";
        }

        elsif ( $value !~ /^$field_name{$field}$/ ) {
            die "value $value is not known for", $field_name{$field}, ". exit";
        }

        push @ops, [ $field, $value, $op ];
    }
}
else {
    $info_p = 0;
}


# chromosome
my %chr_h;
if ( $chr ne 'ALL' ) {
    my @a = split /,/, $chr;
    map { $chr_h{$_}++ } @a;
}
else {
    my @chr_a =
      qw/chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM/;
    map { $chr_h{$_}++ } @chr_a;
}

# samples
my %sample_h;
my @samples_to_output = ( 0 .. 8 );
my @samples_to_keep;

if ( $sample ne 'ALL' ) {
    my @b = split /,/, $sample;
    map { $sample_h{$_}++ } @b;
}

my $sample_no;

open my $in, $fn or die "error open $fn.";
LOOP: while ( my $line = <$in> ) {

    if ( $line =~ /^##/ ) {
        print $line;
        next;
    }

    # CHROM
    if ( $line =~ /^#CHROM/ ) {
        my @a = split /\t/, $line;
        $sample_no = $#a - 8;
        if ( $sample ne 'ALL' ) {
            if ( $sample_no == 1 ) {
                if ( $sample_h{ $a[9] } ) {
                    push @samples_to_output, 9;
                }
            }
            else {    #more than 1 sample
                foreach my $i ( 9 .. $#a ) {
                    if ( $sample_h{ $a[$i] } ) {
                        push @samples_to_output, $i;
                        push @samples_to_keep,   $i;
                    }
                }
            }
        }
        else {
            @samples_to_output = 0 .. $#a;
            if ( $sample_no == 1 ) {
                push @samples_to_keep, 9;
            }
            else {
                foreach my $i ( 9 .. $#a ) {
                    push @samples_to_keep, $i;
                }
            }
        }
        print join( "\t", @a[@samples_to_output] );
    }


    my @lines = split /\t/, $line;
    my ( $qual, $filter, $info_f ) = @lines[ 5, 6, 7 ];

    my @info1 = split /;/, $info_f;

    #    my %info_h = map { split /=/ } @info1;
    my %info_h;

    # info string
    foreach my $i (@info1) {
        if ( $i =~ /=/ ) {
            my ( $a, $b ) = split /=/, $i;
            if ( defined $a && defined $b ) {
                $info_h{$a} = $b;
            }
            else {
                die "$a or $b is null!";
            }
        }
    }

    # qual
    if ( $qual < $quality ) {
        next LOOP;
    }
    if ( $filter_in ne '.' ) {
        if ( $filter ne $filter_in ) {
            next LOOP;
        }
    }

    # het, hom or both
    # for both: do nothing
    # WARNING: 1. did not exclude those ./.;
    #          2. did not work with phased 0|0;
    #          3. did not work with 1/0;
    my @s = map { ( split( /:/, $lines[$_] ) )[0] } @samples_to_keep;
    foreach my $i (@s) {
        if ( $het eq 'het' && $i ne '0/1' ) {
            next LOOP;
        }
        elsif ( $het eq 'hom' && $i ne '1/1' ) {
            next LOOP;
        }
    }

    foreach my $opl (@ops) {
        my ( $f, $v, $o ) = @{$opl};
        if ( defined $info_h{$f} ) {
            if ( $o eq 'lt' ) {
                if ( $info_h{$f} > $v ) {
                    next LOOP;
                }
            }
            elsif ( $o eq 'gt' ) {
                if ( $info_h{$f} < $v ) {
                    next LOOP;
                }

            }

            # 1. not-matching records will be excluded
            elsif ( $o eq 'eq' ) {
                if ( $info_h{$f} !~ /^$v/ ) {

                    #print $info_h{$f},"\t",$v,"\n";
                    next LOOP;
                }

            }
            elsif ( $o eq 'ne' ) {
                if ( $info_h{$f} =~ /^$v/ ) {
                    next LOOP;
                }
            }
        }    #print "$f passed $info_h{$f} over standard $v by $o\n";
        elsif ( $o eq 'eq' || $o eq 'ne' ) {
            next LOOP;
        }
    }
    chomp(@lines);
    if ( $sample eq 'ALL' ) {
        print join( "\t", @lines ), "\n";
    }
    else {
        print join( "\t", @lines[@samples_to_output] ), "\n";
    }
}

close $in

