package hiseq_pipeline;

use strict;

#use hiseq.pm;  #parameters
## global variabes 
my $hiseq=$ENV{'hiseq'};
my $FlowCellInfo=$hiseq.'/FlowCellInfo/';
my $BclToFastq=''; #CASAVA 1.8.2
my $FASTQ=$hiseq.'/FASTQ';

my $usage="perl hiseq_pipeline [flow cell ID 1] ... \n";

foreach my $fid ( @ARGV ){
    #assuming working on current $fid
    &run_pipeline($fid);
    
}


sub run_pipeline {
    my ($fid) = @_;

    # step 1
    demultiplex($fid); 
    bwa_aln($fid);
    bwa_sampe($fid);

    # step 2
    dedup($fid);
    bamstats($fid);
    indel_realign($fid);
    base_quality_recal($fid);
    cleanup_and_report($fid);
    return 1;
}

## Illumina HiSeq and CASAVA 1.8.2 requires a samplesheet file
## format of samplesheet

## demultiplex need
## 1. input_dir: 
## 2. output_dir:
## 3. SampleSheet file

sub demultiplex {
    my $FlowCellID=shift @_;   
    my $FlowCellCSV=$FlowCellInfo.'/'.$FlowCellID.'.csv';
    
    # input_dir is retrieved from FlowCellInfo/ID.csv
    if ( ! -f $FlowCellCSV ) {
        die "$FlowCellCSV not exist";
    }

    # since multiple folders may exist in a run (due to a stop etc.) 
    # it might not be safe to run it automatically for now
    # it requires an absolute path in the first line of this file in the FlowCellInfo folder
    my $input_dir=`head -1 $FlowCellCSV`;

    #
    my $output_dir=$FASTQ.'/'.$FlowCellID;

    if ( -d $output_dir ) {
        die  "$output_dir exist";
    }

    # SampleSheet
    my $SampleSheet=$hiseq.'/SampleSheet/'.$FlowCellID.'.csv';

    if ( ! -f $SampleSheet) {
        die "$SampleSheet not exist";
    }

    my $bcl2fastq="$BclToFastq --fastq-cluster-count 0 --input-dir  $input_dir  --output-dir $output_dir --sample-sheet $SampleSheet ";
    my $sh_out = $hiseq.'/scripts/'.$FlowCellID.'.step1.sh';
    open my $output, ">".$sh_out or die "Error open $sh_out";
    print $output " $bcl2fastq \n";
    print $output " cd $output_dir \n";
    print $output " nohup make -j 4 \n";
    close $output;
    #qsub $queue -V -e $HISEQ/logs/$FlowCellID.log -o $HISEQ/logs/$FlowCellID.log $HISEQ/FASTQ/scripts/$FlowCellID.sh

}


sub bwa_aln() {
        my $command
 }    
    
    
sub bwa_sampe() {}

# step 2
sub dedup() {}
sub bamstats() {}
sub indel_realign() {}
sub base_quality_recal(){}
sub cleanup_and_report(){}

