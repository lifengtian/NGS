#!/bin/bash


USAGE="demultiplex.sh [FlowCellID]"

source $HISEQ/setup.sh


FlowCellID=$1
FlowCellCSV=$FlowCellInfo/$FlowCellID.csv 
# input_dir is retrieved from FlowCellInfo/ID.csv
if [ ! -f $FlowCellCSV ]; then
	echo $FlowCellCSV not exist
	exit
fi


input_dir=`head -1 $FlowCellCSV`

#
output_dir=$HISEQ/FASTQ/$FlowCellID

if [ -d $output_dir ]; then
	echo $output_dir exist
	exit
fi

# SampleSheet
SampleSheet=$HISEQ/SampleSheet/$FlowCellID.csv

if [ ! -f $SampleSheet ]; then
        echo $SampleSheet not exist
        exit
fi


#
echo "

$BclToFastq --fastq-cluster-count 0 --input-dir  $input_dir  --output-dir $output_dir --sample-sheet $SampleSheet
cd $output_dir
nohup make -j 4  
" > $HISEQ/FASTQ/scripts/$FlowCellID.sh
qsub $queue -V -e $HISEQ/logs/$FlowCellID.log -o $HISEQ/logs/$FlowCellID.log $HISEQ/FASTQ/scripts/$FlowCellID.sh

