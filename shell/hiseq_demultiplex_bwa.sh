#!/bin/bash

#############################################
# Name:   hiseq_demultiplex_bwa.sh
# Input:  flow cell ID
# Output: *.dedup.indelrealigner.recal.bam
# version:1.0
# 
# Note:   This script only works on resolute
#         
##############################################

source $HISEQ/NGS/shell/setup.sh
source $HISEQ/NGS/shell/function.sh

hold_jid=''
platform=' -dP illumina '
ref=$HG19
dbsnp=$GATK/b37/dbsnp_132.b37.vcf

Target=$GATK/bed/SureSelect50mbclean
Target_bed=$Target.bed
Target_picard=$Target.picard


####### FUNCTION demultiplex ############
function demultiplex {

USAGE="demultiplex.sh [FlowCellID]"


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
job=qsub $queue -V -e $HISEQ/logs/$FlowCellID.log -o $HISEQ/logs/$FlowCellID.log $HISEQ/FASTQ/scripts/$FlowCellID.sh
jobs=`echo $job | awk '{print $3}'`
hold_jid=$jobs

}






######### FUNCTION bwa ###########
# run BWA
##################################

function bwa {
#$bam/Project_$SampleProject/Sample_$SampleID
p=$1

#forward and reverse reads FASTQ(gzipped FASTQ) file name
r1=$2
r2=$3

if [ ! -f $r1.fastq.gz ]; then
  echo $r1.fastq.gz not exist!
  exit
fi
if [ ! -f $r2.fastq.gz ]; then
  echo $r2.fastq.gz not exist!
  exit
fi



#sampleID 
sm=$4
#flowcellID
fid=$5
lane=$6
index=$7

i="$sm-$lane-$index" #for gatk cleaning-up

mkdir -p $p/log
mkdir -p $p/scripts
mkdir -p $p/temp
mkdir -p $p/QCreport

#align
timestamp=`date +%s`
echo "
$BWA aln -n 3 -t 4 $HG19 $r1.fastq.gz > $r1.sai
$BWA aln -n 3 -t 4 $HG19 $r2.fastq.gz > $r2.sai
" > $p/scripts/$i.aln.sh


job=`qsub $queue -pe smp 4 -hold_jid=$hold_jid -e $p/log/$i.aln.err.log -o $p/log/$i.aln.log $p/scripts/$i.aln.sh`
jobs=`echo $job | awk '{print $3}'`

hold_jid=$jobs

echo "
$BWA sampe -r \"@RG\\tID:$fid-$lane-$index-$sm\\tSM:$sm\" $HG19 $r1.sai $r2.sai $r1.fastq.gz $r2.fastq.gz  | $SAMTOOLS view -S -b -o $p/$i.bam -
$JAVA_BIN -jar $PICARD/SortSam.jar I=$p/$i.bam O=$p/$i.sorted.bam SO=coordinate VALIDATION_STRINGENCY=SILENT TMP_DIR=$p/temp CREATE_INDEX=true
$SAMTOOLS flagstat $p/$i.sorted.bam > $p/$i.sorted.bam.flagstat
$SAMTOOLS depth $p/$i.sorted.bam | perl $GENOMECOVERAGE > $p/$i.sorted.bam.genomecoverage
$FASTQC/fastqc -o $p/QCreport  $p/$i.sorted.bam

$JAVA_BIN -jar $PICARD/MarkDuplicates.jar I=$p/$i.sorted.bam O=$p/$i.dedup.bam M=$p/$i.metric VALIDATION_STRINGENCY=SILENT TMP_DIR=$p/temp CREATE_INDEX=true REMOVE_DUPLICATES=true
$SAMTOOLS flagstat $p/$i.dedup.bam > $p/$i.dedup.bam.flagstat
$SAMTOOLS depth $p/$i.dedup.bam | perl $GENOMECOVERAGE > $p/$i.dedup.bam.genomecoverage

#coverage code can be inserted here?

### gatk BAM cleaning-up
    # Local realignment around Indels
    # 1. Prepare Realign Interval file
    $JAVA_GATK \
        -R $ref \
        -o $p/$i.dedup.intervals \
        -I $p/$i.dedup.bam \
        -T RealignerTargetCreator

    # 2. Indel Realigner
    $JAVA_GATK  \
        -I $p/$i.dedup.bam \
        -R $ref \
        -targetIntervals $p/$i.dedup.intervals \
        -o $p/$i.dedup.indelrealigner.bam \
        -T IndelRealigner


    # Base quality score recalibration
    # 1.Count covariates

        # CountCovariates requires indexed BAMs
    $SAMTOOLS index $p/$i.dedup.indelrealigner.bam

    $JAVA_GATK -R $ref  \
        -knownSites $dbsnp \
        -I $p/$i.dedup.indelrealigner.bam \
        --covariate ReadGroupCovariate \
        --covariate QualityScoreCovariate \
        --covariate CycleCovariate \
        --covariate DinucCovariate \
        -recalFile $p/$i.dedup.indelrealigner.recal.csv \
        $platform \
        -T CountCovariates \


    # 2. Table Recalibration
    $JAVA_GATK -R $ref  \
        -I $p/$i.dedup.indelrealigner.bam \
        --out $p/$i.dedup.indelrealigner.recal.bam \
        -recalFile $p/$i.dedup.indelrealigner.recal.csv \
        $platform \
        --default_read_group defaultRG \
        -T TableRecalibration \

    $SAMTOOLS index $p/$i.dedup.indelrealigner.recal.bam



" > $p/scripts/$i.sampe.sh

if [ "$desc" = "WES" ]; then
	echo "
	$JAVA_BIN -jar $PICARD/CalculateHsMetrics.jar I=$p/$i.dedup.bam O=$p/$i.dedup.bam.target_coverage BI=$Target_picard TI=$Target_picard TMP_DIR=$p/temp VALIDATION_STRINGENCY=SILENT
	$SAMTOOLS depth -b $Target_bed $p/$i.sorted.bam | perl $GENOMECOVERAGE > $p/$i.sorted.bam.target_depth
	$SAMTOOLS depth -b $Target_bed $p/$i.dedup.bam | perl $GENOMECOVERAGE > $p/$i.dedup.bam.target_depth

" >> $p/scripts/$i.sampe.sh
fi

job=`qsub $queue -hold_jid $hold_jid -e $p/log/$i.sampe.err.log -o $p/log/$i.sample.log $p/scripts/$i.sampe.sh`
} #end_of_bwa

######## FUNCTION feed_bwa #########
# find all fastq.gz pairs 
# create BAM folder
# prepare calls for bwa 
####################################

function feed_bwa {
flow_cell_id=$1
sample_sheet=$HISEQ/SampleSheet/$flow_cell_id.csv
fastq=$HISEQ/FASTQ/$flow_cell_id
bam=$HISEQ/BAM/$flow_cell_id



#check files
check_files

if [ ! -d $bam ]; then
    mkdir -p $bam
else
#
	echo "$bam exist. Exit"
	exit
fi

# parse sample_sheet
#FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject
#C0995ACXX,1,lane1,Unknown,,'DefaultSample',N,,,C0995ACXX
count=0
cat $sample_sheet | while read LINE
do
let count++
if [ $count -gt 1 ]; then
echo "$count $LINE"
FCID=`echo $LINE | cut -f1 -d','`
Lane=`echo $LINE | cut -f2 -d','`
SampleID=`echo $LINE | cut -f3 -d','`
Index=`echo $LINE | cut -f5 -d','`
SampleProject=`echo $LINE | cut -f10 -d','`
description=`echo $LINE | cut -f6 -d','`

if [ ! $Index ]; then
Index="NoIndex"
fi

#check existence of fastq.gz
r1=$fastq/Project_$SampleProject/Sample_$SampleID/${SampleID}_${Index}_L00${Lane}_R1_001
r2=$fastq/Project_$SampleProject/Sample_$SampleID/${SampleID}_${Index}_L00${Lane}_R2_001
if [ ! -f $r1.fastq.gz ]; then
	echo $r1.fastq.gz not exisit. Exit.
#exit
fi

if [ ! -f $r2.fastq.gz ]; then
	echo $r2.fastq.gz not exisit. Exit.
#exit
fi

# mkdir BAM folder
    mkdir -p $bam/Project_$SampleProject/Sample_$SampleID
    newr1=$bam/Project_$SampleProject/Sample_$SampleID/${SampleID}_${Index}_L00${Lane}_R1_001
    newr2=$bam/Project_$SampleProject/Sample_$SampleID/${SampleID}_${Index}_L00${Lane}_R2_001
    ln -s $r1.fastq.gz $newr1.fastq.gz
    ln -s $r2.fastq.gz $newr2.fastq.gz

bwa $bam/Project_$SampleProject/Sample_$SampleID $newr1 $newr2 $SampleID $FCID $Lane $Index $description
#sleep 10
fi # end_of_count -gt 1

done # end_of_while read sample_sheet

} # end_of_function_feed_bwa

demultiplex $1
#feed_bwa FlowCellID
feed_bwa    $1

