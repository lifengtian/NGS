#!/bin/sh

source $HISEQ/NGS/shell/setup.sh

###### FUNCTION    check_files ######
# check the locations of required files
#
function check_files {

check_list=( $GENOMECOVERAGE \
             $SAMTOOLS \
             $BWA \
             $HG19 \
             $JAVA_BIN \
             $FASTQC/fastqc \ 
             $PICARD/MarkDuplicates.jar
)

bye='no'
for i in ${check_list[*]}; do
if [ ! -f $i ]; then
    echo $i not exist
        bye='yes'
fi
done

echo bye is $bye

if [ "$bye" = "yes" ]; then
        echo See you next time
        exit
fi

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

mkdir -p $p/log
mkdir -p $p/scripts
mkdir -p $p/temp
mkdir -p $p/QCreport

#align
timestamp=`date +%s`
echo "
$BWA aln  -t 4 $HG19 $r1.fastq.gz > $r1.sai
$BWA aln  -t 4 $HG19 $r2.fastq.gz > $r2.sai
" > $p/scripts/$sm-$lane-$index.aln.sh


job=`qsub $queue -pe smp 4 -e $p/log/$sm-$lane-$index.aln.err.log -o $p/log/$sm-$lane-$index.aln.log $p/scripts/$sm-$lane-$index.aln.sh`
jobs=`echo $job | awk '{print $3}'`

hold_jid=$jobs

echo "
$BWA sampe -r \"@RG\\tID:$fid-$lane-$index-$sm\\tSM:$sm\" $HG19 $r1.sai $r2.sai $r1.fastq.gz $r2.fastq.gz  | $SAMTOOLS view -S -b -o $p/$sm-$lane-$index.bam -
$JAVA_BIN -jar $PICARD/SortSam.jar I=$p/$sm-$lane-$index.bam O=$p/$sm-$lane-$index.sorted.bam SO=coordinate VALIDATION_STRINGENCY=SILENT TMP_DIR=$p/temp CREATE_INDEX=true
$SAMTOOLS flagstat $p/$sm-$lane-$index.sorted.bam > $p/$sm-$lane-$index.sorted.bam.flagstat
$SAMTOOLS depth $p/$sm-$lane-$index.sorted.bam | perl $GENOMECOVERAGE > $p/$sm-$lane-$index.sorted.bam.genomecoverage
$FASTQC/fastqc -o $p/QCreport  $p/$sm-$lane-$index.sorted.bam

$JAVA_BIN -jar $PICARD/MarkDuplicates.jar I=$p/$sm-$lane-$index.sorted.bam O=$p/$sm-$lane-$index.dedup.bam M=$p/$sm-$lane-$index.metric VALIDATION_STRINGENCY=SILENT TMP_DIR=$p/temp CREATE_INDEX=true REMOVE_DUPLICATES=true
$SAMTOOLS flagstat $p/$sm-$lane-$index.dedup.bam > $p/$sm-$lane-$index.dedup.bam.flagstat
$SAMTOOLS depth $p/$sm-$lane-$index.dedup.bam | perl $GENOMECOVERAGE > $p/$sm-$lane-$index.dedup.bam.genomecoverage
" > $p/scripts/$sm-$lane-$index.sampe.sh

job=`qsub $queue -hold_jid $hold_jid -e $p/log/$sm-$lane-$index.sampe.err.log -o $p/log/$sm-$lane-$index.sample.log $p/scripts/$sm-$lane-$index.sampe.sh`

}

######## FUNCTION feed_bwa #########
# find all fastq.gz pairs 
# create Aligned folder
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

# mkdir Aligned folder
    mkdir -p $bam/Project_$SampleProject/Sample_$SampleID
    newr1=$bam/Project_$SampleProject/Sample_$SampleID/${SampleID}_${Index}_L00${Lane}_R1_001
    newr2=$bam/Project_$SampleProject/Sample_$SampleID/${SampleID}_${Index}_L00${Lane}_R2_001
    ln -s $r1.fastq.gz $newr1.fastq.gz
    ln -s $r2.fastq.gz $newr2.fastq.gz

bwa $bam/Project_$SampleProject/Sample_$SampleID $newr1 $newr2 $SampleID $FCID $Lane $Index
#sleep 10
fi # end_of_count -gt 1

done # end_of_while read sample_sheet

}


#feed_bwa FlowCellID
feed_bwa    $1

