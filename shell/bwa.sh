#current path
p=`pwd`

# 
GENOMECOVERAGE=$GATK/NGS/perl/genomecoverage.pl


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
#forward and reverse reads FASTQ(gzipped FASTQ) file name
r1=$1
r2=$2

if [ ! -f $r1.fastq.gz ]; then
  echo $r1.fastq.gz not exist!
  exit
fi
if [ ! -f $r2.fastq.gz ]; then
  echo $r2.fastq.gz not exist!
  exit
fi



#sample name
sm=$3

mkdir -p $p/scripts
mkdir -p $p/log


#align
timestamp=`date +%s`
echo "
$BWA aln  -t 4 $HG19 $r1.fastq.gz > $r1.sai
$BWA aln  -t 4 $HG19 $r2.fastq.gz > $r2.sai

#$BWA sampe -r "@RG\tID:$sm-$timestamp\tSM:$sm" $HG19 $r1.sai $r2.sai $r1.fastq.gz $r2.fastq.gz  | $SAMTOOLS view -S -b -o $p/$sm.bam -
" > $p/scripts/bwa1_$sm.sh

job=`qsub -pe mpi 4 -e $p/log/bwa1_$sm.log -o $p/log/bwa1_$sm.log $p/scripts/bwa1_$sm.sh`
jobs[$count]=`echo $job | awk '{print $3}'`

        for i in ${jobs[*]}; do
                hold_jid=$hold_jid${i},
        done

echo "
$BWA sampe -r \"@RG\\tID:$sm-$timestamp\\tSM:$sm\" $HG19 $r1.sai $r2.sai $r1.fastq.gz $r2.fastq.gz  | $SAMTOOLS view -S -b -o $p/$sm.bam -

$JAVA_BIN -jar $PICARD/SortSam.jar I=$p/$sm.bam O=$p/$sm.sorted.bam SO=coordinate VALIDATION_STRINGENCY=SILENT TMP_DIR=. CREATE_INDEX=true
$SAMTOOLS flagstat $p/$sm.sorted.bam > $p/$sm.sorted.bam.flagstat
$SAMTOOLS depth $p/$sm.sorted.bam | perl $GENOMECOVERAGE > $p/$sm.sorted.bam.genomecoverage
$FASTQC -o QCreport -f $p/$sm.sorted.bam
$JAVA_BIN -jar $PICARD/MarkDuplicates.jar I=$p/$sm.sorted.bam O=$p/$sm.dedup.bam M=$p/$sm.metric VALIDATION_STRINGENCY=SILENT TMP_DIR=. CREATE_INDEX=true
$SAMTOOLS flagstat $p/$sm.dedup.bam > $p/$sm.dedup.bam.flagstat
$SAMTOOLS depth $p/$sm.dedup.bam | perl $GENOMECOVERAGE > $p/$sm.dedup.bam.genomecoverage
" > $p/scripts/bwa2_$sm.sh

job=`qsub -hold_jid $hold_jid -e $p/log/bwa2_$sm.log -o $p/log/bwa2_$sm.log $p/scripts/bwa2_$sm.sh`

jobs[$count]=`echo $job | awk '{print $3}'`
        for i in ${jobs[*]}; do
                hold_jid=$hold_jid${i},
        done

}

######## FUNCTION feed_bwa #########
# find all fastq.gz pairs 
# create Aligned folder
# prepare calls for bwa 
####################################

function feed_bwa {
sample_sheet=$1
run_base=$2  # /mnt/isilon/cag/ngs/hiseq/111019_SN1089_0052_AC0995ACXX


#check files
check_files


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
r1=$run_base/Data/Intensities/BaseCalls/Unaligned/Project_$SampleProject/Sample_$SampleID/${SampleID}_${Index}_L00${Lane}_R1_001
r2=$run_base/Data/Intensities/BaseCalls/Unaligned/Project_$SampleProject/Sample_$SampleID/${SampleID}_${Index}_L00${Lane}_R2_001
if [ ! -f $r1.fastq.gz ]; then
echo $r1.fastq.gz not exisit. Exit.
exit
fi

if [ ! -f $r2.fastq.gz ]; then
echo $r2.fastq.gz not exisit. Exit.
exit
fi

echo bwa $r1 $r2 $SampleID
fi
done

}


feed_bwa samplesheet.csv /mnt/isilon/cag/ngs/hiseq/111019_SN1089_0052_AC0995ACXX

