# calculate target coverage 

source setup.sh

Target=$GATK/bed/SureSelect50mbclean
Target_bed=$TARGET.bed

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
             $PICARD/MarkDuplicates.jar \
	     $Target_bed \
	     $HISEQ/NGS/perl/CalcTargetCoverage.pl \	
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
function calc_coverage {
p=$1
sm=$2
fid=$3
lane=$4
index=$5

if [ ! -f $p/$sm-$lane-$index.dedup.bam ]; then
  echo $p/$sm-$lane-$index.dedup.bam not exist!
  exit
fi

if [ ! -f $p/$sm-$lane-$index.sorted.bam ]; then
	echo $p/$sm-$lane-$index.sorted.bam not exist!
	exit
fi

echo "
$JAVA_BIN -jar $PICARD/CalculateHsMetrics.jar I=$p/$sm-$lane-$index.dedup.bam O=$p/$sm-$lane-$index.dedup.bam.target_coverage BI=$cov/SureSelect50mbclean.picard TI=$cov/SureSelect50mbclean.picard TMP_DIR=/scratch/local
mkdir $p/QCreport
$SAMTOOLS depth -b $Target_bed $p/$sm-$lane-$index.sorted.bam | perl $GENOMECOVERAGE > $p/$sm-$lane-$index.sorted.bam.target_depth
$SAMTOOLS depth -b $Target_bed $p/$sm-$lane-$index.dedup.bam | perl $GENOMECOVERAGE > $p/$sm-$lane-$index.dedup.bam.target_depth
$HISEQ/NGS/perl/CalcTargetCoverage.pl $p/$sm-$lane-$index.dedup.bam 
/mnt/isilon/cag/ngs/hiseq/gatk/FastQC/fastqc -o $p/QCreport $p/$sm-$lane-$index.dedup.bam
" > $p/scripts/$sm-$lane-$index.target_coverage.sh

#job=`qsub $queue -V -e $p/log/$sm-$lane-$index.target_coverage.err.log -o $p/log/$sm-$lane-$index.target_coverage.log $p/scripts/$sm-$lane-$index.target_coverage.sh`

}

######## FUNCTION feed_bwa #########
# find all fastq.gz pairs 
# create Aligned folder
# prepare calls for bwa 
####################################

function feed_bwa {
flow_cell_id=$1
sample_sheet=$HISEQ_ANALYSIS/SampleSheet/$flow_cell_id.csv
fastq=$HISEQ_ANALYSIS/FASTQ/$flow_cell_id
bam=$HISEQ_ANALYSIS/BAM/$flow_cell_id



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

calc_coverage $bam/Project_$SampleProject/Sample_$SampleID $SampleID $FCID $Lane $Index
#sleep 10
fi # end_of_count -gt 1

done # end_of_while read sample_sheet

}


#feed_bwa FlowCellID
feed_bwa    $1

