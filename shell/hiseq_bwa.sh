
source setup.sh

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
#  exit
fi
if [ ! -f $r2.fastq.gz ]; then
  echo $r2.fastq.gz not exist!
#  exit
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

#align
timestamp=`date +%s`
echo "
$BWA aln  -t 4 $HG19 $r1.fastq.gz > $r1.sai
$BWA aln  -t 4 $HG19 $r2.fastq.gz > $r2.sai
" > $p/scripts/$sm.aln.sh

echo $p/scripts/$sm.aln.sh
cat $p/scripts/$sm.aln.sh
echo
echo

#job=`qsub -pe smp 4 -e $p/log/$sm.aln.err.log -o $p/log/$sm.aln.log $p/scripts/$sm.aln.sh`
#jobs=`echo $job | awk '{print $3}'`

#hold_jid=$jobs

echo "
$BWA sampe -r \"@RG\\tID:$fid-$lane-$index-$sm\\tSM:$sm\" $HG19 $r1.sai $r2.sai $r1.fastq.gz $r2.fastq.gz  | $SAMTOOLS view -S -b -o $p/$sm.bam -
$JAVA_BIN -jar $PICARD/SortSam.jar I=$p/$sm.bam O=$p/$sm.sorted.bam SO=coordinate VALIDATION_STRINGENCY=SILENT TMP_DIR=$p/temp CREATE_INDEX=true
$SAMTOOLS flagstat $p/$sm.sorted.bam > $p/$sm.sorted.bam.flagstat
$SAMTOOLS depth $p/$sm.sorted.bam | perl $GENOMECOVERAGE > $p/$sm.sorted.bam.genomecoverage
$FASTQC -o QCreport -f $p/$sm.sorted.bam
$JAVA_BIN -jar $PICARD/MarkDuplicates.jar I=$p/$sm.sorted.bam O=$p/$sm.dedup.bam M=$p/$sm.metric VALIDATION_STRINGENCY=SILENT TMP_DIR=$p/temp CREATE_INDEX=true REMOVE_DUPLICATES=true
$SAMTOOLS flagstat $p/$sm.dedup.bam > $p/$sm.dedup.bam.flagstat
$SAMTOOLS depth $p/$sm.dedup.bam | perl $GENOMECOVERAGE > $p/$sm.dedup.bam.genomecoverage
" > $p/scripts/$sm.sampe.sh

echo $p/scripts/$sm.sampe.sh
cat $p/scripts/$sm.sampe.sh

#job=`qsub -hold_jid $hold_jid -e $p/log/$sm.sampe.err.log -o $p/log/$sm.sample.log $p/scripts/$sm.sampe.sh`

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

