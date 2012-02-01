#!/bin/bash

#############################################
# Name:   hiseq_demultiplex_bwa.sh
# Input:  FlowcellID,eg., D0EF0ACXX
# Output: $HISEQ/BAM/FlowCellID/*.dedup.indelrealigner.recal.bam
# version:1.0
# 
# Note:   This script only works on resolute
#         
##############################################

source $HISEQ/NGS/shell/setup.sh
source $HISEQ/NGS/shell/function.sh

### GATK 
hapmap=$GATK/b37/hapmap_3.3.b37.sites.vcf
kg=$GATK/b37/1000G_omni2.5.b37.sites.vcf
dbsnp=$GATK/b37/dbsnp_132.b37.vcf
kg_indels=$GATK/b37/1000G_indels_for_realignment.b37.vcf
Mills_Devine=$GATK/b37/Mills_Devine_2hit.indels.b37.sites.vcf

hold_demul_jid='1000000,'
hold_bwa_jid='1000000,'
hold_sampe_jid='1000000,'
hold_gatk_jid='1000000,'

platform=' -dP illumina '
ref=$HG19

### Target
Target=$GATK/bed/SureSelect50mbclean
Target_bed=$Target.bed
Target_picard=$Target.picard


####### FUNCTION demultiplex #############
# produce zipped FASTQ using CASAVA 1.8.2
##########################################
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

# exit if demultiplexed.
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
job=`qsub $queue -pe smp 4 -e $HISEQ/logs/$FlowCellID.log -o $HISEQ/logs/$FlowCellID.log $HISEQ/FASTQ/scripts/$FlowCellID.sh`
jobs=`echo $job | awk '{print $3}'`
hold_demul_jid=$jobs

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
desc=$8
sequencer=$9

#i="$sm-$lane-$index" #for gatk cleaning-up
i=$sm

mkdir -p $p/log
mkdir -p $p/scripts
mkdir -p $p/temp
mkdir -p $p/QCreport

#align
timestamp=`date +%s`
echo "
$BWA aln -n 3 -t 4 -q 20 $HG19 $r1.fastq.gz > $r1.sai
$BWA aln -n 3 -t 4 -q 20 $HG19 $r2.fastq.gz > $r2.sai
" > $p/scripts/$i.aln.sh

job=`qsub $queue -pe smp 4 -hold_jid $hold_demul_jid -e $p/log/$i.aln.err.log -o $p/log/$i.aln.log $p/scripts/$i.aln.sh`
jobs=`echo $job | awk '{print $3}'`
echo $p $sm bwa aln jobs: $jobs 

hold_bwa_jid=$jobs

# enforce correct RG
id=$fid-$lane-$index-$sm
lb=2x101
pl=illumina
pu=$sequencer-$fid
cn=CAG

echo "
$BWA sampe -r \"@RG\\tID:$fid-$lane-$index-$sm\\tSM:$sm\" $HG19 $r1.sai $r2.sai $r1.fastq.gz $r2.fastq.gz  | $SAMTOOLS view -S -b -o $p/$i.bam -
$JAVA_BIN -jar $PICARD/SortSam.jar I=$p/$i.bam O=$p/$i.sorted.temp.bam SO=coordinate VALIDATION_STRINGENCY=SILENT TMP_DIR=$p/temp CREATE_INDEX=true
$JAVA_BIN -jar $PICARD/AddOrReplaceReadGroups.jar   I=$p/$i.sorted.temp.bam O=$p/$i.sorted.bam SO=coordinate RGID=$id  RGLB=$lb RGPL=$pl RGPU=$pu RGCN=$cn RGSM=$sm CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT TMP_DIR=$p/temp
$SAMTOOLS flagstat $p/$i.sorted.bam > $p/$i.sorted.bam.flagstat
$SAMTOOLS depth $p/$i.sorted.bam | perl $GENOMECOVERAGE > $p/$i.sorted.bam.genomecoverage
$FASTQC/fastqc -o $p/QCreport  $p/$i.sorted.bam

$JAVA_BIN -jar $PICARD/MarkDuplicates.jar I=$p/$i.sorted.bam O=$p/$i.dedup.bam M=$p/$i.metric VALIDATION_STRINGENCY=SILENT TMP_DIR=$p/temp CREATE_INDEX=true REMOVE_DUPLICATES=true
$SAMTOOLS flagstat $p/$i.dedup.bam > $p/$i.dedup.bam.flagstat
$SAMTOOLS depth $p/$i.dedup.bam | perl $GENOMECOVERAGE > $p/$i.dedup.bam.genomecoverage
" > $p/scripts/$i.sampe.sh


if [ "$desc" = "WES" ]; then
echo "
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
    #$SAMTOOLS index $p/$i.dedup.indelrealigner.bam

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

    #$SAMTOOLS index $p/$i.dedup.indelrealigner.recal.bam


	    $SAMTOOLS mpileup -uf $ref $p/$i.dedup.indelrealigner.recal.bam | $BCFTOOLS view -bvcg - > $p/$i.bcf  
	    $BCFTOOLS view $p/$i.bcf   > $p/$i.samtools.vcf

        $JAVA_BIN -jar $PICARD/CalculateHsMetrics.jar I=$p/$i.dedup.bam O=$p/$i.dedup.bam.target_coverage BI=$Target_picard TI=$Target_picard TMP_DIR=$p/temp VALIDATION_STRINGENCY=SILENT
        $SAMTOOLS depth -b $Target_bed $p/$i.sorted.bam | perl $GENOMECOVERAGE > $p/$i.sorted.bam.target_depth
        $SAMTOOLS depth -b $Target_bed $p/$i.dedup.bam | perl $GENOMECOVERAGE > $p/$i.dedup.bam.target_depth

" >> $p/scripts/$i.sampe.sh




fi

job=`qsub $queue -pe smp 1 -hold_jid $hold_bwa_jid -e $p/log/$i.sampe.err.log -o $p/log/$i.sample.log $p/scripts/$i.sampe.sh`
jobs=`echo $job | awk '{print $3}'`
hold_sampe_jid=$jobs

echo $p $i bwa_sampe jobs: $jobs  waiting on $hold_bwa_jid
}


function gatk {
######################### GATK UnifiedGenotyper ########################		
p=$1  #path
i=$2  #sample name

##output prefix
run=$i

#### Make sure the bamlist contains the actual sample names in the bam file @RG SM tag!!!!!
bamlist=(
$i
)
#### samples list
smlist=${bamlist[*]}


## Default UnifiedGenotyper options;
opt="-stand_call_conf 10 -stand_emit_conf 30"
# SOLiD
#platform="-dP solid -solid_nocall_strategy PURGE_READ"

chr=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM)


#### pipeline output files
vcf=$p/$run.vcf
snp=$p/$run.SNP.vcf
indel=$p/$run.INDEL.vcf
clusterFile=$p/$run.clusters
tabl=$p/$run.tabl
VRR=$p/$run.VRR.report
snp_recal=$p/$run.snp.recal
indel_recal=$p/$run.indel.recal
snp_tranches=$p/$run.snp.tranches
indel_tranches=$p/$run.indel.tranches

#### create folders
mkdir -p $p/log
mkdir -p $p/gatkscripts


count=0                                    # temp variable to count jobs submitted

for c in ${chr[*]}; do
                echo "$JAVA_GATK -T UnifiedGenotyper -glm BOTH  -R $ref -I $p/$i.dedup.indelrealigner.recal.bam -o $vcf.$c.vcf $opt  -L $c  -A AlleleBalance -A DepthOfCoverage -A FisherStrand
                " > $p/gatkscripts/$run.GATK1.$c.sh
	#echo qsub $queue -hold_jid $hold_jid -e $p/log/$run.GATK1.$c.log -o $p/log/$run.GATK1.$c.log $p/gatkscripts/$run.GATK1.$c.sh
        job=`qsub $queue -pe smp 1 -hold_jid $hold_sampe_jid -e $p/log/$run.GATK1.$c.log -o $p/log/$run.GATK1.$c.log $p/gatkscripts/$run.GATK1.$c.sh`
        ug_jobs[$count]=`echo $job | awk '{print $3}'`
        echo UnifiedGenotyper job for $c: ${ug_jobs[$count]} submitted
        let count=count+1
done

hold_gatk_jid=''

        for ji in ${ug_jobs[*]}; do
                hold_gatk_jid=$hold_gatk_jid${ji},
        done

for c in ${chr[*]}; do
                vcf_chr=$vcf_chr" -V $vcf.$c.vcf"
        done


        echo "
                #### combine VCFs
                $JAVA_GATK -T CombineVariants $vcf_chr -o $vcf -R $ref

                #### produce SNP and INDEL VCFs 
                $JAVA_GATK -T SelectVariants \
                        -R $ref \
                        -V $vcf \
                        -o $snp  \
                        -env   \
                        -selectType SNP  

                $JAVA_GATK -T SelectVariants \
                        -R $ref \
                        -V $vcf \
                        -o $indel  \
                        -env  \
                        -selectType INDEL 


        ### for whole exome, we have several choices here
        ### Here in this pipeline, we do GATK
        ### followed by
        ### 1. Normal QUAL and QD filtering
        ### 2. VQSR


        ## Oct 27, 2011 based on Best Practice 3
        ## QD < 2.0, MQ < 40.0, FS > 60.0, HaplotypeScore > 13.0, MQRankSum < -12.5, ReadPosRankSum < -8.0.
                $JAVA_GATK -T VariantFiltration \
                        -V $snp \
                        -o $snp.temp.vcf \
                        --clusterSize 3 \
                        --clusterWindowSize 10 \
                        --filterExpression \"QD < 2.0 || MQ < 40.0 || FS > 60.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" \
                        --filterName GATKSNPStandard \
                        -R $ref \

                $JAVA_GATK -T SelectVariants \
                        -R $ref \
                        -V $snp.temp.vcf \
                        -o $snp.gatkstandard.vcf  \
                        -env \
                        -ef
	
	 $JAVA_GATK -T VariantFiltration \
                -V $indel \
                -o $indel.temp.vcf \
                --clusterSize 3 \
                --clusterWindowSize 10 \
                --filterExpression \"QD < 2.0 ||  ReadPosRankSum < -20.0 || FS > 200.0 || InbreedingCoeff < -0.8  \" \
                --filterName GATKINDELStandard \
                -R $ref \

        
        $JAVA_GATK -T SelectVariants \
                -R $ref \
                -V $indel.temp.vcf \
                -o $indel.gatkstandard.vcf  \
                -env \
                -ef

	$JAVA_GATK -T SelectVariants -sn $i -R $ref -V $snp.gatkstandard.vcf -o $p/$i.gatkstandard.snp.vcf  -env -ef 
        $JAVA_GATK -T SelectVariants -sn $i -R $ref -V $indel.gatkstandard.vcf -o $p/$i.gatkstandard.indel.vcf  -env -ef 
        $JAVA_GATK -T VariantEval  -eval $p/$i.gatkstandard.snp.vcf -D $dbsnp  -R $ref -o $p/$i.gatkstandard.report
		" >  $p/gatkscripts/$run.GATK2.sh 

	### Annovar
	annovar=$GATK/annovar
annovar_db=$annovar/humandb/
ver=" --genetype knowngene --ver1000g 1000g2011may --buildver hg19 "

mkdir -p $p/annovar

echo "
$annovar/convert2annovar.pl -format vcf4 $p/$i.gatkstandard.snp.vcf > $p/annovar/$i.snp.av
$annovar/auto_annovar.pl $ver -step 1,4,7-9 -model recessive --verdbsnp 132 $p/annovar/$i.snp.av $annovar_db 
$annovar/summarize_annovar.pl  --verdbsnp 132 $ver --outfile $p/annovar/$i.snp.sum $p/annovar/$i.snp.av $annovar_db 

$annovar/convert2annovar.pl -format vcf4 $p/$i.gatkstandard.indel.vcf > $p/annovar/$i.indel.av
$annovar/auto_annovar.pl $ver -step 1,4,7-9 -model recessive --verdbsnp 132 $p/annovar/$i.indel.av $annovar_db 
$annovar/summarize_annovar.pl  --verdbsnp 132 $ver --outfile $p/annovar/$i.indelsum $p/annovar/$i.indel.av $annovar_db

" >> $p/gatkscripts/$run.GATK2.sh



	job=`qsub $queue -pe smp 1 -hold_jid $hold_gatk_jid -e $p/log/$run.GATK2.log -o $p/log/$run.GATK2.log $p/gatkscripts/$run.GATK2.sh`
}	### end_of_GATK



######## FUNCTION feed_bwa #########
# find all fastq.gz pairs 
# create BAM folder
# prepare calls for bwa 
####################################

function feed_bwa {
flow_cell_id=$1
sequencer=$2
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
#	exit
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
#if [ ! -f $r1.fastq.gz ]; then
#	echo $r1.fastq.gz not exisit. Exit.
##exit
#fi

#if [ ! -f $r2.fastq.gz ]; then
#	echo $r2.fastq.gz not exisit. Exit.
##exit
#fi

# mkdir BAM folder
    mkdir -p $bam/Project_$SampleProject/Sample_$SampleID
    newr1=$bam/Project_$SampleProject/Sample_$SampleID/${SampleID}_${Index}_L00${Lane}_R1_001
    newr2=$bam/Project_$SampleProject/Sample_$SampleID/${SampleID}_${Index}_L00${Lane}_R2_001
    ln -s $r1.fastq.gz $newr1.fastq.gz
    ln -s $r2.fastq.gz $newr2.fastq.gz

bwa $bam/Project_$SampleProject/Sample_$SampleID $newr1 $newr2 $SampleID $FCID $Lane $Index $description $sequencer
#gatk $bam/Project_$SampleProject/Sample_$SampleID $SampleID
#sleep 10
fi # end_of_count -gt 1

done # end_of_while read sample_sheet

} # end_of_function_feed_bwa

demultiplex $1
#feed_bwa FlowCellID Sequencer
feed_bwa    $1 $2

