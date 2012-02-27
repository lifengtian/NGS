
function gatk_pipeline {

# GATK variant caller for  NGS sequencing
# INPUT : a set of BAM files
# OUTPUT: vcf
# VERSION: Nov 23, 2011
# NOTE: GATK v1.2-59-gd7367c1 changed: 
#       1. -B: binding --> -V for variants; -input
#       2. -comp,VCF   --> -comp x.vcf ; -eval y.vcf

# step1: cleanup BAMs, produce realigned, recalibrated BAMs
# step2: call variants using UG -glm BOTH (both SNP and INDEL)
# step3: filter variants by standard hard-filtering or VQSR (both SNP and INDEL)
# step4: generate individual VCFs for each sample

# to skip a step
# create a stepx.log in the current dir
p=$1

######################### How to Run the Pipeline #########################
#
# modify the following parameters
## run folder, default is current dir
bam=$p/BAM

## reference assembly
## genomes were aligned to hg19, all GATK b37 files were converted from "1" to "chr1"
ref_version=b37

#### Make sure the bamlist contains the actual sample names in the bam file @RG SM tag!!!!!
#### samples list
smlist=${bamlist[*]}

### Variant Quality Score Recalibration? ###


## Default UnifiedGenotyper options;
opt="-stand_call_conf 10 -stand_emit_conf 30"

# SOLiD
#platform="-dP solid -solid_nocall_strategy PURGE_READ"
# Illumina
#platform="-dP illumina "



## Please check the storage space
# make sure you have enough tmp space!!!!
# -Djava.io.tmpdir=/somewhere/gatktmp
# I use /scratch/local if I qsub

# For imputation, Beagle requires huge amount of memory
# 24gb is the default size in this script


#### 
# GATK
#export JAVA_BIN=$JAVA_HOME/bin/JAVA
#export GATK=/scr2/caglab/gatk
#export GATK_JAR=$GATK/GenomeAnalysisTK.jar
#export JAVA_GATK="$JAVA_BIN -Djava.io.tmpdir=. -Xmx6g -jar $GATK_JAR"
#export HG19=$GATK/hg19/hg19.fa
#export PICARD=$GATK/picard
#export SAMTOOLS=$GATK/bin/SAMTOOLS
#export BWA=$GATK/bin/bwa
#export FASTQC=$GATK/FastQC

temp=$p

R=$GATK/R
Rscript=$R/2.7.0/bin/Rscript
Rresources=$GATK/git/public/R

beagle="/usr/JAVA/latest/bin/JAVA -Djava.io.tmpdir=$temp -Xmx24g -jar $gatk/beagle/beagle.jar"



chr=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM)


if [ "$ref_version" = "hg18" ]; then
#### Global variables for gatk hg18
ref=$gatk/hg18/human_hg18_bioscope.fasta
hapmap=$gatk/hg18/hapmap_3.3.hg18.sites.vcf
kg=$gatk/hg18/1000G_omni2.5.hg18.sites.vcf
dbsnp=$gatk/hg18/dbsnp_132.hg18.vcf
kg_indels=1000G_indels_for_realignment.hg18.vcf

elif [ "$ref_version" = "b37" ]; then

ref=$HG19
hapmap=$GATK/b37/hapmap_3.3.b37.sites.vcf
kg=$GATK/b37/1000G_omni2.5.b37.sites.vcf
dbsnp=$GATK/b37/dbsnp_132.b37.vcf
kg_indels=$GATK/b37/1000G_indels_for_realignment.b37.vcf
Mills_Devine=$GATK/b37/Mills_Devine_2hit.indels.b37.sites.vcf

else
	echo "$ref_version has to be hg18 or b37"
	exit
fi




#### pipeline output files
VCFdir=$p/VCF
vcf=$VCFdir/$run.vcf
snp=$VCFdir/$run.SNP.vcf
indel=$VCFdir/$run.INDEL.vcf
clusterFile=$VCFdir/$run.clusters
tabl=$$VCFdir/$run.tabl
VRR=$VCFdir/$run.VRR.report
snp_recal=$VCFdir/$run.snp.recal
indel_recal=$VCFdir/$run.indel.recal
snp_tranches=$VCFdir/$run.snp.tranches
indel_tranches=$VCFdir/$run.indel.tranches

#### create folders
mkdir -p $p/gatkscripts


count=0                                    # temp variable to count jobs submitted



#### check required files
check_list=($ref $hapmap $kg $dbsnp \
$SAMTOOLS \
$GATK_JAR \
$GATK/beagle/beagle.jar \
$PICARD/MarkDuplicates.jar \
$PICARD/CollectAlignmentSummaryMetrics.jar \
$PICARD/MeanQualityByCycle.jar \
$Rscript \
$Mills_Devine \
$HG19 )

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

#### DON'T forget the ','!
hold_jid=100000000,

date
echo $bamlist
#### Step 1. Cleanup BAM files

if [ ! -f $p/step1.log ]; then
echo step1 started
date > $p/step1.log

for i in ${bamlist[*]}; do

    echo "
    ###############################################################################
    #cleanup bams (index, rmdup, index) for $i
    ###############################################################################

    
    # remove duplicates
    #$JAVA_BIN -jar $PICARD/MarkDuplicates.jar REMOVE_DUPLICATES=true I=$bam/$i.bam O=$bam/$i.dedup.bam M=$bam/$i.metrics VALIDATION_STRINGENCY=SILENT  CREATE_INDEX=true ASSUME_SORTED=true
    

    # Local realignment around Indels
    # 1. Prepare Realign Interval file
    $JAVA_GATK \
        -R $ref \
        -o $bam/$i.dedup.intervals \
        -I $bam/$i.dedup.bam \
	-T RealignerTargetCreator

    # 2. Indel Realigner
    $JAVA_GATK  \
        -I $bam/$i.dedup.bam \
        -R $ref \
        -targetIntervals $bam/$i.dedup.intervals \
        -o $bam/$i.dedup.indelrealigner.bam \
	-T IndelRealigner
        

    # Base quality score recalibration
    # 1.Count covariates 

	# CountCovariates requires indexed BAMs
    $SAMTOOLS index $p/$i.dedup.indelrealigner.bam

    $JAVA_GATK -R $ref  \
        -knownSites $dbsnp \
        -I $bam/$i.dedup.indelrealigner.bam \
        --covariate ReadGroupCovariate \
        --covariate QualityScoreCovariate \
        --covariate CycleCovariate \
        --covariate DinucCovariate \
        -recalFile $bam/$i.dedup.indelrealigner.recal.csv \
        $platform \
        -T CountCovariates \
    

    # 2. Table Recalibration 
    $JAVA_GATK -R $ref  \
        -I $bam/$i.dedup.indelrealigner.bam \
        --out $bam/$i.dedup.indelrealigner.recal.bam \
        -recalFile $bam/$i.dedup.indelrealigner.recal.csv \
        $platform \
        --default_read_group defaultRG \
	-T TableRecalibration \

    $SAMTOOLS index $bam/$i.dedup.indelrealigner.recal.bam
    
    " > $p/gatkscripts/$run.$i.step1.sh

    job=`qsub $queue -pe smp 1 -hold_jid $hold_jid  -e $p/logs/$run.$i.step1.log -o $p/logs/$run.$i.step1.log $p/gatkscripts/$run.$i.step1.sh`
    sleep 5  
    jobs[$count]=`echo $job | awk '{print $3}'`
    echo [1 BAM cleanup] ${jobs[$count]} submitted
    let count=$count+1

done


for i in ${jobs[*]}; do
    hold_jid=$hold_jid,${i}
done

echo "done" >> $p/step1.log
date >> $p/step1.log
fi # step1 done
        
             

################################## Step 2 Variant Calling #############################################

if [ ! -f $p/step2.log ]; then
	echo step2 started > $p/step2.log
	date >> $p/step2.log

	for i in ${bamlist[*]};do
    		bams=$bams" -I $bam/"$i.dedup.indelrealigner.recal.bam
	done

	count=0

	for c in ${chr[*]}; do
    		echo "$JAVA_GATK -T UnifiedGenotyper -glm BOTH  -R $ref $bams -o $vcf.$c.vcf $opt  -L $c  -A AlleleBalance -A DepthOfCoverage -A FisherStrand
    		" > $p/gatkscripts/$run.step2.$c.sh
    	
	job=`qsub $queue -pe smp 1 -hold_jid $hold_jid -e $p/logs/$run.step2.$c.log -o $p/logs/$run.step2.$c.log $p/gatkscripts/$run.step2.$c.sh`
	sleep 2    	
	ug_jobs[$count]=`echo $job | awk '{print $3}'`
    		echo UnifiedGenotyper job for $c: ${ug_jobs[$count]} submitted
    		let count=count+1
	done


	for i in ${ug_jobs[*]}; do
		hold_jid=$hold_jid${i},
	done

	echo step2 ended >> $p/step2.log
	date >> $p/step2.log

fi #end of step2




if [ ! -f $p/step3.log ]; then
	echo step3 started > $p/step3.log
	date >> $p/step3.log

############################### Step 4 combine Variants, Annotation, VQSR recalibration ###############################################
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

	" > $p/gatkscripts/$run.step3.sh


	### for whole exome, we have several choices here
	### Here in this pipeline, we do GATK
	### followed by
	### 1. Normal QUAL and QD filtering
	### 2. VQSR


	## Oct 27, 2011 based on Best Practice 3
	## QD < 2.0", "MQ < 40.0", "FS > 60.0", "HaplotypeScore > 13.0", "MQRankSum < -12.5", "ReadPosRankSum < -8.0".
	echo "
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


	" >> $p/gatkscripts/$run.step3.sh

	### Indel
	### DATA_TYPE_SPECIFIC_FILTERS should be "QD < 2.0", "ReadPosRankSum < -20.0", "InbreedingCoeff < -0.8", "FS > 200.0".
	echo "
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


	" >> $p/gatkscripts/$run.step3.sh

	if [ "$VQSR" = "yes" ]; then
		echo "
    		###########################################################################################
    		# VQSR
    		# IN: variants as input, hapmap, 1kg, dbsnp as training
    		# OUT: cluster file
    		###########################################################################################

   		$JAVA_GATK \
   			-T VariantRecalibrator \
   			-R $ref \
   			-input $snp \
   			-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
   			-resource:omni,known=false,training=true,truth=false,prior=12.0 $kg \
   			-resource:dbsnp,known=true,training=false,truth=false,prior=8.0 $dbsnp \
   			-an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an MQ -an FS -an DP \
   			-mode SNP \
   			-recalFile $snp_recal \
   			-tranchesFile $snp_tranches \
   			-Rscript $Rscript \
   			-resources $git/public/R \
   			-rscriptFile $p/VCF/rscript.snp.r \

   		$JAVA_GATK \
   			-T ApplyRecalibration \
   			-R $ref \
   			-input $snp \
   			--ts_filter_level 99.0 \
   			-tranchesFile $snp_tranches \
   			-recalFile $snp_recal \
   			-o $snp.vqsr.vcf

    		#indel -an InbreedingCoeff is not included
   		$JAVA_GATK \
   			-T VariantRecalibrator \
   			-R $ref \
   			-input $indel \
   			-resource:mills,VCF,known=true,training=true,truth=true,prior=12.0 $Mills_Devine \
   			-an QD -an FS -an HaplotypeScore -an ReadPosRankSum  \
   			-mode INDEL \
   			-recalFile $indel_recal \
   			-tranchesFile $indel_tranches \
   			-Rscript $Rscript \
   			-resources $Rresources \
   			-rscriptFile $p/VCF/rscript.indel.r \

   		$JAVA_GATK \
   			-T ApplyRecalibration \
   			-R $ref \
   			-input $indel \
   			--ts_filter_level 99.0 \
   			-tranchesFile $indel_tranches \
   			-recalFile $indel_recal \
   			-o $indel.vqsr.vcf

		" >> $p/gatkscripts/$run.step3.sh
	fi

	job=`qsub $queue -pe smp 1 -hold_jid $hold_jid -e $p/logs/$run.step3.log -o $p/logs/$run.step3.log $p/gatkscripts/$run.step3.sh`
	filter_jid=`echo $job | awk '{print $3}'`
	echo [3 SelectVariants: $filter_job ] waiting for UnifiedGenotyper: $hold_jid 
    hold_jid=$filter_jid

	echo step3 ended >> $p/step3.log
	date >> $p/step3.log
fi #end_of_step4





######################################### step 4 #########################################
if [ ! -f $p/step4.log ]; then
echo step4 started  > $p/step4.log
date >> $p/step4.log


#### Produce individual VCFs  
for i in ${smlist[*]};do

    echo "
    	$JAVA_GATK -T SelectVariants -sn $i -R $ref -V $snp.gatkstandard.vcf -o $p/VCF/$i.gatkstandard.snp.vcf  -env -ef 
    	$JAVA_GATK -T SelectVariants -sn $i -R $ref -V $indel.gatkstandard.vcf -o $p/VCF/$i.gatkstandard.indel.vcf  -env -ef 
    	$JAVA_GATK -T VariantEval  -eval $p/$i.gatkstandard.snp.vcf -D $dbsnp  -R $ref -o $p/VCF/$i.gatkstandard.report
    "     > $p/gatkscripts/$run.step4.VariantEval.$i.sh

   if [ "$VQSR" = "yes" ]; then
        echo "
        $JAVA_GATK -T SelectVariants -sn $i -R $ref -V $snp.vqsr.vcf -o $p/VCF/$i.gatkvqsr.snp.vcf  -env -ef 
       	$JAVA_GATK -T VariantEval  -eval $p/VCF/$i.gatkvqsr.snp.vcf -comp $p/VCF/$i.gatkstandard.snp.vcf -D $dbsnp  -R $ref -o $p/VCF/$i.gatkvqsr.report
    	"     >> $p/gatkscripts/$run.step4.VariantEval.$i.sh
		
	fi
	job=`qsub $queue -pe smp 1 -hold_jid $hold_jid -e $p/logs/$run.step4.VariantEval.$i.log -o $p/logs/$run.step4.VariantEval.$i.log $p/gatkscripts/$run.step4.VariantEval.$i.sh`
	tempjob=`echo $job | awk '{print $3}'`
    echo [3 VariantEval: $tempjob]  waiting for : $hold_job 

    

done

echo step4 ended >> $p/step4.log
date >> $p/step4.log
fi # end of step4 

}



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


