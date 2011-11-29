# GATK variant caller for whole genome NGS sequencing
# INPUT : a set of BAM files
# OUTPUT: vcf
# VERSION: Nov 23, 2011
# NOTE: GATK v1.2-59-gd7367c1 changed: 
#       1. -B: binding --> -V for variants; -input
#       2. -comp,VCF   --> -comp x.vcf ; -eval y.vcf

# step1: cleanup BAMs, produce realigned, recalibrated BAMs
# step2: call variants using UG -glm BOTH (both SNP and INDEL)
# step3: collect mapping stats
# step4: filter variants by standard hard-filtering or VQSR (both SNP and INDEL)
# step5: generate individual VCFs for each sample
# step6: Impute with BEAGLE (optional)

# to skip a step
# create a stepx.log in the current dir



######################### How to Run the Pipeline #########################
#
# modify the following parameters
## run folder, default is current dir
p=`pwd`   
bam=$p

## reference assembly
## genomes were aligned to hg19, all GATK b37 files were converted from "1" to "chr1"
ref_version=b37

##output prefix
run=SN1109_0079_2

#### Make sure the bamlist contains the actual sample names in the bam file @RG SM tag!!!!!
bamlist=(
SN1109_0079_lane2
)
#### samples list
smlist=${bamlist[*]}

### Variant Quality Score Recalibration? ###
VQSR=yes


## Default UnifiedGenotyper options;
opt="-stand_call_conf 10 -stand_emit_conf 30"

# SOLiD
#platform="-dP solid -solid_nocall_strategy PURGE_READ"
# Illumina
platform="-dP illumina "


## Please check the storage space
# make sure you have enough tmp space!!!!
# -Djava.io.tmpdir=/somewhere/gatktmp
# I use /scratch/local if I qsub

# For imputation, Beagle requires huge amount of memory
# 24gb is the default size in this script


#### 
# GATK
#export JAVA_BIN=$JAVA_HOME/bin/java
#export GATK=/scr2/caglab/gatk
#export GATK_JAR=$GATK/GenomeAnalysisTK.jar
#export JAVA_GATK="$JAVA_BIN -Djava.io.tmpdir=. -Xmx6g -jar $GATK_JAR"
#export HG19=$GATK/hg19/hg19.fa
#export PICARD=$GATK/picard
#export SAMTOOLS=$GATK/bin/samtools
#export BWA=$GATK/bin/bwa
#export FASTQC=$GATK/FastQC

temp=$p

R=$GATK/R
Rscript=$R/2.7.0/bin/Rscript
Rresources=$GATK/git/public/R

beagle="/usr/java/latest/bin/java -Djava.io.tmpdir=$temp -Xmx24g -jar $gatk/beagle/beagle.jar"



chr=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)


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
    $java -jar $picard/MarkDuplicates.jar REMOVE_DUPLICATES=true I=$bam/$i.bam O=$p/$i.dedup.bam M=$p/$i.metrics VALIDATION_STRINGENCY=SILENT  CREATE_INDEX=true ASSUME_SORTED=true
    

    # Local realignment around Indels
    # 1. Prepare Realign Interval file
    $javagatk \
        -R $ref \
        -o $p/$i.dedup.intervals \
        -I $p/$i.dedup.bam \
	-T RealignerTargetCreator

    # 2. Indel Realigner
    $javagatk  \
        -I $p/$i.dedup.bam \
        -R $ref \
        -targetIntervals $p/$i.dedup.intervals \
        -o $p/$i.dedup.indelrealigner.bam \
	-T IndelRealigner
        

    # Base quality score recalibration
    # 1.Count covariates 

	# CountCovariates requires indexed BAMs
    $samtools index $p/$i.dedup.indelrealigner.bam

    $javagatk -R $ref  \
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
    $javagatk -R $ref  \
        -I $p/$i.dedup.indelrealigner.bam \
        --out $p/$i.dedup.indelrealigner.recal.bam \
        -recalFile $p/$i.dedup.indelrealigner.recal.csv \
        $platform \
        -T TableRecalibration \

    $samtools index $p/$i.dedup.indelrealigner.recal.bam
    
    " > $p/gatkscripts/$run.$i.step1.sh

    job=`qsub -e $p/log/$run.$i.step1.log -o $p/log/$run.$i.step1.log $p/gatkscripts/$run.$i.step1.sh`
    jobs[$count]=`echo $job | awk '{print $3}'`
    echo [1 BAM cleanup] ${jobs[$count]} submitted
    let count=$count+1

done


for i in ${jobs[*]}; do
    hold_jid=$hold_jid${i},
done

echo "done" >> $p/step1.log
date >> $p/step1.log
fi # step1 done
        
             

################################## Step 2 Variant Calling #############################################

if [ ! -f $p/step2.log ]; then
	echo step2 started > $p/step2.log
	date >> $p/step2.log

	for i in ${bamlist[*]};do
    		bams=$bams" -I $p/"$i.dedup.indelrealigner.recal.bam
	done

	count=0

	for c in ${chr[*]}; do
    		echo "$javagatk -T UnifiedGenotyper -glm BOTH  -R $ref $bams -o $vcf.$c.vcf $opt  -L $c  -A AlleleBalance -A DepthOfCoverage -A FisherStrand
    		" > $p/gatkscripts/$run.step2.$c.sh
    	
		job=`qsub -hold_jid $hold_jid -e $p/log/$run.step2.$c.log -o $p/log/$run.step2.$c.log $p/gatkscripts/$run.step2.$c.sh`
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

################################### Step 3 Collect Alignement and Depth of Coverage information  ########################################################

for i in ${bamlist[*]}; do
    echo "
    # Summary report
    $java -jar $picard/CollectAlignmentSummaryMetrics.jar I=$p/$i.dedup.indelrealigner.recal.bam O=$p/$i.dedup.indelrealigner.recal.picard.summary R=$ref VALIDATION_STRINGENCY=SILENT
    $java -jar $picard/MeanQualityByCycle.jar I=$p/$i.dedup.indelrealigner.recal.bam O=$p/$i.dedup.indelrealigner.recal.meanquality \
        CHART=$p/$i.dedup.indelrealigner.recal.meanquality.pdf VALIDATION_STRINGENCY=SILENT
    $samtools flagstat $p/$i.dedup.indelrealigner.recal.bam > $p/$i.dedup.indelrealigner.recal.flagstat
    " > $p/gatkscripts/$run.$i.step3.dedup.summary.sh
	
    job=`qsub -hold_jid $hold_jid -o $p/log/$run.$i.step3.dedup.log $p/gatkscripts/$run.$i.step3.dedup.summary.sh`
    tempjob=`echo $job | awk '{print $3}'`
    echo [1 SummaryMetrics: $tempjob ] waiting for UnifiedGenotyper: $hold_jid
                       
    echo "
	$javagatk -T DepthOfCoverage -R $ref -I $p/$i.dedup.indelrealigner.recal.bam -o $p/$i.recal.doc -omitBaseOutput -ct 1 -ct 2 -ct 3 -ct 4 -ct 5 -ct 6 -ct 7 -ct 8
	" > $p/gatkscripts/$run.$i.step3.recal.doc.sh
   job=`qsub -hold_jid $hold_jid  -o $p/log/$run.$i.step3.recal.doc.log $p/gatkscripts/$run.$i.step3.recal.doc.sh`
   tempjob=`echo $job | awk '{print $3}'`
   echo [1 DepthOfCoverage: $tempjob ] waiting for UnifiedGenotyper: $hold_jid
	
   summarize original BAM statistics
    echo "
       $java -jar $picard/CollectAlignmentSummaryMetrics.jar I=$p/$i.bam O=$p/$i.picard.summary R=$ref VALIDATION_STRINGENCY=SILENT
       $java -jar $picard/MeanQualityByCycle.jar I=$p/$i.bam O=$p/$i.meanquality CHART=$p/$i.meanquality.pdf VALIDATION_STRINGENCY=SILENT
       $samtools flagstat $p/$i.bam > $p/$i.flagstat
       $javagatk -T DepthOfCoverage -R $ref -I $p/$i.bam -o $p/$i.doc -omitBaseOutput -ct 1 -ct 2 -ct 3 -ct 4 -ct 5 -ct 6 -ct 7 -ct 8
       " > $p/gatkscripts/$run.$i.step3.doc.sh
    job=`qsub -hold_jid $hold_jid -o $p/log/$run.$i.step3.doc.log $p/gatkscripts/$run.$i.step3.doc.sh`
    tempjob=`echo $job | awk '{print $3}'`
    echo [1 Summary of Original BAM: $tempjob ] waiting for UnifiedGenotyper: $hold_jid

done


if [ ! -f $p/step4.log ]; then
	echo step4 started > $p/step4.log
	date >> $p/step4.log

############################### Step 4 combine Variants, Annotation, VQSR recalibration ###############################################
	for c in ${chr[*]}; do
    		vcf_chr=$vcf_chr" -V $vcf.$c.vcf"
	done


	echo "
    		#### combine VCFs
    		$javagatk -T CombineVariants $vcf_chr -o $vcf -R $ref

    		#### produce SNP and INDEL VCFs 
    		$javagatk -T SelectVariants \
        		-R $ref \
        		-V $vcf \
			-o $snp  \
			-env   \
			-selectType SNP  

    		$javagatk -T SelectVariants \
			-R $ref \
        		-V $vcf \
			-o $indel  \
			-env  \
			-selectType INDEL 

	" > $p/gatkscripts/$run.step4.sh


	### for whole exome, we have several choices here
	### Here in this pipeline, we do GATK
	### followed by
	### 1. Normal QUAL and QD filtering
	### 2. VQSR


	## Oct 27, 2011 based on Best Practice 3
	## QD < 2.0", "MQ < 40.0", "FS > 60.0", "HaplotypeScore > 13.0", "MQRankSum < -12.5", "ReadPosRankSum < -8.0".
	echo "
		$javagatk -T VariantFiltration \
			-V $snp \
			-o $snp.temp.vcf \
			--clusterSize 3 \
			--clusterWindowSize 10 \
			--filterExpression \"QD < 2.0 || MQ < 40.0 || FS > 60.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" \
			--filterName GATKSNPStandard \
			-R $ref \

		$javagatk -T SelectVariants \
			-R $ref \
			-V $snp.temp.vcf \
			-o $snp.gatkstandard.vcf  \
			-env \
			-ef


	" >> $p/gatkscripts/$run.step4.sh

	### Indel
	### DATA_TYPE_SPECIFIC_FILTERS should be "QD < 2.0", "ReadPosRankSum < -20.0", "InbreedingCoeff < -0.8", "FS > 200.0".
	echo "
	    $javagatk -T VariantFiltration \
        	-V $indel \
        	-o $indel.temp.vcf \
        	--clusterSize 3 \
        	--clusterWindowSize 10 \
        	--filterExpression \"QD < 2.0 ||  ReadPosRankSum < -20.0 || FS > 200.0 || InbreedingCoeff < -0.8  \" \
        	--filterName GATKINDELStandard \
        	-R $ref \

	
    	$javagatk -T SelectVariants \
        	-R $ref \
        	-V $indel.temp.vcf \
        	-o $indel.gatkstandard.vcf  \
        	-env \
        	-ef


	" >> $p/gatkscripts/$run.step4.sh

	if [ "$VQSR" = "yes" ]; then
		echo "
    		###########################################################################################
    		# VQSR
    		# IN: variants as input, hapmap, 1kg, dbsnp as training
    		# OUT: cluster file
    		###########################################################################################

   		$javagatk \
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
   			-rscriptFile $p/rscript.snp.r \

   		$javagatk \
   			-T ApplyRecalibration \
   			-R $ref \
   			-input $snp \
   			--ts_filter_level 99.0 \
   			-tranchesFile $snp_tranches \
   			-recalFile $snp_recal \
   			-o $snp.vqsr.vcf

    		#indel -an InbreedingCoeff is not included
   		$javagatk \
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
   			-rscriptFile $p/rscript.indel.r \

   		$javagatk \
   			-T ApplyRecalibration \
   			-R $ref \
   			-input $indel \
   			--ts_filter_level 99.0 \
   			-tranchesFile $indel_tranches \
   			-recalFile $indel_recal \
   			-o $indel.vqsr.vcf

		" >> $p/gatkscripts/$run.step4.sh
	fi

	job=`qsub -hold_jid $hold_jid -e $p/log/$run.step4.log -o $p/log/$run.step4.log $p/gatkscripts/$run.step4.sh`
	hold_jid=`echo $job | awk '{print $3}'`
	echo [3 SelectVariants: $filter_job ] waiting for UnifiedGenotyper: $hold_jid 

	echo step4 ended >> $p/step4.log
	date >> $p/step4.log
fi #end_of_step4





######################################### step 5 #########################################
if [ ! -f $p/step5.log ]; then
echo step5 started  > $p/step5.log
date >> $p/step5.log


#### Produce individual VCFs  
for i in ${smlist[*]};do

    echo "
    	$javagatk -T SelectVariants -sn $i -R $ref -V $snp.gatkstandard.vcf -o $p/$i.gatkstandard.snp.vcf  -env -ef 
    	$javagatk -T SelectVariants -sn $i -R $ref -V $indel.gatkstandard.vcf -o $p/$i.gatkstandard.indel.vcf  -env -ef 
    	$javagatk -T VariantEval  -eval $p/$i.gatkstandard.snp.vcf -D $dbsnp  -R $ref -o $p/$i.gatkstandard.report
    "     > $p/gatkscripts/$run.step5.VariantEval.$i.sh

   if [ "$VQSR" = "yes" ]; then
        echo "
        $javagatk -T SelectVariants -sn $i -R $ref -V $snp.vqsr.vcf -o $p/$i.gatkvqsr.snp.vcf  -env -ef 
       	$javagatk -T VariantEval  -eval $p/$i.gatkvqsr.snp.vcf -comp $p/$i.gatkstandard.snp.vcf -D $dbsnp  -R $ref -o $p/$i.gatkvqsr.report
    	"     >> $p/gatkscripts/$run.step5.VariantEval.$i.sh
		
	fi
	job=`qsub -hold_jid $hold_jid -e $p/log/$run.step5.VariantEval.$i.log -o $p/log/$run.step5.VariantEval.$i.log $p/gatkscripts/$run.step5.VariantEval.$i.sh`
	tempjob=`echo $job | awk '{print $3}'`
    echo [3 VariantEval: $tempjob]  waiting for : $filter_job 

    

done

echo step5 ended >> $p/step5.log
date >> $p/step5.log
fi # end of step5


if [ ! -f $p/step6.log ]; then
echo step6 started  > $p/step6.log 
date >> $p/step6.log

###################################### step 5: Impute with Beagle #############################################

echo "
    date
    $javagatk   -R $ref -T ProduceBeagleInput  -V $snp.gatkstandard.vcf -o $p/$run.beagle

    $beagle like=$p/$run.beagle out=

    gunzip -c $p/$run.beagle.dose.gz > $p/$run.beagle.dose
    gunzip -c $p/$run.beagle.gprobs.gz > $p/$run.beagle.gprobs
    gunzip -c $p/$run.beagle.phased.gz > $p/$run.beagle.phased

    $javagatk -T BeagleOutputToVCF \
        -R $ref \
		-V $recal.filtered.vcf \
        -beagleR2 $p/$run.beagle.r2 \
        -beaglePhased $p/$run.beagle.phased \
        -beagleProbs $p/$run.beagle.gprobs \
        --out $p/$run.beagle.vcf
    
   date
" > $p/gatkscripts/$run.step5.beagle.sh

job=`qsub -hold_jid $vqsr_job -e $p/log/$run.step5.beagle.log -o $p/log/$run.step5.beagle.log  $p/gatkscripts/$run.step5.beagle.sh` 
beagle_job=`echo $job | awk '{print $3}'`
echo [4 Beagle imputation: $beagle_job]  waiting for VQSR: $vqsr_job 

#### Produce individual VCFs from Beagle output
#### Compare imputed and non-imputed results against each other and dbSNP
for i in ${smlist[*]};do
    echo "
    $javagatk -T SelectVariants \
	-sn $i \
	-R $ref \
	-V $p/$run.beagle.vcf \
	-o $p/$i.beagle.vcf  \
	-env \
	-ef 
  
    $javagatk -T VariantEval  \
	-eval $p/$i.recal.filtered.vcf \
	-comp $p/$i.beagle.vcf \
	-D $dbsnp  \
	-R $ref \
	-o $p/$i.normal_vs_beagle.gatkreport 
    
    $javagatk -T VariantEval  \
	-eval $p/$i.beagle.vcf \
	-D $dbsnp  \
	-R $ref \
	-o $p/$i.beagle.gatkreport
    "     > $p/gatkscripts/$run.$i.beagle.VariantEval.sh

    job=`qsub -e $p/log/$run.$i.beagle.VariantEval.log -o $p/log/$run.$i.beagle.VariantEval.log -hold_jid $beagle_job $p/gatkscripts/$run.$i.beagle.VariantEval.sh`
	tempjob=`echo $job | awk '{print $3}'`
	echo [4 VariantEval: $tempjob] waiting for Beagle Imputation: $beagle_job
done


echo step6 ended >> $p/step6.log
date >> $p/step6.log
fi # end of step6

