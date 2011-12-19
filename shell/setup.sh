 
 HISEQ_ANALYSIS=/mnt/isilon/cag/ngs/hiseq
 GENOMECOVERAGE=$GATK/NGS/perl/genomecoverage.pl

 JAVA_BIN=$JAVA_HOME/bin/java
 GATK=/mnt/isilon/cag/ngs/hiseq/gatk
 GATK_JAR=$GATK/GenomeAnalysisTK.jar
 JAVA_GATK="$JAVA_BIN -Djava.io.tmpdir=. -Xmx12g -jar $GATK_JAR"
 HG19=$GATK/hg19/hg19.fa
 PICARD=$GATK/picard
 SAMTOOLS=$GATK/bin/samtools
 BWA=$GATK/bin/bwa
 FASTQC=$GATK/FastQC

 SEL_VARIANT="$JAVA_GATK  -T SelectVariants -R $HG19 "
 VE="$JAVA_GATK -T VariantEval -R $HG19 "
 SUM="$GATK/annovar/summarize_annovar.pl "
 CVT="$GATK/annovar/convert2annovar.pl "
 HDB="$GATK/annovar/humandb/"

FlowCellInfo=$HISEQ_ANALYSIS/FlowCellInfo
BclToFastq=$HISEQ_ANALYSIS/CASAVA1.8.2/bin/configureBclToFastq.pl

qsub_opt="-q  \
all.q@n-0-1.local,all.q@n-0-10.local,all.q@n-0-11.local,all.q@n-0-12.local,all.q@n-0-13.local,all.q@n-0-14.local,all.q@n-0-15.local,all.q@n-0-16.local,all.q@n-0-17.local,all.q@n-0-18.local,all.q@n-0-19.local,all.q@n-0-2.local,all.q@n-0-20.local,all.q@n-0-21.local,all.q@n-0-22.local,all.q@n-0-23.local,all.q@n-0-24.local,all.q@n-0-25.local,all.q@n-0-26.local,all.q@n-0-27.local,all.q@n-0-28.local,all.q@n-0-29.local,all.q@n-0-3.local,all.q@n-0-30.local,all.q@n-0-31.local,all.q@n-0-4.local,all.q@n-0-5.local,all.q@n-0-6.local,all.q@n-0-7.local,all.q@n-0-8.local,all.q@n-0-9.local \
    -m abe -M lifeng4209@gmail.com"
