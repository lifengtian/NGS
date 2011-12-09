 
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


