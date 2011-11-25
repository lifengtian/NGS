
#current path
p=`pwd`

#sample name
sm=??

#forward and reverse reads FASTQ(gzipped FASTQ) file name
r1=$p/??
r2=$p/??

date
echo $sm 

#align
$BWA aln  -t 12 $REF $r1.fastq.gz > $r1.sai
$BWA aln  -t 12 $REF $r2.fastq.gz > $r2.sai

$BWA sampe -r "@RG\tID:$sm\tSM:$sm" $REF $r1.sai $r2.sai $r1.fastq.gz $r2.fastq.gz  | $SAMTOOLS view -S -b -o $p/$sm.bam -


$JAVA -jar $PICARD/SortSam.jar I=$p/$sm.bam O=$p/$sm.sorted.bam SO=coordinate VALIDATION_STRINGENCY=SILENT TMP_DIR=. CREATE_INDEX=true
$SAMTOOLS flagstat $p/$sm.sorted.bam > $p/$sm.sorted.bam.flagstat
$SAMTOOLS depth $p/$sm.sorted.bam | perl GENOMECOVERAGE.PL > $p/$sm.sorted.bam.genomecoverage
$FASTQC -o QCreport -f $p/$sm.sorted.bam
$JAVA -jar $PICARD/MarkDuplicates.jar I=$p/$sm.sorted.bam O=$p/$sm.dedup.bam M=$p/$sm.metric VALIDATION_STRINGENCY=SILENT TMP_DIR=. CREATE_INDEX=true
$SAMTOOLS flagstat $p/$sm.dedup.bam > $p/$sm.dedup.bam.flagstat
$SAMTOOLS depth $p/$sm.dedup.bam | perl GENOMECOVERAGE.PL > $p/$sm.dedup.bam.genomecoverage

date

