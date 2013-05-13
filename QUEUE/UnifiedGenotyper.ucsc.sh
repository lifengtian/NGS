JAVA_QUEUE="/usr/bin/java -XX:ParallelGCThreads=1 -Djava.io.tmpdir=/scr1/cag_lab2/tianl/ -Xmx6g -Xms6g -jar /mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/Queue-2.4-9-g532efad/Queue.jar"

HG19=/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/bundle/2.2/hg19ucsc/ucsc.hg19.fasta
dbsnp=/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/bundle/2.2/hg19ucsc/dbsnp_137.hg19.vcf

$JAVA_QUEUE -S $PWD/UnifiedGenotyper.scala -R $HG19 -I bam.list --scatter_gather 50  -qsub   -jobParaEnv smp   -gv graph.gv  -l DEBUG  -jobQueue long.q -memLimit 28  -jobResReq h_vmem=36g   -keepIntermediates  -L /mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/bed/ensgene.ucsc.exons.bed  --interval_padding 100  -retry 3 -pwd $PWD -run -D $dbsnp 
 #-startFromScratch -run

