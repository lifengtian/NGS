#JAVA_QUEUE="/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/java/jdk1.6.0_35/bin/java -Djava.io.tmpdir=/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/temp -Xmx6g -Xms6g -jar /mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/Queue-2.3-9-ge5ebf34/Queue.jar"
JAVA_QUEUE="/usr/bin/java -XX:ParallelGCThreads=1 -Djava.io.tmpdir=/scr1/cag_lab2/tianl/ -Xmx6g -Xms6g -jar /mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/Queue-2.3-9-ge5ebf34/Queue.jar"


$JAVA_QUEUE -S $PWD/UnifiedGenotyper.vqsr.scala -R $HG19 -I bam.list --scatter_gather 200  -qsub   -jobParaEnv smp   -gv graph.gv  -l DEBUG  -jobQueue long.q -memLimit 28  -jobResReq h_vmem=36g   -keepIntermediates  -L /mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/bed/ensgene.ucsc.exons.bed  --interval_padding 100  -retry 3 -pwd $PWD -run
 #-startFromScratch -run

