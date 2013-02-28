JAVA_QUEUE="/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/java/jdk1.6.0_35/bin/java -Djava.io.tmpdir=/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/temp -Xmx6g -Xms6g -jar /mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/Queue-2.3-4-g57ea19f/Queue.jar"


$JAVA_QUEUE -S $PWD/UnifiedGenotyper.scala -R $HG19 -I bam.list --scatter_gather 100  -qsub   -jobParaEnv smp   -gv graph.gv  -l DEBUG  -jobQueue long.q -memLimit 28  -jobResReq h_vmem=36g   -keepIntermediates  -L /mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/bed/ucschg19knownexons.clean.bed  --interval_padding 100  -retry 5 -pwd $PWD -run
 #-startFromScratch -run

