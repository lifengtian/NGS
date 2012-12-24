JAVA_QUEUE="/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/java/jdk1.6.0_35/bin/java -Djava.io.tmpdir=/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/temp -Xmx6g -Xms6g -jar /mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/Queue-2.2-16-g9f648cb/Queue.jar"

# I get rid of sampleName
$JAVA_QUEUE -S $PWD/UnifiedGenotyper.scala -R $HG19 -I bam.list --scatter_gather 80  -qsub   -jobParaEnv smp  -l DEBUG  -jobQueue long.q -memLimit 20  -jobResReq h_vmem=22g   -keepIntermediates

