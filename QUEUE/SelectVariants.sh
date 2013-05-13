JAVA_QUEUE="/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/java/jdk1.6.0_35/bin/java -Djava.io.tmpdir=/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/temp -Xmx6g -Xms6g -jar /mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/Queue-2.2-16-g9f648cb/Queue.jar"

#JAVA_QUEUE="/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/java/jdk1.6.0_35/bin/java -Djava.io.tmpdir=/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/temp -Xmx6g -Xms6g -jar /mnt/isilon/cag/ngs/hiseq/respublica/pipeline/git/gatk/dist/Queue.jar" 

vcf=$1
type=$2
sn=$3

if [ -n "$sn" ] ; then

echo sample $3
out="SID"$3.$type.vcf

# I get rid of sampleName
$JAVA_QUEUE -S $PWD/SelectVariants.scala -R $HG19 -V $vcf  -qsub   -jobParaEnv smp   -l DEBUG  -jobQueue long.q -memLimit 28  -jobResReq h_vmem=36g   -keepIntermediates -run --out $out --sampleName $3

else
for i in $(grep ^#CHROM $vcf | cut -f10-1000) ; do
out="SID"$i.$type.vcf
echo $out
# I get rid of sampleName
$JAVA_QUEUE -S $PWD/SelectVariants.scala -R $HG19 -V $vcf  -qsub   -jobParaEnv smp   -l DEBUG  -jobQueue long.q -memLimit 28  -jobResReq h_vmem=36g   -keepIntermediates -run --out $out --sampleName $i 
#-L /mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/bed/SureSelectV2.picard --interval_padding 100
 #-startFromScratch -run
#--scatter_gather 400


done

fi
