JAVA_QUEUE="/usr/bin/java -XX:ParallelGCThreads=1 -Djava.io.tmpdir=/scr1/cag_lab2/tianl/ -Xmx12g -Xms12g -jar /mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/Queue-2.5-2-gf57256b/Queue.jar"

#cag chr1 to chr22, chrX, chrY, chrM
HG19=/mnt/isilon/cag/ngs/hiseq/tianl/cagpipelines/cag.hg19/cag.hg19.fasta
dbsnp=/mnt/isilon/cag/ngs/hiseq/tianl/cagpipelines/cag.hg19/dbsnp_135.hg19.vcf
mills=/mnt/isilon/cag/ngs/hiseq/tianl/cagpipelines/cag.hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf
omni=/mnt/isilon/cag/ngs/hiseq/tianl/cagpipelines/cag.hg19/1000G_omni2.5.hg19.sites.vcf
hapmap=/mnt/isilon/cag/ngs/hiseq/tianl/cagpipelines/cag.hg19/hapmap_3.3.hg19.vcf

$JAVA_QUEUE -S $PWD/UnifiedGenotyper.scala -R $HG19 -I bam.list --scatter_gather 50  -qsub   -jobParaEnv smp   -gv graph.gv  -l DEBUG  -jobQueue long.q -memLimit 28  -jobResReq h_vmem=36g   -keepIntermediates  -L /mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/bed/ensgene.ucsc.exons.bed  --interval_padding 100  -retry 3 -pwd $PWD -run  -D $dbsnp \
-Omni $omni \
-Hapmap $hapmap \
-Mills $mills 
#-startFromScratch -run

