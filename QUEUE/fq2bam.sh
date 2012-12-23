smlists=(
4361218280
#9235681606
#2513613612
)

JAVA_QUEUE="/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/java/jdk1.6.0_35/bin/java -Djava.io.tmpdir=/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/temp -Xmx6g -Xms6g -jar /mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/Queue-2.2-16-g9f648cb/Queue.jar"

for sm in ${smlists[*]}; do
        fq1=/mnt/isilon/cag/ngs/hiseq/Miami/FASTQ.cagid/$sm/pe_1.fq.gz
        fq2=/mnt/isilon/cag/ngs/hiseq/Miami/FASTQ.cagid/$sm/pe_2.fq.gz
        out=$PWD/$sm.bam

	mkdir -p $sm/fq1
	mkdir -p $sm/fq2
	mkdir -p $sm/bam
	mkdir -p $sm/log

	echo "
	zcat $fq1 | split -l 1000000 - $sm/fq1/fq
	" | qsub -cwd -N split$sm

	echo "
	zcat $fq2 | split -l 1000000 - $sm/fq2/fq
	" | qsub -cwd -N split$sm

	echo "sleep 10" | qsub -cwd -hold_jid split$sm -sync y

	cd $sm/fq1
	for j in fq* ; do 
	echo "
		date
		$JAVA_BIN -Xmx4g -jar $PICARD/FastqToSam.jar QUALITY_FORMAT=Standard FASTQ=$PWD/$j  FASTQ2=$PWD/../fq2/$j  OUTPUT=$PWD/../bam/$sm.$j.bam SAMPLE_NAME=$sm LB=LB PU=PU PL=Illumina CN=CAG SORT_ORDER=unsorted CREATE_INDEX=TRUE
		date
	" | qsub -cwd -S /bin/sh -N f$sm -l mem_free=10G,h_vmem=12G -o $PWD/../log/f$sm.$j.o.log -e $PWD/../log/f$sm.$j.err.log

	done
	cd ../..

	
	echo " ls -lt $PWD/$sm/bam/*bam | awk '{print \$9}' > $sm/bam.list 

	$JAVA_QUEUE -S $PWD/align_reads.scala -i $sm/bam.list --scatter_gather 40  -qsub   -jobParaEnv smp   --bwa_threads 4  -gv $PWD/graph.gv  -l DEBUG  -jobQueue short.q -memLimit 20  -jobResReq h_vmem=22g -run  --samplename $sm -outputDir $sm
	" > $sm.sh
	 qsub -cwd -S /bin/sh -hold_jid f$sm -N pipeline -V -l mem_free=10G,h_vmem=12G  $sm.sh
done


