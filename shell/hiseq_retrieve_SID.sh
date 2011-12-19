

## sid, the only argument required
sids=$@
p=`pwd`

mkdir -p $p/scripts
mkdir -p $p/logs
mkdir -p $p/BAM
mkdir -p $p/VCF

for sid in ${sids[*]}; do
list=`find $HISEQ/BAM  -name "*.dedup.bam" -print | grep Sample_$sid`

for i in ${list[*]}; do
	I=$I" I=$i"
	del=$del" $i "
done

echo mergeing $I
echo
echo delete $del
echo
out=$p/BAM/$sid.bam 
dedup=$p/BAM/$sid.dedup.bam
M=$p/BAM/$sid.metrics

timestamp=`date +%s`

if [ ! -f $out ]; then
	echo "
		$JAVA_BIN -jar $PICARD/MergeSamFiles.jar $I O=$out  VALIDATION_STRINGENCY=SILENT TMP_DIR=$p/temp CREATE_INDEX=true
		$JAVA_BIN -jar $PICARD/MarkDuplicates.jar REMOVE_DUPLICATES=true I=$out O=$dedup M=$M VALIDATION_STRINGENCY=SILENT  CREATE_INDEX=true ASSUME_SORTED=true
		$SAMTOOLS flagstat $out > $out.flagstat
		$SAMTOOLS flagstat $dedup > $dedup.flagstat
		$SAMTOOLS depth $dedup | perl $GENOMECOVERAGE > $dedup.genomecoverage

	     " > $p/scripts/$sid.merge.$timestamp.sh
	qsub $queue -V -e $p/logs/sid/$sid.merge.$timestamp.log -o $p/logs/sid/$sid.merge.$timestamp.log $p/scripts/$sid.merge.$timestamp.sh
else
	echo "$out exist. Exit!"
fi
sleep 2
done

