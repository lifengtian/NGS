#!/bin/sh

# THIS script is called by create_vcfs.pl, which passes parameters (target, bamlist, etc.) via gatk.sh
# It retrieves all BAMs for each sample (by SIDs)
#    >1 BAMs: Merge and dedup them
#    =1 BAM: make links to bam and bai
# It runs GATK on the dedupped bams. Variants are called together.




################################################## FUNCTION ########################################
# include functions including
# gatk_pipeline

source $HISEQ/NGS/shell/function.sh







################################################### MAIN ############################################

# input: running path
# e.g.:   bash hiseq_gatk_resolute.sh /path/to/config.yml
p=$1


source $p/gatk.sh
# Resolute queue
source $HISEQ/NGS/shell/setup.sh


mkdir -p $p/scripts
mkdir -p $p/logs
mkdir -p $p/BAM
mkdir -p $p/VCF


for sid in ${bamlist[*]}; do
    # this is a cheap way to find all processed BAMs for a sample
    # a better way is to request it through the SampleSheet object. will do it later
    # SampleSheet->get_bams_by_sid('sid')
    list=`find $HISEQ/BAM  -name "*.dedup.bam" -print | grep Sample_$sid`

    I=""
    count=0

    for i in ${list[*]}; do
    I=$I" I=$i"
    let count=$count+1
    done

    echo ---------------------- $sid  from $I to $p ----------------------------
    date
    out=$p/BAM/$sid.bam
    dedup=$p/BAM/$sid.dedup.bam
    M=$p/BAM/$sid.metrics

    timestamp=`date +%s`

    if [ ! -f $out ]; then
    # if there is more than 1 bam for this sample
    # we will merge and dedup them
    # the newly created BAM will live in the Project space.
        if [ "$count" -gt 1 ]; then
            echo "
            $JAVA_BIN -jar $PICARD/MergeSamFiles.jar $I O=$out  VALIDATION_STRINGENCY=SILENT TMP_DIR=$p/temp CREATE_INDEX=true
            $JAVA_BIN -jar $PICARD/MarkDuplicates.jar REMOVE_DUPLICATES=true I=$out O=$dedup M=$M VALIDATION_STRINGENCY=SILENT TMP_DIR=$p/temp CREATE_INDEX=true ASSUME_SORTED=true
            $SAMTOOLS flagstat $out > $out.flagstat
            $SAMTOOLS flagstat $dedup > $dedup.flagstat
            $SAMTOOLS depth $dedup | perl $GENOMECOVERAGE > $dedup.genomecoverage
            " > $p/scripts/$sid.merge.$timestamp.sh

            if [ "$target" != "genome" ]; then
            echo "
            $JAVA_BIN -jar $PICARD/CalculateHsMetrics.jar I=$dedup O=$dedup.target.coverage BI=$GATK/bed/$target.picard TI=$GATK/bed/$target.picard TMP_DIR=$p/temp VALIDATION_STRINGENCY=SILENT
            " >> $p/scripts/$sid.merge.$timestamp.sh
            fi

            job=`qsub $queue  -e $p/logs/$sid.merge.$timestamp.log -o $p/logs/$sid.merge.$timestamp.log $p/scripts/$sid.merge.$timestamp.sh`
            hold_jid=`echo $job | awk '{print $3}'`
            d=`date`
            echo $d: Merge-Dedup-flagstat $hold_jid
        else 
            ln -s ${list[0]} $dedup
            in_bai=`echo ${list[0]} | sed 's/bam$/bai/'`
            out_bai=`echo $dedup | sed 's/bam$/bai/'`
            ln -s $in_bai $out_bai 
        fi # -gt 1
    else
        echo $out exist! Exit!
        exit
    fi

    sleep 1
done



gatk_pipeline

