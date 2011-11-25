
a=(1 )
# 2 3)
for i in ${a[*]}; do
/work/1a/cagapps/jdk/jdk1.6.0_22/bin/java -jar /work/1a/cagapps/gatk/picard/AddOrReplaceReadGroups.jar   I=$i.bam O=$i.new.bam SO=coordinate RGID=20110829$i RGLB=50x35RRBC RGPL=SOLiD RGPU=bioscope-pairing RGSM=KS01 CREATE_INDEX=true

done
