
source setup.sh

list=`ls $HISEQ/SampleSheet/*.csv`
#sum=0
#for i in ${list[*]}; do
#	sids=$sids" "`cat $i | grep -v FCID | cut -f3 -d, | sort | uniq` 		
#done

#echo $sids | uniq

cat $HISEQ/SampleSheet/*.csv | grep -v FCID | cut -f3 -d, | sort | uniq
