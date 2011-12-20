
source setup.sh

list=`ls $HISEQ/FlowCellInfo/*.csv`
sum=0
for i in ${list[*]}; do
	echo $i
	let sum=$sum+1
done

echo total: $sum
