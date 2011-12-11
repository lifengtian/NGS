#Date: Oct 27, 2011
#Proj: Saudi 1 to 13
#Path: /abings/2nd1a/results/lifeng/Saudi
p=`pwd`

annovar=$GATK/annovar
annovar_db=$GATK/annovar/humandb

### calculate the dbSNP132%
function dbsnp {

a=`ls *.gatkstandard.report`
for i in ${a[*]}; do
sn=`echo $i | cut -f1 -d '.' `
db=`head -3 $i | tail -1 | awk '{print $9}'`
nvars=`grep -v ^# $sn.gatkstandard.snp.vcf | wc -l`

echo $i $db $nvars
done

}


function annovar {
mkdir -p $p/scripts
mkdir -p $p/log
mkdir -p $p/annovar

a=`ls *.gatkstandard.snp.vcf`
for i in ${a[*]}; do
sn=`echo $i | cut -f1 -d '.' `
echo $sn
echo "
$annovar/convert2annovar.pl -format vcf4 $p/$i > $p/annovar/$sn.av
$annovar/auto_annovar.pl --ver1000g 1000g2011may --buildver hg19  -step 1,4,7-9 -model recessive --verdbsnp 132 $p/annovar/$sn.av $annovar_db 
#$annovar/summarize_annovar.pl  --ver1000g 1000g2011may --verdbsnp 132 --buildver hg19 --outfile $p/annovar/$sn.step8.sum $p/annovar/$sn.av.step8.varlist $annovar_db 
$annovar/summarize_annovar.pl  --ver1000g 1000g2011may --verdbsnp 132 --buildver hg19 --outfile $p/annovar/$sn.step2.sum $p/annovar/$sn.av.step2.varlist $annovar_db 

" > $p/scripts/$i.annovar.sh
qsub -o $p/log/$i.annovar.log $p/scripts/$i.annovar.sh

done

}

function annovar_indel {
#KS01.gatkstandard.indel.vcf
a=`ls *.gatkstandard.indel.vcf`
for i in ${a[*]}; do
s=`echo $i | cut -f1 -d '.' `
sn=$s.indel
echo $sn
echo "
$annovar/convert2annovar.pl -format vcf4 $p/$i > $p/annovar/$sn.av
$annovar/auto_annovar.pl --ver1000g 1000g2011may --buildver hg19  -step 1,4,7-9 -model recessive --verdbsnp 132 $p/annovar/$sn.av $annovar_db 
#$annovar/summarize_annovar.pl  --ver1000g 1000g2011may --verdbsnp 132 --buildver hg19 --outfile $p/annovar/$sn.step8.sum $p/annovar/$sn.av.step8.varlist $annovar_db 
$annovar/summarize_annovar.pl  --ver1000g 1000g2011may --verdbsnp 132 --buildver hg19 --outfile $p/annovar/$sn.step2.sum $p/annovar/$sn.av.step2.varlist $annovar_db 

" > $p/scripts/$i.annovar.indel.sh
qsub -o $p/log/$i.annovar.indel.log $p/scripts/$i.annovar.indel.sh

done
}



function analyze_annovar {

p=`pwd`
annovar=$p/annovar
cd $annovar
a=`ls *.step9.varlist`
echo hom
for i in ${a[*]}; do
   hom=`grep hom $i| cut -f2 | cut -f1 -d ':'`
   echo $i $hom
done


echo comp het
for i in ${a[*]}; do
	chet=`grep het $i | cut -f2 | cut -f1 -d ':' | sort | uniq -c | sort -k1,1n | awk '{if ( $1>1 ) print $2}'`
   	echo $i $chet

done


}

function coverage {
a=`ls *.dedup.bam.coverage`

first=0
for i in ${a[*]}; do
	header=`grep BAIT_SET $i| cut -f1,12,22,29,30`
	cov=`grep ^SureSelect50mbclean $i|  cut -f1,12,22,29,30  `
	if [ $first -eq  0 ]; then
		echo sample $header
		let first=$first+1
	fi
	
    echo $i $cov
done
}


function find_sharedvars {
a=(
"SAU_08"
"SAU_09"
"SAU_10"
"SAU_11"
"SAU_12"
)


b=(
"SAU_22"
"SAU_23"
"SAU_24"
"SAU_25"
"SAU_26"
)
echo hom
for i in ${!a[*]}; do
	variants=`perl SharedVars.pl hom.txt ${a[$i]} ${b[$i]}`
	echo ============== ${a[$i]} ${b[$i]} $variants
	for j in ${variants[*]}; do
		grep $j ../Saudi/annovar/${a[$i]}*step9*
		grep $j ../Saudi3/annovar/${b[$i]}*step9*
		echo --------------
	done
done

echo comphet
for i in ${!a[*]}; do
    variants=`perl SharedVars.pl comphet.txt ${a[$i]} ${b[$i]}`
    echo ============== ${a[$i]} ${b[$i]} $variants
    for j in ${variants[*]}; do
        grep $j ../Saudi/annovar/${a[$i]}*step9*
        grep $j ../Saudi3/annovar/${b[$i]}*step9*
        echo --------------
    done
done
}

annovar
annovar_indel
#analyze_annovar
#dbsnp
#find_sharedvars
