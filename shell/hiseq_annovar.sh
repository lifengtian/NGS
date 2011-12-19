function annovar {
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
qsub -o $p/logs/$i.annovar.log $p/scripts/$i.annovar.sh

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
qsub -o $p/logs/$i.annovar.indel.log $p/scripts/$i.annovar.indel.sh

done
}

annovar
annovar_indel
