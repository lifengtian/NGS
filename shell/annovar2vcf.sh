in_csv=$1
in_vcf=$2
out=$3

perl csv_tab.pl $in_csv  | sed -e 's/;/:/g' -e 's/=/:/g'  | sort -k1,1 -k2,2n | bgzip -c > $in_csv.gz
tabix -b 2 -e 3 -s 1 -f $in_csv.gz 
cat $in_vcf | $GATK/vcftools/bin/vcf-annotate -a $in_csv.gz -c CHROM,FROM,TO,INFO/Func,INFO/exon,INFO/Gene,INFO/ExonicFunc,INFO/AAChange,INFO/Conserved,INFO/SegDup,INFO/1000g,-,-,INFO/dbSNP132,INFO/SIFT,INFO/PolyPhen2,INFO/LJB_PhyloP,INFO/LJB_MutationTaster,INFO/LJB_LRT,INFO/Ref,INFO/Obs,-,-,-,-,- -d annovar_description.txt > $out 
#bgzip -c AS02.gatkstandard.snp.annovar.vcf > AS02.gatkstandard.snp.annovar.vcf.gz
#tabix -p vcf AS02.gatkstandard.snp.annovar.vcf.gz 
