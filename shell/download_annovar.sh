# download Annovar databases 
# default is by openbioinformatics.org, NOT ucsc
# Lifeng Tian
# 2011


list=(
gene
1000g2011may
snp132
snp129
avsift
ljb_pp2
ljb_mt
ljb_phylop
ljb_lrt
cg69
knowngene
ensgene
omimgene
dgv
segdup
mce44way
gwascatalog
phastConsElements46way
)

for i in ${list[*]};do
        #annotate_variation.pl -downdb  -buildver hg19 $i  humandb/
        $GATK/annovar/annotate_variation.pl -buildver hg19 --downdb --webfrom annovar $i $GATK/annovar/humandb/
done

