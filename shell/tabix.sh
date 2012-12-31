
vcf=$1
bgzip -c $vcf > $vcf.gz
tabix -p vcf $vcf.gz

