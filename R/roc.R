#roc for BWA
#no of reads: 200,000
roc.ag50<-read.table("/Users/tianl/lec3/r12.bwa.ag50.roc",header=FALSE)
tp.ag50<-roc.ag50$V2/200000
fp.ag50<-roc.ag50$V3/200000
plot(fp.ag50,tp.ag50,col="blue",main="ROC of BWA",xlab="False positive rate (%)",ylab="True positive rate (%)",type="b",xlim=c(0,0.002))

roc.ag1<-read.table("/Users/tianl/lec3/r12.bwa.ag1.roc",header=FALSE)
tp.ag1<-roc.ag1$V2/200000
fp.ag1<-roc.ag1$V3/200000

par(new=TRUE)
points(fp.ag1,tp.ag1,col="green",type="b",xlab="",ylab="")




### VCF ROC

## GATK
## No of variants: grep -v ^# r12.bwa.gatk.vcf | wc -l
##         17140
roc.gatk<-read.table("/Users/tianl/lec3/r12.bwa.gatk.roc",header=FALSE)

tp.gatk<-roc.gatk$V2/17140
fp.gatk<-roc.gatk$V3/17140

plot(fp.gatk,tp.gatk,col="blue",main="ROC of GATK vs SAMtools",xlab="False positive rate (%)",ylab="True positive rate (%)",type="b",xlim=c(0,0.15))
par(new=TRUE)
  
## No of variants: grep -v ^# r12.bwa.samtools.vcf | wc -l
##          17654
  
roc.samtools<-read.table("/Users/tianl/lec3/r12.bwa.samtools.roc",header=FALSE)

tp.samtools<-roc.samtools$V2/17654
fp.samtools<-roc.samtools$V3/17654

points(fp.samtools,tp.samtools,col="green")
lines(fp.samtools,tp.samtools,col="green")


