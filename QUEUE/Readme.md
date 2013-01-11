GATK QUEUE scripts live here.


VariantAnnotator:    this.scatterCount = nContigs

VariantEval: scatterCount NOT work; this.num_threads = 4 works.

BQSR plot: http://gatkforums.broadinstitute.org/discussion/1297/no-plots-generated-by-the-baserecalibrator-walker
	git/gatk/public/R/scripts/org/broadinstitute/sting/utils/recalibration/BQSR.R

Rscript: Rscript /mnt/isilon/cag/ngs/hiseq/respublica/pipeline/git/gatk/public/R/scripts/org/broadinstitute/sting/queue/util/queueJobReport.R  report.txt report.pdf
