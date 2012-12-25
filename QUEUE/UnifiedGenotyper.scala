package org.broadinstitute.sting.queue.qscripts.examples

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._

import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.gatk.phonehome._

import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel

  
import org.broadinstitute.sting.queue.extensions.snpeff._

/**
 * UnifiedGenotyper 
 */
class Genotyper extends QScript {
  // 'qscript' is now the same as 'ExampleUnifiedGenotyper.this'
  qscript =>


  // Required arguments.  All initialized to empty values.

  @Input(doc="The reference file for the bam files.", shortName="R")
  var referenceFile: File = _ // _ is scala shorthand for null

  @Input(doc="Bam file to genotype.", shortName="I")
  var bamFile: Seq[File] = _

  // The following arguments are all optional.

  @Input(doc="An optional file with a list of intervals to proccess.", shortName="L", required=false)
  var intervals: File = _

  @Argument(doc="A optional list of filter names.", shortName="snpfilter", required=false)
  var snpfilterNames: List[String] =  List("LowQualityDepth","MappingQuality","StrandBias","HaplotypeScoreHigh","MQRankSumLow","ReadPosRankSumLow" )

  @Argument(doc="An optional list of filter expressions.", shortName="snpfilterExpression", required=false)
  var snpfilterExpressions: List[String] = List ("QD < 2.0", "MQ < 40.0", "FS > 60.0", "HaplotypeScore > 13.0", "MQRankSum < -12.5", "ReadPosRankSum < -8.0")

  @Argument(doc="A optional list of filter names.", shortName="indelfilter", required=false)
  var indelfilterNames: List[String] =  List("LowQualityDepth","StrandBias","ReadPosRankSumLow","InbreedingCoeffLow" )

  @Argument(doc="An optional list of filter expressions.", shortName="indelfilterExpression", required=false)
  var indelfilterExpressions: List[String] = List ("QD < 2.0", "FS > 200.0",  "ReadPosRankSum < -20.0", "InbreedingCoeff < -0.8")


  // This trait allows us set the variables below in one place,
  // and then reuse this trait on each CommandLineGATK function below.
  trait UnifiedGenotyperArguments extends CommandLineGATK {
    this.reference_sequence = qscript.referenceFile
    this.intervals = if (qscript.intervals == null) Nil else List(qscript.intervals)
    // Set the memory limit to 4 gigabytes on each command.
    this.memoryLimit = 4
    this.phone_home = GATKRunReport.PhoneHomeOption.NO_ET
    this.gatk_key = "/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/TianL_email.chop.edu.key"
  }
  
  
  val queueLogDir: String = ".qlog/" // Gracefully hide Queue's output
  
  @Hidden
  @Argument(doc="How many ways to scatter/gather", fullName="scatter_gather", shortName="sg", required=false)
  var nContigs: Int = -1

  /*  Functions  */
  def call_genotypes_indel(bams: Seq[File])  {
    val genotyper = new UnifiedGenotyper with UnifiedGenotyperArguments
    val variantFilter = new VariantFiltration with UnifiedGenotyperArguments
    val evalUnfiltered = new VariantEval with UnifiedGenotyperArguments
    val evalFiltered = new VariantEval with UnifiedGenotyperArguments

    genotyper.scatterCount = nContigs 

    genotyper.input_file = bams
    genotyper.out = "unifiedgenotyper.indel.vcf"

    genotyper.genotype_likelihoods_model = GenotypeLikelihoodsCalculationModel.Model.INDEL

    evalUnfiltered.eval :+= genotyper.out
    evalUnfiltered.out = swapExt(genotyper.out, "vcf", "eval")

    variantFilter.variant = genotyper.out
    variantFilter.out = swapExt(genotyper.out, "vcf", "filtered.vcf")
    variantFilter.filterName = indelfilterNames
    variantFilter.filterExpression = indelfilterExpressions

    evalFiltered.eval :+= variantFilter.out
    evalFiltered.out = swapExt(variantFilter.out, "vcf", "eval")

    add(genotyper, evalUnfiltered)
    // Only add variant filtration to the pipeline if filters were passed in
    if (indelfilterNames.size > 0)
      add(variantFilter, evalFiltered)
}

def call_genotypes_snp(bams: Seq[File])  {
    val genotyper = new UnifiedGenotyper with UnifiedGenotyperArguments
    val variantFilter = new VariantFiltration with UnifiedGenotyperArguments
    val evalUnfiltered = new VariantEval with UnifiedGenotyperArguments
    val evalFiltered = new VariantEval with UnifiedGenotyperArguments

    genotyper.scatterCount = nContigs

    genotyper.input_file = bams
    genotyper.out = "unifiedgenotyper.snp.vcf"

    genotyper.genotype_likelihoods_model = GenotypeLikelihoodsCalculationModel.Model.SNP

    evalUnfiltered.eval :+= genotyper.out
    evalUnfiltered.out = swapExt(genotyper.out, "vcf", "eval")
    evalUnfiltered.dbsnp = "/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/hg19/dbsnp_135.hg19.vcf"

    variantFilter.variant = genotyper.out
    variantFilter.out = swapExt(genotyper.out, "vcf", "filtered.vcf")
    variantFilter.filterName = snpfilterNames
    variantFilter.filterExpression = snpfilterExpressions

    evalFiltered.eval :+= variantFilter.out
    evalFiltered.out = swapExt(variantFilter.out, "vcf", "eval")
    evalFiltered.dbsnp = "/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/hg19/dbsnp_135.hg19.vcf"

    add(genotyper, evalUnfiltered)
    // Only add variant filtration to the pipeline if filters were passed in
    //if (snpfilterNames.size > 0)
      add(variantFilter, evalFiltered)
}


  def annotate_snpEff(inVcf: File) {
        val eff = new SnpEff
        eff.config = new File("/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/snpEff_2_0_5/snpEff.config")
        eff.genomeVersion = "GRCh37.64"

	eff.inVcf = inVcf
	var snpEffout: File  = swapExt(eff.inVcf, "vcf", "snpEff.vcf")
        eff.outVcf = swapExt(eff.inVcf, "vcf", "snpEff.out")
       	 
        eff.javaClasspath = List("/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/snpEff_2_0_5/")
        eff.jarFile = "/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/snpEff_2_0_5/snpEff.jar"
    
        add(eff)
	add(varannotator(eff.inVcf,eff.outVcf,snpEffout))
}

  def script() {

	call_genotypes_snp(bamFile)
	call_genotypes_indel(bamFile)
	annotate_snpEff("unifiedgenotyper.indel.filtered.vcf")
	annotate_snpEff("unifiedgenotyper.snp.filtered.vcf")


  }

case class varannotator (inVcf: File, inSnpEffFile: File, outVcf: File) extends VariantAnnotator  {
    this.variant = inVcf
    this.snpEffFile = inSnpEffFile
    this.out = outVcf
    this.alwaysAppendDbsnpId = true
    this.D = "/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/hg19/dbsnp_135.hg19.vcf"
    this.R = "/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/hg19/hg19.fa"
    this.A = Seq("SnpEff")
    this.isIntermediate = false
    this.analysisName = queueLogDir + outVcf + ".varannotator"
    this.jobName = queueLogDir + outVcf + ".varannotator"
    this.scatterCount = nContigs
  }
}

