package org.broadinstitute.sting.queue.qscripts.examples

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.snpeff._

import org.broadinstitute.sting.gatk.phonehome._


class SnpEffWalker extends QScript {
  // Create an alias 'qscript' to be able to access variables
  // 'qscript' is now the same as 'SnpEffWalker.this'

  qscript =>


  // Required arguments.  All initialized to empty values.

//  @Input(doc="The reference file for the bam files.", shortName="R")
//  var referenceFile: File = _ // _ is scala shorthand for null

  @Input(doc="VCF file to annotate.", shortName="I")
  var vcfFile: File = _

  // The following arguments are all optional.

  @Input(doc="An optional file with a list of intervals to proccess.", shortName="L", required=false)
  var intervals: File = _

  @Argument(doc="A optional list of filter names.", shortName="filter", required=false)
  var filterNames: List[String] = Nil // Nil is an empty List, versus null which means a non-existent List.

  @Argument(doc="An optional list of filter expressions.", shortName="filterExpression", required=false)
  var filterExpressions: List[String] = Nil

  // This trait allows us set the variables below in one place,
  // and then reuse this trait on each CommandLineGATK function below.
  trait UnifiedGenotyperArguments extends CommandLineGATK {
    this.phone_home = GATKRunReport.PhoneHomeOption.NO_ET
   this.gatk_key = "/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/TianL_email.chop.edu.key"
 }

   // global
   val queueLogDir: String = ".qlog/" // Gracefully hide Queue's output

  def script() {
	val eff = new SnpEff 
	eff.config = new File("/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/snpEff_2_0_5/snpEff.config")
	eff.genomeVersion = "GRCh37.64"
	eff.inVcf = vcfFile 
	var snpEffout: File  = swapExt(eff.inVcf, "vcf", "snpEff.out")
	eff.outVcf = swapExt(eff.inVcf, "vcf", "snpEff.vcf")
  	eff.javaClasspath = List("/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/snpEff_2_0_5/")
	eff.jarFile = "/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/snpEff_2_0_5/snpEff.jar"  
  
	add(eff)
	add(varannotator(eff.inVcf, snpEffout, eff.outVcf) )
	// $JAVA_GATK  -T VariantAnnotator -R $HG19 -A SnpEff -V $input_vcf --snpEffFile $TEMP/$snpEff_out_vcf -o $output_vcf -D $PIPELINE/gatk/hg19/dbsnp_135.hg19.vcf -alwaysAppendDbsnpId


  }

  case class varannotator (inVcf: File, inSnpEffFile: File, outVcf: File) extends VariantAnnotator with UnifiedGenotyperArguments  {
    this.variant = inVcf
    this.snpEffFile = inSnpEffFile
    this.out = outVcf
    this.alwaysAppendDbsnpId = true
    this.D = "/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/hg19/dbsnp_135.hg19.vcf"
    this.R = "/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/hg19/hg19.fa"
    this.isIntermediate = false
    this.analysisName = queueLogDir + outVcf + ".varannotator"
    this.jobName = queueLogDir + outVcf + ".varannotator"
  }

}

