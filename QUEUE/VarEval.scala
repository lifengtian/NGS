package org.broadinstitute.sting.queue.qscripts.examples

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._

import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.gatk.phonehome._

import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel

  
import org.broadinstitute.sting.queue.extensions.snpeff._

class VariantEvalScript extends QScript {
  qscript =>


  // Required arguments.  All initialized to empty values.

  @Input(doc="The reference file for the bam files.", shortName="R")
  var referenceFile: File = _ // _ is scala shorthand for null

  @Input(doc="eval vcf ", shortName="Eval")
  var evalFile: File = _

  @Input(doc="comp vcf ", shortName="Comp")
  var compFile: File = _

  // The following arguments are all optional.

  @Input(doc="An optional file with a list of intervals to proccess.", shortName="L", required=false)
  var intervals: File = _






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


  def script() {

	add ( evalfunc(evalFile, compFile) )

  }

case class evalfunc (evalVcf: File, compVcf: File) extends VariantEval { 
    this.eval :+= evalVcf 
    this.comp :+= compVcf
    this.out = swapExt(evalVcf, "vcf", "eval")
    this.D = "/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/hg19/dbsnp_135.hg19.vcf"
    this.R = "/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/hg19/hg19.fa"
    this.isIntermediate = false
    this.analysisName = queueLogDir + evalVcf + ".varannotator"
    this.jobName = queueLogDir + evalVcf + ".varannotator"
    //this.scatterCount = nContigs
	this.num_threads = 4
  }
}

