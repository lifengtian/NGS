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

  @Input(doc="VCF file ", shortName="V")
  var vcfFile: File = _

  @Output(doc="output VCF file" , fullName="out")
  var outFile: File = _

  @Argument(doc="sample name", fullName="sampleName")
  var sampleName: String = _
  // The following arguments are all optional.

  @Input(doc="An optional file with a list of intervals to proccess.", shortName="L", required=false)
  var intervals: File = _

  @Argument(doc="interval padding ", fullName="interval_padding", shortName="ip",required=false)
  var ip: Int = _

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
//    this.intervals = if (qscript.intervals == null) Nil else List(qscript.intervals)
    this.memoryLimit = 8 
    this.phone_home = GATKRunReport.PhoneHomeOption.NO_ET
    this.gatk_key = "/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/TianL_email.chop.edu.key"
//    this.interval_padding = ip
  }
  
  
  val queueLogDir: String = ".qlog/" // Gracefully hide Queue's output
  
  @Hidden
  @Argument(doc="How many ways to scatter/gather", fullName="scatter_gather", shortName="sg", required=false)
  var nContigs: Int = -1

  /*  Functions  */

  def script() {

	add ( selectsample ( vcfFile, outFile, sampleName ) )
  }

case class selectsample (inVcf: File, outVcf: File, samplename: String) extends SelectVariants {
    this.variant = inVcf
    this.out = outVcf
    this.sample_name = Seq(samplename)
    this.excludeNonVariants = true
    this.excludeFiltered = true
    this.R = qscript.referenceFile 
    this.isIntermediate = false
    this.analysisName = queueLogDir + outVcf + ".varannotator"
    this.jobName = queueLogDir + outVcf + ".varannotator"
  }
}

