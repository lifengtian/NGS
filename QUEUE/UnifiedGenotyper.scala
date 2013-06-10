package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.function
import org.broadinstitute.sting.queue.extensions.gatk._

import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.gatk.phonehome._

import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel

  
import org.broadinstitute.sting.queue.extensions.snpeff._

import org.broadinstitute.sting.queue.util.VCF_BAM_utilities
import collection.JavaConversions._   
import net.sf.samtools.SAMFileReader
import net.sf.samtools.SAMFileHeader.SortOrder

import org.broadinstitute.sting.queue.util.QScriptUtils
import org.broadinstitute.sting.queue.function.ListWriterFunction

import org.broadinstitute.sting.gatk.arguments.ValidationExclusion

class StandardUnifiedGenotyper extends QScript {
  qscript =>

  // Required arguments.  All initialized to empty values.

  @Input(doc="The reference file for the bam files.", shortName="R")
  var referenceFile: File = _ 

  @Input(doc="The reference file for the dbSNP files.", shortName="D")
  var dbsnp_file: File = _

  @Input(doc="A file contains Bam files to genotype.", shortName="I")
  var bamFile: File = _

  // The following arguments are all optional.

  @Input(doc="An optional file with a list of intervals to proccess.", shortName="L", required=false)
  var intervals: File = _

  @Argument(shortName="pwd",required=true)
  var pwd: String = _

  @Argument(doc="interval padding ", fullName="interval_padding", shortName="ip",required=false)
  var ip: Int = 100 

  @Argument(doc="A optional list of filter names.", shortName="snpfilter", required=false)
  var snpfilterNames: List[String] =  List("LowQualityDepth","MappingQuality","StrandBias","HaplotypeScoreHigh","MQRankSumLow","ReadPosRankSumLow" )

  @Argument(doc="An optional list of filter expressions.", shortName="snpfilterExpression", required=false)
  var snpfilterExpressions: List[String] = List ("QD < 2.0", "MQ < 40.0", "FS > 60.0", "HaplotypeScore > 13.0", "MQRankSum < -12.5", "ReadPosRankSum < -8.0")

  @Argument(doc="A optional list of filter names.", shortName="indelfilter", required=false)
  var indelfilterNames: List[String] =  List("LowQualityDepth","StrandBias","ReadPosRankSumLow","InbreedingCoeffLow" )

  @Argument(doc="An optional list of filter expressions.", shortName="indelfilterExpression", required=false)
  var indelfilterExpressions: List[String] = List ("QD < 2.0", "FS > 200.0",  "ReadPosRankSum < -20.0", "InbreedingCoeff < -0.8")

  @Argument(doc="Do you want the combined VCF?", shortName="combine", required=false)
  var combineFlag : Boolean = false 
 
  @Argument(doc="Do you want to selectsamples into individual VCFs?", shortName="selectsamples", required=false)
  var selectsamplesFlag : Boolean = true


  // This trait allows us set the variables below in one place,
  // and then reuse this trait on each CommandLineGATK function below.
  trait UnifiedGenotyperArguments extends CommandLineGATK {
    this.reference_sequence = qscript.referenceFile
    this.intervals = if (qscript.intervals == null) Nil else List(qscript.intervals)
    this.memoryLimit = 8 
    this.phone_home = GATKRunReport.PhoneHomeOption.NO_ET
    this.gatk_key = "/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/TianL_email.chop.edu.key"
    this.interval_padding = ip
  }
  
  // Global  
  val queueLogDir: String = ".qlog/" // Gracefully hide Queue's output

  @Hidden
  @Argument(doc="How many ways to scatter/gather", fullName="scatter_gather", shortName="sg", required=false)
  var nContigs: Int = -1

  @Hidden
  @Argument(doc="only run snp", shortName="onlysnp", required=false)
  var onlysnp: Boolean = false 


  /*  Functions  */
  def call_genotypes_indel(bams: Seq[File], outVcf: File, waitforVcf: File)  {
    val genotyper = new myUGindel(bams, outVcf, waitforVcf) with UnifiedGenotyperArguments
    val variantFilter = new VariantFiltration with UnifiedGenotyperArguments
    val evalUnfiltered = new VariantEval with UnifiedGenotyperArguments
    val evalFiltered = new VariantEval with UnifiedGenotyperArguments

//    genotyper.scatterCount = nContigs 
 //   genotyper.num_cpu_threads_per_data_thread = 2
   // genotyper.input_file = bams
   // genotyper.out = outVcf 

    //genotyper.genotype_likelihoods_model = GenotypeLikelihoodsCalculationModel.Model.INDEL

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

def call_genotypes_snp(bams: Seq[File], outVcf: File)  {
    val genotyper = new myUGsnp(bams, outVcf) with UnifiedGenotyperArguments
    val variantFilter = new VariantFiltration with UnifiedGenotyperArguments
    val evalUnfiltered = new VariantEval with UnifiedGenotyperArguments
    val evalFiltered = new VariantEval with UnifiedGenotyperArguments

    //genotyper.scatterCount = nContigs
    //genotyper.num_cpu_threads_per_data_thread = 2
    //genotyper.input_file = bams
    //genotyper.out = outVcf

    //genotyper.genotype_likelihoods_model = GenotypeLikelihoodsCalculationModel.Model.SNP

    evalUnfiltered.eval :+= genotyper.out
    evalUnfiltered.out = swapExt(genotyper.out, "vcf", "eval")
    evalUnfiltered.dbsnp = dbsnp_file 

    variantFilter.variant = genotyper.out
    variantFilter.out = swapExt(genotyper.out, "vcf", "filtered.vcf")
    variantFilter.filterName = snpfilterNames
    variantFilter.filterExpression = snpfilterExpressions

    evalFiltered.eval :+= variantFilter.out
    evalFiltered.out = swapExt(variantFilter.out, "vcf", "eval")
    evalFiltered.dbsnp = dbsnp_file 

    add(genotyper, evalUnfiltered)
    // Only add variant filtration to the pipeline if filters were passed in
    //if (snpfilterNames.size > 0)
      add(variantFilter, evalFiltered)
}


  def annotate_snp (inVcf: File) {
        val eff = new SnpEff
        eff.config = new File("/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/snpEff_2_0_5/snpEff.config")
        eff.genomeVersion = "GRCh37.64"


	eff.inVcf = inVcf
	val snpEffout  = swapExt(eff.inVcf, "vcf", "snpEff.vcf") 
        eff.outVcf = swapExt(eff.inVcf, "vcf", "snpEff.out") 
       	 
        eff.javaClasspath = List("/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/snpEff_2_0_5/")
        eff.jarFile = "/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/snpEff_2_0_5/snpEff.jar"
   
	val annovarout = swapExt(snpEffout, "vcf", "annovar.vcf") 
	add(eff, varannotator(eff.inVcf,eff.outVcf,snpEffout), annovar_snp (snpEffout, annovarout))



}

  def annotate_indel (inVcf: File) {
        val eff = new SnpEff
        eff.config = new File("/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/snpEff_2_0_5/snpEff.config")
        eff.genomeVersion = "GRCh37.64"

        eff.inVcf = inVcf
        val snpEffout  = swapExt(eff.inVcf, "vcf", "snpEff.vcf")
	eff.outVcf = swapExt(eff.inVcf, "vcf", "snpEff.out")

        eff.javaClasspath = List("/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/snpEff_2_0_5/")
        eff.jarFile = "/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/gatk/snpEff_2_0_5/snpEff.jar"

        val annovarout = swapExt(snpEffout, "vcf", "annovar.vcf")
        add(eff, varannotator(eff.inVcf,eff.outVcf,snpEffout), annovar_indel (snpEffout, annovarout))
}



def findSampleIDsFromBAMs ( bams: Seq[File] ) : List[String] = {
    val sampleTable = scala.collection.mutable.Map.empty[String, Seq[File]]

    for (bam <- bams) {
      val samReader = new SAMFileReader(bam)
      val header = samReader.getFileHeader
      val readGroups = header.getReadGroups
	
      // only allow one sample per file. Bam files with multiple samples would require pre-processing of the file
      // with PrintReads to separate the samples. Tell user to do it himself!
      assert(!QScriptUtils.hasMultipleSamples(readGroups), "The pipeline requires that only one sample is present in a BAM file. Please separate the samples in " + bam)

      // Fill out the sample table with the readgroups in this file
      for (rg <- readGroups) {
        val sample = rg.getSample
        if (!sampleTable.contains(sample))
         {
		 sampleTable(sample) = Seq(bam) 
      	 } 
	}
    }

	var sampleList : List[String] = Nil

	for ((sample, bamList) <- sampleTable.toMap ) {
		sampleList = sampleList :+ sample
	}

	//var cohortFile = new File("/mnt/isilon/cag/ngs/hiseq/Miami/RUN/IMF/VCF/test2/sampleName.list")
        //add(writeList(sampleList, cohortFile))

	sampleList
  }




  def script() {
	val snpFile="unifiedgenotyper.snp.vcf"
	val indelFile="unifiedgenotyper.indel.vcf"

	var indel_filtered_vcfFile=swapExt(indelFile,"vcf","filtered.vcf")
	var snp_filtered_vcfFile=swapExt(snpFile, "vcf", "filtered.vcf")

	var snp_filtered_snpEff_annovar_vcfFile = swapExt(snpFile, "vcf", "filtered.snpEff.annovar.vcf")
	var indel_filtered_snpEff_annovar_vcfFile = swapExt(indelFile, "vcf", "filtered.snpEff.annovar.vcf")
	var combined_vcfFile="unifiedgenotyper.combined.vcf"

	val bams = QScriptUtils.createSeqFromFile(bamFile)


	val sampleNames = findSampleIDsFromBAMs(bams )

        //var cohortFile = new File("/mnt/isilon/cag/ngs/hiseq/Miami/RUN/IMF/VCF/test2/sampleName2.list")
        //add(writeList(sampleNames, cohortFile))

	if ( onlysnp ) {
	        call_genotypes_snp(bams, snpFile)
		annotate_snp(snp_filtered_vcfFile)
	} else {
	        call_genotypes_snp(bams, snpFile)
        	call_genotypes_indel(bams, indelFile, snpFile)	
		annotate_indel(indel_filtered_vcfFile)
		annotate_snp(snp_filtered_vcfFile)
	
		// combine SNP and INDEL
		if ( combineFlag ) {
                	add ( combineSNPandINDEL(snp_filtered_snpEff_annovar_vcfFile, indel_filtered_snpEff_annovar_vcfFile, combined_vcfFile) )
        	}

	}

        if ( selectsamplesFlag ) {
                for ( sam <- sampleNames ) {
                                if ( onlysnp ) {
					add ( selectsample ( snp_filtered_snpEff_annovar_vcfFile, "SID" + sam + ".snp.vcf", sam) )
				} else {
	
					if ( combineFlag ) {
                	                        add ( selectsample ( combined_vcfFile, "SID" + sam + ".vcf" , sam ) )
                        	        } else
                                	{
                                        	add ( selectsample ( indel_filtered_snpEff_annovar_vcfFile, "SID" + sam + ".indel.vcf", sam) )
                                        	add ( selectsample ( snp_filtered_snpEff_annovar_vcfFile, "SID" + sam + ".snp.vcf", sam) )
                                	}
				}
                }
        }


}


case class combineSNPandINDEL (snp: File, indel: File, outFile: File) extends CombineVariants with UnifiedGenotyperArguments {
	@Input var snpFile: File = snp
	@Input var indelFile: File = indel
	@Output var output: File = outFile 
        this.variant = Seq(snp, indel)
        this.out = outFile
	this.isIntermediate = false
}

case class selectsample (inVcf: File, outVcf: File, samplename: String) extends SelectVariants {
    this.variant = inVcf
    this.out = outVcf
    this.sample_name = Seq(samplename)
    //this.excludeNonVariants = true
    //this.excludeFiltered = true
    this.R = qscript.referenceFile
    this.isIntermediate = false
    this.analysisName = queueLogDir + outVcf + ".varannotator"
    this.jobName = queueLogDir + outVcf + ".varannotator"
	this.U = ValidationExclusion.TYPE.LENIENT_VCF_PROCESSING
  }

case class varannotator (inVcf: File, inSnpEffFile: File, outVcf: File) extends VariantAnnotator  {
    this.variant = inVcf
    this.snpEffFile = inSnpEffFile
    this.out = outVcf
    this.alwaysAppendDbsnpId = true
    this.D = dbsnp_file
    this.R = referenceFile
    this.A = Seq("SnpEff")
    this.isIntermediate = false
    this.analysisName = queueLogDir + outVcf + ".varannotator"
    this.jobName = queueLogDir + outVcf + ".varannotator"
    this.scatterCount = nContigs
  }


  // General arguments to non-GATK tools
  trait ExternalCommonArgs extends CommandLineFunction {
    this.memoryLimit = 6
    this.isIntermediate = false
    this.nCoresRequest = 2
  }

case class annovar_indel (inVCF: File, outVCF: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input var vcf = inVCF
	@Output var annovarVcf = outVCF
	def commandLine = "/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/annotation_pipeline/annotate_indel_queue.sh " + pwd + " UG " +  vcf + " " + annovarVcf 
    this.analysisName = queueLogDir + annovarVcf + ".annovar_indel"
    this.jobName = queueLogDir + annovarVcf + ".annovar_indel"
    this.memoryLimit=6
}

case class annovar_snp (inVCF: File, outVCF: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input var vcf = inVCF
        @Output var annovarVcf = outVCF
        def commandLine = "/mnt/isilon/cag/ngs/hiseq/respublica/pipeline/annotation_pipeline/annotate_snp_queue.sh " + pwd + " UG " +  vcf + " " + annovarVcf
    this.analysisName = queueLogDir + annovarVcf + ".annovar_snp"
    this.jobName = queueLogDir + annovarVcf + ".annovar_snp"
    this.memoryLimit=6
}

  case class writeList(inBams: Seq[File], outBamList: File) extends ListWriterFunction {
    this.inputFiles = inBams
    this.listFile = outBamList
    this.analysisName = queueLogDir + outBamList + ".bamList"
    this.jobName = queueLogDir + outBamList + ".bamList"
  }

case class myUGindel(inBams: Seq[File], outVcf: File, waitforVcf: File)  extends UnifiedGenotyper { 
    @Input var vcf = waitforVcf
    this.scatterCount = nContigs
    this.num_cpu_threads_per_data_thread = 2
    this.input_file = inBams
    this.out = outVcf
    this.genotype_likelihoods_model = GenotypeLikelihoodsCalculationModel.Model.INDEL
    this.memoryLimit = 14
    //this.min_mapping_quality_score = 0
}

case class myUGsnp (inBams: Seq[File], outVcf: File)  extends UnifiedGenotyper {
    this.scatterCount = nContigs
    this.num_cpu_threads_per_data_thread = 2
    this.input_file = inBams
    this.out = outVcf
    this.genotype_likelihoods_model = GenotypeLikelihoodsCalculationModel.Model.SNP
    this.memoryLimit = 14
    //this.min_mapping_quality_score = 0 //min_mapping_quality_score
}


}
