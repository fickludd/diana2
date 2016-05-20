package se.lth.immun

import java.io.File
import se.jt.Params

class DianaIRTMapperParams extends Params {
	import Params._

	val target 	= ReqString("csv-file with target DIANA results")
	val decoy 	= ReqString("csv-file with decoy DIANA results")
	
	val scoreCol = "fragmentMarkovAllRatioProb"	## "Column to interpret as score"
	val rtApexCol = "rtApex"					## "Column to interpret as the apex rt"
	val rtAssayCol = "rtApexAssay"				## "Column to interpret as the assay rt (iRT)"
	val assayIDCol = "assayID"					## "Column to interpret as the assay ID"
	val scoreType = "pvalue"					## "how to interpret scores (pvalue, positive or negative)"
	
	val qValue = 0.01				## "q-value threshold for target data points"
	
	val irtR2 = 0.9					## "r2 threshold required for regression iRT maps. Will fail if this level is not met."
	val irtMode = "weight-mean-ipolate"		## "The iRT mapping algorithm to use: simple-reg, robust-reg, median-ipolate or weight-mean-ipolate"
	val irtNAnchors = 5				## "The number of anchor points (iRT peptides) required for interpolating iRT maps. Will fail if this number is not met."
	
	val outDir			= ""		## "output directory (by default same as input mzML)"
	val outName			= ""		## "basename for output files (by default same as input mzML)"
	val verbose 		= false		## "set to enable a lot of output"
	
	def scoreTransform:Double => Double = {
		scoreType.value match {
			case "pvalue" => 1 - _
			case "positive" => x => x
			case "negative" => x => -x
		}
	}
	
	def computeIRTMap = 
		irtMode.value match {
			case "median-ipolate" => IRT.medianInterpolationMap _
			case "weight-mean-ipolate" => IRT.weightedMeanInterpolationMap _
			case "simple-reg" => IRT.simpleRegressionMap _
			case "robust-reg" => IRT.robustRegressionMap _
		}
	
	def outIrtMap = {
		val (dir, name) = outBase
		new File(dir, name + ".irtmap")
	}
	
	def outIrtMapPlot = {
		val (dir, name) = outBase
		new File(dir, name + ".irtmap.png")
	}
	
	def outBase = {
		val targetFile = new File(target)
		val dir = 
			if (outDir.value != "") outDir.value
			else targetFile.getParent
		val name =
			if (outName.value != "") stripExt(outName.value, ".irtmap")
			else stripExts(targetFile.getName) + "." + irtMode.value
		(dir, name) 
	}
	
	def stripExt(path:String, ext:String) =
		if (path.toLowerCase.endsWith(ext.toLowerCase))
			path.dropRight(ext.length)
		else path
	
	def stripExts(path:String) =
		stripExt(stripExt(path, ".csv"), ".tsv")
}