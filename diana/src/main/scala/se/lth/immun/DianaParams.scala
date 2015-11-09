package se.lth.immun

import se.jt.Params
import java.io.File
import se.lth.immun.diana.DianaLibParams

class DianaParams(val name:String, val version:String) extends Params {

	import Params._
	
	val ip 		= ReqString("ip-address:port of panther server")
	val traML 	= ReqString("TraML file with assays to analyze")
	
	val extractWidthPPM = 20.0		## "PPM with of chromatogram extraction"
	
	val verbose 		= false 	## "increase details in output"
	val verboseFreq 	= 200		## "output stats every n:th assay"
	val concurrency 	= 1 		## "the number of assays to analyze in parallel"
	val profiling 		= false 	## "set to enable CPU profiling"
	
	val outDir			= ""		## "output directory (by default same as input mzML)"
	val outName			= ""		## "basename for output files (by default same as input mzML)"
	
	val nReport			= 10		## "number of random assay to export control figure for"
	val zipQcFolder = false		## "set to zip the entire qc folder on algorithm completion"
	val reportSeed		= -1L		## "seed to use for report assay selection (<0 means random)"
	
	
	val adv = new DianaLibParams
	
	def outBase = {
		val traMLFile = new File(traML)
		val dir = 
			if (outDir.value != "") outDir.value
			else traMLFile.getParent
		val name =
			if (outName.value != "") outName.value
			else stripExt(traMLFile.getName)
		(dir, name) 
	}
	
	
	def stripExt(path:String) =
		if (path.toLowerCase.endsWith(".traml"))
			path.dropRight(6)
		else
			path
	
			
	def outCsv = {
		val (dir, name) = outBase
		new File(dir, name+".csv")
		}
			
	var t0 = 0L
}