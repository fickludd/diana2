package se.lth.immun

import se.jt.CLIApp
import java.io.File
import java.util.Properties

object Diana extends CLIApp {

	var properties = new Properties
	properties.load(this.getClass.getResourceAsStream("/pom.properties"))
	val name 		= properties.getProperty("pom.artifactId")
	val version 	= properties.getProperty("pom.version")
	val params 		= new DianaParams(name, version)
	

	def main(args:Array[String]):Unit = {
		
		params.t0 = System.currentTimeMillis
		
		failOnError(parseArgs(name, version, args, params, List("mzML"), None))
		
		
		val (outDir, outName) = params.outBase
		
		println(name + " "+version)
    	println(" traML file: " + params.traML.value)
    	println("    out dir: " + outDir)
    	println("   out name: " + outName)
		println()
		
		val streamer 	= ReportStreamer(params, outQcZipFile(params.outBase))
		
		println(" parsing traML...")
		
		val assays = TraML.parse(new File(params.traML))
		
		println("   parsed %d assays".format(assays.length))
		
		val scorer = new DiaAssayScorer(params)
		val results = scorer.analyze(assays)
		
		Csv.write(params.outCsv, results)
	}
	
	def toOutFile(ext:String)(base:(String, String)):File = 
		new File(base._1, base._2+"."+ext)
	
	def outQcZipFile = toOutFile("qc.zip") _
}