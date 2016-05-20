package se.lth.immun

import se.jt.CLIApp
import se.lth.immun.IRT._

import java.util.Properties
import java.io.File

import scala.io.Source
import scala.util.{Either, Left, Right, Try, Failure, Success}

object DianaIRTMapper extends CLIApp {

	var properties = new Properties
	properties.load(this.getClass.getResourceAsStream("/pom.properties"))
	val name 		= properties.getProperty("pom.artifactId")
	val version 	= properties.getProperty("pom.version")
	
	val params = new DianaIRTMapperParams
	
	
	case class DataPoint(
			score:Double,
			apexRT:Double,
			assayRT:Double,
			assayID:String,
			peptide:String
		)
	
	case class StatsDataPoint(
			pValue:Double,
			pValueSum:Double,
			qValue:Double,
			dataPoint:DataPoint
		)
	
	
	val t0 = System.currentTimeMillis
	def timeStamp =
		"{%8.2fs} ".format((System.currentTimeMillis - t0) / 1000.0)
	
	def status(msg:String) = 
		println(timeStamp + msg)
		
	def info(msg:String) =
		print("".padTo(timeStamp.length, ' ') + "  " + msg)
		
	def fail(msg:String) =
		new Exception("Failure creating iRT map: "+msg)
	
		
	
	def main(args:Array[String]):Unit = {
		failOnError(parseArgs(name, version, args, params, List("target", "decoy"), None))
		
		status("reading input files...")
		val targets = readDIANAcsv(params.target)
		val decoys = readDIANAcsv(params.decoy)
		
		status("selecting best peak per assay...")
		val tBest = targets.groupBy(_.assayID).values.map(_.maxBy(_.score)).toSeq
		val dBest = decoys.groupBy(_.assayID).values.map(_.maxBy(_.score)).toSeq
		info("target assays: %d\n".format(tBest.length))
		info(" decoy assays: %d\n".format(dBest.length))
		
		status("computing qValues and applying FDR=%f cutoff...".format(params.qValue.value))
		val niceTargets = computeStats(tBest, dBest).takeWhile(_.qValue < params.qValue)
		
		val irtDataPoints = niceTargets.map(x => IRTDataPoint(x.dataPoint.peptide, x.dataPoint.assayRT, x.dataPoint.apexRT, x.dataPoint.score))

		status("building rt -> iRT map...")
		val irtMap = 
			Try(params.computeIRTMap(irtDataPoints)) match {
				case Success(map) => map
				case Failure(e) =>
					e match {
						case nde:org.apache.commons.math3.exception.NoDataException =>
							throw fail("not enough iRT datapoints, n=%d\n  IRT DATAPOINTS:\n%s".format(
									irtDataPoints.length,
									irtDataPoints.mkString("\n  ")
								))
						case e:Exception =>
							throw fail(e.getMessage)
					}
			}
		
		info("IRT MAP:\n")
		irtMap.writeMapParams(info)
		
		status("writing irtMap plot to '%s'...".format(params.outIrtMapPlot))
		irtMap.plot(params.outIrtMapPlot)
		
		for (errMsg <- irtMap.anythingWrong(params.irtR2, params.irtNAnchors))
			throw fail(errMsg)
		
		status("writing irtMap to '%s'...".format(params.outIrtMap))
		irtMap.toFile(params.outIrtMap)	
	}
	
	
	
	
	def computeStats(targets:Seq[DataPoint], decoys:Seq[DataPoint]) = {
		val tSort = targets.sortBy(-_.score)
		val dSort = decoys.sortBy(-_.score)
		
		var dIndex = 0
		var pValueSum = 0.0
		for (tIndex <- 0 until tSort.length) yield {
			val t = tSort(tIndex)
			dIndex = dSort.indexWhere(_.score < t.score, dIndex)
			val pvalue = math.max(dIndex, 0.5) / dSort.length
			pValueSum += pvalue
			val qvalue = pValueSum.toDouble / (tIndex+1)
			StatsDataPoint(pvalue, pValueSum, qvalue, t)
		}
	}
	
	
	
	def readDIANAcsv(path:String) = {
		/*rawArea correctedArea   isotopeArea     fragmentRankAllRatioProb        fragmentRankPcsRatioProb        fragmentMarkovAllRatioProb      fragmentMarkovPcsRat
ioProb      fragmentCorrScore       isotopeRankAllRatioProb isotopeRankPcsRatioProb isotopeMarkovAllRatioProb       isotopeMarkovPcsRatioProb       isotopeC
orrScore        rtStart rtApex  rtEnd   rtApexAssay     q1      charge  estimateApexIntensity   alternatives    assayID protein
		*/
		var headerParsed = false
		
		var iSCORE = -1
		var iAPEX_RT = -1
		var iASSAY_APEX_RT = -1
		var iASSAY_ID = -1
		
		val fScore = params.scoreTransform
		val ret = Seq.newBuilder[DataPoint]
		for (line <- Source.fromFile(new File(path)).getLines) {
			val p = line.split("\t").map(_.trim)
			if (!headerParsed) {
				iSCORE = p.indexOf(params.scoreCol.value)
				iAPEX_RT = p.indexOf(params.rtApexCol.value)
				iASSAY_APEX_RT = p.indexOf(params.rtAssayCol.value)
				iASSAY_ID = p.indexOf(params.assayIDCol.value)
				headerParsed = true
			} else {
				val assayID = p(iASSAY_ID) 
				ret += new DataPoint(
						fScore(p(iSCORE).toDouble),
						p(iAPEX_RT).toDouble,
						p(iASSAY_APEX_RT).toDouble,
						assayID,
						assayID.split("_").head
					)
			}
		}
		ret.result
	}
}