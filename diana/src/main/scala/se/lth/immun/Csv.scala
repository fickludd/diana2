package se.lth.immun

import java.io.File
import java.io.FileWriter
import java.io.BufferedWriter
import java.util.Locale

import se.lth.immun.DianaActor.AnalysisCompleteSmall
import se.lth.immun.DianaActor.ResultSmall
import se.lth.immun.DianaActor.MSLevelResultsSmall

object Csv {

	val MISSING_STRING = "_"
	
	def write(f:File, completes:Seq[AnalysisCompleteSmall]) = {
		
		val w = new BufferedWriter(new FileWriter(f))
		
		def row(vals:Seq[Any]) =
			w.write(vals.map(_.toString).mkString("\t")+"\n")
			
		row(Array(
				"rawArea",
				"correctedArea",
				"isotopeArea",
				"fragmentRankAllRatioProb","fragmentRankPcsRatioProb","fragmentMarkovAllRatioProb","fragmentMarkovPcsRatioProb","fragmentCorrScore",
				"isotopeRankAllRatioProb","isotopeRankPcsRatioProb","isotopeMarkovAllRatioProb","isotopeMarkovPcsRatioProb","isotopeCorrScore",
				"rtStart","rtApex","rtEnd",
				"rtApexAssay",
				"q1","charge",
				"estimateApexIntensity",
				"alternatives",
				"assayID","protein"))
				
		val uk = Locale.UK
		def format(str:String, a:Any) =
			str.formatLocal(uk, a)
		def formatPrecursor(str:String, xOpt:Option[MSLevelResultsSmall], f:MSLevelResultsSmall => Any) =
			xOpt.map(x => format(str, f(x))).getOrElse(MISSING_STRING)
		for {
			AnalysisCompleteSmall(assay, results) <- completes
			ResultSmall(fragment, precursor, peakTimes) <- results
		} {
			row(List(
					format("%10.1f", fragment.estimates.rawArea), 
					format("%10.1f", fragment.estimates.correctedArea), 
					formatPrecursor("%10.1f", precursor, _.estimates.rawArea), 
					format("%.2e", fragment.ratioProbs.rankAll), 
					format("%.2e", fragment.ratioProbs.rankPCs), 
					format("%.2e", fragment.ratioProbs.markovAll), 
					format("%.2e", fragment.ratioProbs.markovPCs), 
					format("%.4f", fragment.corrScore), 
					formatPrecursor("%.2e", precursor, _.ratioProbs.rankAll),
					formatPrecursor("%.2e", precursor, _.ratioProbs.rankPCs),
					formatPrecursor("%.2e", precursor, _.ratioProbs.markovAll),
					formatPrecursor("%.2e", precursor, _.ratioProbs.markovPCs),
					formatPrecursor("%.4f", precursor, _.corrScore), 
					format("%.3f", peakTimes.start), 
					format("%.3f", peakTimes.apex), 
					format("%.3f", peakTimes.end), 
					format("%5.1f", assay.expectedRt),
					format("%.2f", assay.mz),
					assay.z,
					format("%.1f", fragment.estimates.estimateApex),
					results.length,
					assay.id, 
					assay.protein))
		}
		
		w.close()
	}
}