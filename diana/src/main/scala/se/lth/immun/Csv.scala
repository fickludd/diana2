package se.lth.immun

import java.io.File
import java.io.FileWriter
import java.io.BufferedWriter
import java.util.Locale

import se.lth.immun.DianaAnalysisActor.AnalysisComplete
import se.lth.immun.diana.DianaAnalysis.AssayResults
import se.lth.immun.diana.DianaAnalysis.MSLevelResults

object Csv {

	val MISSING_STRING = "_"
	
	def write(f:File, completes:Seq[AnalysisComplete]) = {
		
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
		def formatPrecursor(str:String, xOpt:Option[MSLevelResults], f:MSLevelResults => Any) =
			xOpt.map(x => format(str, f(x))).getOrElse(MISSING_STRING)
		for {
			AnalysisComplete(at, AssayResults(timings, results)) <- completes
			r <- results
		} {
			row(List(
					format("%10.1f", r.fragment.estimates.rawArea), 
					format("%10.1f", r.fragment.estimates.correctedArea), 
					formatPrecursor("%10.1f", r.precursor, _.estimates.rawArea), 
					format("%.2e", r.fragment.ratioProbs.rankAll), 
					format("%.2e", r.fragment.ratioProbs.rankPCs), 
					format("%.2e", r.fragment.ratioProbs.markovAll), 
					format("%.2e", r.fragment.ratioProbs.markovPCs), 
					format("%.4f", r.fragment.corrScore), 
					formatPrecursor("%.2e", r.precursor, _.ratioProbs.rankAll),
					formatPrecursor("%.2e", r.precursor, _.ratioProbs.rankPCs),
					formatPrecursor("%.2e", r.precursor, _.ratioProbs.markovAll),
					formatPrecursor("%.2e", r.precursor, _.ratioProbs.markovPCs),
					formatPrecursor("%.4f", r.precursor, _.corrScore), 
					format("%.1f", at.times(r.g.istart)), 
					format("%.1f", at.times(r.fragment.estimates.iEstimateApex)), 
					format("%.1f", at.times(r.g.iend)), 
					format("%5.1f", at.assay.expectedRt),
					format("%.2f", at.assay.mz),
					at.assay.z,
					format("%.1f", r.fragment.estimates.estimateApex),
					results.length,
					at.assay.id, 
					at.assay.protein))
		}
	}
}