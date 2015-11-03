package se.lth.immun

import java.io.File
import java.io.FileWriter
import java.io.BufferedWriter
import java.util.Locale

import se.lth.immun.DianaAnalysisActor.AnalysisComplete

object Csv {

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
		for {
			AnalysisComplete(at, results) <- completes
			r <- results
		} {
			row(List(
					"%10.1f".formatLocal(uk, r.fragment.estimates.rawArea), 
					"%10.1f".formatLocal(uk, r.fragment.estimates.correctedArea), 
					"%10.1f".formatLocal(uk, r.precursor.estimates.rawArea), 
					"%.2e".formatLocal(uk, r.fragment.ratioProbs.rankAll), 
					"%.2e".formatLocal(uk, r.fragment.ratioProbs.rankPCs), 
					"%.2e".formatLocal(uk, r.fragment.ratioProbs.markovAll), 
					"%.2e".formatLocal(uk, r.fragment.ratioProbs.markovPCs), 
					"%.4f".formatLocal(uk, r.fragment.corrScore), 
					"%.2e".formatLocal(uk, r.precursor.ratioProbs.rankAll),
					"%.2e".formatLocal(uk, r.precursor.ratioProbs.rankPCs),
					"%.2e".formatLocal(uk, r.precursor.ratioProbs.markovAll),
					"%.2e".formatLocal(uk, r.precursor.ratioProbs.markovPCs),
					"%.4f".formatLocal(uk, r.precursor.corrScore), 
					"%.1f".formatLocal(uk, at.times(r.g.istart)), 
					"%.1f".formatLocal(uk, at.times(r.fragment.estimates.iEstimateApex)), 
					"%.1f".formatLocal(uk, at.times(r.g.iend)), 
					"%5.1f".formatLocal(uk, at.assay.expectedRt),
					"%.2f".formatLocal(uk, at.assay.mz),
					at.assay.z,
					"%.1f".formatLocal(uk, r.fragment.estimates.estimateApex),
					results.length,
					at.assay.id, 
					at.assay.protein))
		}
	}
}