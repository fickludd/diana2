package se.lth.immun.diana

import se.lth.immun.signal.Filter
import scala.collection.mutable.ArrayBuffer

object DianaAnalysis {

	import DianaLib._
	import DianaPeakCandidate.PCGroup
	import DianaPeakCandidate.GroupValidation
	import DianaPeakCandidate.GroupEstimation
	
	
	case class RatioProbs(
			rankAll:Double,
			rankPCs:Double,
			markovAll:Double,
			markovPCs:Double
	) {
		def isRemotelyUnlikely(cutoff:Double) = 
			rankAll < cutoff ||
			rankPCs < cutoff ||
			markovAll < cutoff ||
			markovPCs < cutoff
	}
		
	class MSLevelResults {
		var gv:GroupValidation = _
		var ratioProbs:RatioProbs = _
		var corrScore:Double = _
		var estimates:GroupEstimation = _
	}
		
	class Result(
			val g:PCGroup,
			hasMS1Data:Boolean
	) {
		val fragment = new MSLevelResults
		val precursor = if (hasMS1Data) Some(new MSLevelResults) else None
	}
	
	
	case class AnalysisTimings(peakDetection:Long, fragmentScoring:Long, precursorScoring:Long, total:Long)	
	case class AssayResults(timings:ArrayBuffer[(String, Long)], results:Seq[Result])
	
	
	def run(at:AssayTrace, params:DianaLibParams):AssayResults = {
		
		val t = new Timer
		val timings = new ArrayBuffer[(String, Long)]
		
		val signalProcessor = params.getSignalProcessor
		val allTraces 	= at.ms1Traces ++ at.ms2Traces
		val useMS1 		= at.ms1Traces.length > 1
	
		def getPCs(trace:Trace, y:DianaSignalProcessor.SmoothAndBase) = {	
			var dy 		= signalProcessor.getDerivate(y.smooth)
			var ddy		= signalProcessor.getDerivate(dy)
			
			DianaPeakCandidate.findPCs(y.smooth, dy, ddy, y.base, params.binSize, trace)
		}
		
		def nonRetarded(pcGroup:PCGroup) = {
			val nFrags = pcGroup.pcs.count(_.trace.channel.msLevel == 2)
			nFrags > 1 || 
			(pcGroup.pcs.find(_.trace.channel.msLevel == 2) match {
				case Some(pc) =>
					pc.trace.intensity.slice(pc.istart, pc.iend).sum > params.singleFragmentMinArea
				case None => false
			})
		}
		
		val smooths 	= allTraces.map(_.map(signalProcessor.getSmoothAndBase))
		val pcs			= smooths.map(traceDat => getPCs(traceDat.trace, traceDat.data))
		val results		= DianaPeakCandidate.groupPCs(pcs, DianaPCEvaluator.findEdges)
							.filter(nonRetarded)
							.map(new Result(_, useMS1))
							
		timings += "tPeakDetection" -> t.click
		if (results.isEmpty) {
			timings += "total" -> t.total
			AssayResults(timings, results)
		
		} else {
			
			val fragmentResults = calculateFragmentScores(
					smooths.filter(_.trace.channel.msLevel == 2), 
					results,
					params,
					timings
				)
				
			timings += "tFragmentScoring" -> t.click
				
			val finalResults =
				if (useMS1) {
					val ms1Smooths = smooths.filter(_.trace.channel.msLevel == 1)
					calculatePrecursorScores(ms1Smooths, fragmentResults, params, timings)
				} else
					fragmentResults
					
			timings += "tPrecursorScoring" -> t.click
			timings += "total" -> t.total
			
			AssayResults(timings, finalResults)
		}
	}
	
	
	def calculateFragmentScores(
			smooths:Seq[TraceData[DianaSignalProcessor.SmoothAndBase]],
			results:Seq[Result],
			params:DianaLibParams,
			timings:ArrayBuffer[(String, Long)]
	):Seq[Result] = {
		if (results.isEmpty) {
			timings += "tSmoothing" -> 0L
			timings += "tRatioCalcs" -> 0L 
			timings += "tValidations" -> 0L
			timings += "tStatsCalcs" -> 0L
			timings += "tRatioProbs" -> 0L
			timings += "tCorrAndEstimate" -> 0L
			return results
		}
		
		val t = new Timer
		
		val reduced 	= Filter.baseLineReduce(smooths.map(_.trace.intensity).toArray)
        val savitzkied 	= DianaLib.zipChannelData(smooths.map(_.trace.channel), reduced).map(_.rawMap(Filter.savitzkyGolay9))
        
        timings += "tSmoothing" -> t.click
        
        val ratioGroup 	= computeRatioTraces(savitzkied)
        val state 		= new DianaChromatogramState(params, ratioGroup)
		
        timings += "tRatioCalcs" -> t.click
        
		for (r <- results) 
			r.fragment.gv = r.g.validate(2, smooths, state, params)
		
        timings += "tValidations" -> t.click
        
		state.calculateChromatogramStatistics(results.map(r => r.g -> r.fragment.gv), params)
		
        timings += "tStatsCalcs" -> t.click
        
		for (r <- results) {
			r.fragment.ratioProbs = 
				RatioProbs(
					DianaPCEvaluator.nullRatioProb(
						r.g, r.fragment.gv, state, state.statsAll, Rank
					),
					DianaPCEvaluator.nullRatioProb(
						r.g, r.fragment.gv, state, state.statsPcs, Rank
					),
					DianaPCEvaluator.nullRatioProb(
						r.g, r.fragment.gv, state, state.statsAll, Markov
					),
					DianaPCEvaluator.nullRatioProb(
						r.g, r.fragment.gv, state, state.statsPcs, Markov
					)
				)
		}
        timings += "tRatioProbs" -> t.click
		
		val finalResults = 
			for (
				r <- results.filter(_.fragment.ratioProbs.isRemotelyUnlikely(params.pCutoff))
			) yield {
				val estimate = r.g.estimateAndIntegrate(state, smooths.map(_.trace.intensity), r.fragment.gv)
				r.fragment.corrScore = DianaPCEvaluator.corrScore(r.g, estimate)
				r.fragment.estimates = estimate
				r
			}
        
        timings += "tCorrAndEstimate" -> t.click
        
        finalResults
	}
	
	
	def calculatePrecursorScores(
			smooths:Seq[TraceData[DianaSignalProcessor.SmoothAndBase]],
			results:Seq[Result],
			params:DianaLibParams,
			timings:ArrayBuffer[(String, Long)]
	):Seq[Result] = {
		if (results.isEmpty) return results
		
		val reduced 	= Filter.baseLineReduce(smooths.map(_.trace.intensity).toArray)
        val filtered 	= zipChannelData(smooths.map(_.trace.channel), reduced).map(_.rawMap(Filter.savitzkyGolay9))
        val ratioGroup 	= computeRatioTraces(filtered)
        val state 		= new DianaChromatogramState(params, ratioGroup)
		
		for (r <- results)
			r.precursor.get.gv = r.g.validate(1, smooths, state, params)
		
		state.calculateChromatogramStatistics(results.map(r => r.g -> r.precursor.get.gv), params)
		
		for (r <- results) {
			r.precursor.get.ratioProbs = 
				RatioProbs(
					DianaPCEvaluator.nullRatioProb(
						r.g, r.precursor.get.gv, state, state.statsAll, Rank
					),
					DianaPCEvaluator.nullRatioProb(
						r.g, r.precursor.get.gv, state, state.statsPcs, Rank
					),
					DianaPCEvaluator.nullRatioProb(
						r.g, r.precursor.get.gv, state, state.statsAll, Markov
					),
					DianaPCEvaluator.nullRatioProb(
						r.g, r.precursor.get.gv, state, state.statsPcs, Markov
					)
				)
			
			val estimate = r.g.estimateAndIntegrate(state, smooths.map(_.trace.intensity), r.precursor.get.gv)
			r.precursor.get.corrScore = DianaPCEvaluator.corrScore(r.g, estimate)
			r.precursor.get.estimates = estimate
		}
		
		results
	}
}