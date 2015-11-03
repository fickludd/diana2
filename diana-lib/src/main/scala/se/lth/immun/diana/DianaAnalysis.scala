package se.lth.immun.diana

import se.lth.immun.signal.Filter

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
			val g:PCGroup
	) {
		val fragment = new MSLevelResults
		val precursor = new MSLevelResults
	}
	
	
	def run(at:AssayTrace, params:DianaLibParams) = {
		
		val signalProcessor = params.getSignalProcessor
		val allTraces 	= at.ms1Traces ++ at.ms2Traces
	
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
							.map(new Result(_))
							
		val fragmentResults = calculateFragmentScores(
				smooths.filter(_.trace.channel.msLevel == 2), 
				results,
				params
			)
			
		val ms1Smooths = smooths.filter(_.trace.channel.msLevel == 1)
		if (ms1Smooths.length > 1) 
			calculatePrecursorScores(ms1Smooths, fragmentResults, params)
		else
			fragmentResults
	}
	
	
	def calculateFragmentScores(
			smooths:Seq[TraceData[DianaSignalProcessor.SmoothAndBase]],
			results:Seq[Result],
			params:DianaLibParams
	) = {
		val reduced 	= Filter.baseLineReduce(smooths.map(_.trace.intensity).toArray)
        val savitzkied 	= DianaLib.zipChannelData(smooths.map(_.trace.channel), reduced).map(_.rawMap(Filter.savitzkyGolay9))
        val ratioGroup 	= computeRatioTraces(savitzkied)
        val state 		= new DianaChromatogramState(params, ratioGroup)
		
		for (r <- results) 
			r.fragment.gv = r.g.validate(2, smooths, state, params)
		
		state.calculateChromatogramStatistics(results.map(r => r.g -> r.fragment.gv), params)
		
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
		
		for (
			r <- results.filter(_.fragment.ratioProbs.isRemotelyUnlikely(params.pCutoff))
		) yield {
			val estimate = r.g.estimateAndIntegrate(state, smooths.map(_.trace.intensity), r.fragment.gv)
			r.fragment.corrScore = DianaPCEvaluator.corrScore(r.g, estimate)
			r.fragment.estimates = estimate
			r
		}
	}
	
	
	def calculatePrecursorScores(
			smooths:Seq[TraceData[DianaSignalProcessor.SmoothAndBase]],
			results:Seq[Result],
			params:DianaLibParams
	) = {
		val reduced 	= Filter.baseLineReduce(smooths.map(_.trace.intensity).toArray)
        val filtered 	= zipChannelData(smooths.map(_.trace.channel), reduced).map(_.rawMap(Filter.savitzkyGolay9))
        val ratioGroup 	= computeRatioTraces(filtered)
        val state 		= new DianaChromatogramState(params, ratioGroup)
		
		for (r <- results)
			r.precursor.gv = r.g.validate(1, smooths, state, params)
		
		state.calculateChromatogramStatistics(results.map(r => r.g -> r.precursor.gv), params)
		
		for (r <- results) {
			r.precursor.ratioProbs = 
				RatioProbs(
					DianaPCEvaluator.nullRatioProb(
						r.g, r.precursor.gv, state, state.statsAll, Rank
					),
					DianaPCEvaluator.nullRatioProb(
						r.g, r.precursor.gv, state, state.statsPcs, Rank
					),
					DianaPCEvaluator.nullRatioProb(
						r.g, r.precursor.gv, state, state.statsAll, Markov
					),
					DianaPCEvaluator.nullRatioProb(
						r.g, r.precursor.gv, state, state.statsPcs, Markov
					)
				)
			
			val estimate = r.g.estimateAndIntegrate(state, smooths.map(_.trace.intensity), r.precursor.gv)
			r.precursor.corrScore = DianaPCEvaluator.corrScore(r.g, estimate)
			r.precursor.estimates = estimate
		}
		
		results
	}
}