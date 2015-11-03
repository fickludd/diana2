package se.lth.immun.diana

import se.lth.immun.math.Ratios
import se.lth.immun.math.Matrix
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation
import org.apache.commons.math3.stat.StatUtils
import org.apache.commons.math3.distribution.NormalDistribution
import se.lth.immun.math.Stats
import DianaPeakCandidate._
import DianaLib._
import se.lth.immun.markov.MarkovChainDistribution




object DianaChromatogramState {
	class Stats(
		val ratioGroup:RatioTraceGroup,
		pcMaxWidth:Int,
		params:DianaLibParams
	) {
		val variance = 
			ratioGroup.traces.map(tr => 
				StatUtils.variance(
					tr.ratios.filter(d => 
						d > 0.0 && d < Double.PositiveInfinity
					))
			).sum / ratioGroup.traces.length
		
		
		val correlations = {
			val rL	= ratioGroup.traces.length
			val x 	= Matrix.get2d[Double](rL)
			for {
				ui <- 0 until rL
				di <- (ui+1) until rL
			} {
				var rr = ratioGroup.traces(ui).ratios.zip(ratioGroup.traces(di).ratios)
							.filter(t => 
										t._1 > 0.0 && t._1 < Double.PositiveInfinity
									&&	t._2 > 0.0 && t._2 < Double.PositiveInfinity
								).unzip
				x(ui)(di) = Stats.pearsonCorrelation(rr._1.toArray, rr._2.toArray)
			}
			x
		}
		
		
		val markovDists = 
			DianaLib.calculateMarkovDists(ratioGroup, pcMaxWidth, params.ratioUppedBound)
	}
}




class DianaChromatogramState(
		val params:DianaLibParams,
		val ratioGroup:RatioTraceGroup
) {
	import DianaChromatogramState._
	
	val upperBound = params.ratioUppedBound.value
	val lowerBound = 1 / upperBound
	
	var npcgs = 0
	var ratioOkDists:Seq[Array[Double]] = Nil
	var statsAll:Stats = _
	var statsPcs:Stats = _
	
	
	/*
	def buildExpRatioTable(n:Int, ratioTraces:Seq[RatioTrace]) = {
		val table = new Array[Array[Double]](n)
		var ri = 0
		for (i <- 1 until n) {
			table(i) = new Array[Double](n)
			for (j <- i+1 until n) {
				table(i)(j) = ratioTraces(ri).expRatio
				ri += 1
			}
		}
		table
	}
	*/
		
		
	def ratioOkProb(qouta:Double, ri:Int) = 
		ratioOkDists(ri).count(_ >= qouta).toDouble / npcgs
	
	
	
	def calculateChromatogramStatistics(pcs:Seq[(PCGroup, GroupValidation)], params:DianaLibParams) = {
		
		val maxWidth 	= pcs.map(_._1.iwidth).max
		npcgs			= pcs.length
		ratioOkDists	= 
			for (ri <- 0 until ratioGroup.traces.length) yield
				pcs.map(t => t._1.getQuota(ri, t._2)).sorted.toArray
		
		val intervals = pcs.map(t => (t._1.istart, t._1.iend))
		val iL	= intervals.map(t => t._2 - t._1).sum
		statsPcs = new Stats(ratioGroup.map(x => slices(x, intervals, iL)), maxWidth, params)
		statsAll = new Stats(ratioGroup, maxWidth, params)
	}
	
	
	
								
	def slices(x:Array[Double], intervals:Seq[(Int, Int)], iLength:Int = -1) = {
		var iL = iLength
		if (iL < 0) iL = intervals.map(t => t._2 - t._1).sum
		val ret = new Array[Double](iL)
		var ir = 0
		for (ival <- intervals) 
			for (ix <- ival._1 until ival._2) {
				ret(ir) = x(ix)
				ir += 1
			}
		ret
	}
}
