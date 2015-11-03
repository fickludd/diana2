package se.lth.immun.diana

import DianaPeakCandidate._
import DianaSignalProcessor._
import se.lth.immun.math.Ratios

import DianaLib._

import se.lth.immun.math.Matrix
import se.lth.immun.math.Stats
import se.lth.immun.signal.Filter
import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.Queue

import org.apache.commons.math3.distribution.ChiSquaredDistribution
import org.apache.commons.math3.stat.StatUtils

object DianaPCEvaluator {
	def findEdges(
			pcs:Seq[PC]
	):Array[List[Int]] = {
		val L 		= pcs.length
		val edges 	= new Array[List[Int]](L)
		
		for (i <- 0 until L) {
            edges(i) = Nil
            val li = pcs(i)
            var k = i + 1
            var ki = li
            while ({
            	if (k < L) {
            		ki = pcs(k)
            		math.abs(li.iapex - ki.iapex) < 2
            	} else false
            }) {
                edges(i) = edges(i) :+ k
                k += 1
            }
        }
		return edges
	}
	
	
	def nullRatioProb(
			group:PCGroup,
			gv:GroupValidation,
			state:DianaChromatogramState,
			stats:DianaChromatogramState.Stats,
			mode:StatsProbMode
	):Double = {
		
		val nChannels	= state.ratioGroup.nChannels
		val nRatios		= state.ratioGroup.traces.length
		var covSum		= 0.0
		
		val pvalues = 
			for (ri <- 0 until nRatios) yield {
				val p = 
					mode match {
						case Rank =>
							val quota = group.getQuota(ri, gv)
							state.ratioOkProb(quota, ri)
						
						case Markov =>
							val ok = gv.rValids(ri).ok
							val count = 
								math.max(
									0, 
									ok.slice(group.nistart, group.niend).count(b => b)-1
								)
							val pdf = stats.markovDists(ri).pdf(group.niwidth)
							pdf.take(pdf.length - count).sum
					}
				
				if (p == 0) java.lang.Double.MIN_NORMAL else p
			}
		
		Ratios.iterate(stats.correlations.length, (ri, ui, di) => 
				covSum += DianaLib.estimateCov(stats.correlations(ui)(di), stats.variance)
			)
		
		val E = 2.0 * nRatios
		val v = 4.0 * nRatios + 2.0 * covSum
		val f = (2.0 * E * E) / v
		val c = v / (2.0 * E)
		if (f <= 0.0 || java.lang.Double.isNaN(f) || java.lang.Double.isNaN(c))
			return 1.0
		
		/*
		println("pvalues: "+pvalues.mkString(" "))
		println("corr: "+corr.mkString(" "))
		println("covs: "+corr.map(c => estimateCov(c, variance)).mkString(" "))
		println("E: "+E)
		println("v: "+v)
		println("f: "+f)
		println("c: "+c)
		println("fishers: "+fishers(pvalues))
		*/
		
		val chi = new ChiSquaredDistribution(f)
		return 1 - chi.cumulativeProbability(DianaLib.fishers(pvalues) / c)
		//group.rtProb 	= state.rtProb((group.iend + group.istart) / 2)
		//return group
	}
	
	
	
	
	/*
	def getFragmentPvalues():Seq[Double] = {
		var ps = new ArrayBuffer[Double]
		for (di <- 0 until state.numChromatograms)
			ps += 1 / state.targetRatioTable(0)(di).mean
		var sum = ps.sum
		return ps.map(_ / sum)
	}
	
	
	def pvalueForFragments(g:PCGroup, intensities:Seq[Array[Double]]):PCGroup = {
		val p0s	= getFragmentPvalues
		val xs	= intensities.map(_.slice(g.istart, g.iend).sum / 10)
		val n 	= xs.sum
		
		val s = new ArrayBuffer[Double]
		g.fragmentPvalues = p0s.zip(xs).map(t => {
			val p0 = t._1
			val x = t._2
			val my 		= n * p0
			val sigma 	= math.max(math.sqrt(n*p0*(1-p0)), 0.05 * my)
			s += sigma
			math.max(0.001, math.min(0.999, NormalDistribution.cdf((x - my) / sigma)))
		})
		
		g.sigmas = s
		g.xs = xs
		
		val chi = new ChiSquaredDistribution(2*state.numChromatograms)
		g.pEstimates = 1 - chi.cumulativeProbability(DianaUtil.fishers(g.fragmentPvalues))
		
		return g
	}
	*/
	
	
	def corrScore(g:PCGroup, est:GroupEstimation):Double = {
		val xs = est.estimates.map(Filter.savitzkyGolay9 _)//intensities.map(_.slice(g.istart, g.iend))
		
		val corrs = 
			for {
				ui <- 0 until xs.length
				di <- (ui+1) until xs.length
			} yield Stats.pearsonCorrelation(xs(ui), xs(di))
		
		//g.corrs = corrs
		val corrScore = StatUtils.mean(corrs.toArray)
		if (java.lang.Double.isNaN(corrScore)) 0.0 else corrScore
	}
}
