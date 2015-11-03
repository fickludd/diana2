package se.lth.immun.diana

import se.lth.immun.markov.MarkovChainDistribution

object DianaLib {

	trait RatioIndexHolder {
		def n:Int
		def indexByChannels(ui:Int, di:Int) = {
			val k = n - ui
			val Rn = (n*(n-1))/2
			(di - ui) + Rn - (k*(k-1))/2 - 1
		}
	}
	
	trait StatsProbMode
	case object Rank extends StatsProbMode
	case object Markov extends StatsProbMode
	
	trait StatsNullMode
	case object Full extends StatsNullMode
	case object PCs extends StatsNullMode
	
	case class Channel(mz:Double, z:Int, id:String, msLevel:Int, expIntensity:Double)
	case class Trace(channel:Channel, intensity:Array[Double]) {
		def rawMap[D](f:Array[Double] => Array[Double]) = Trace(channel, f(intensity))
		def map[D](f:Array[Double] => D) = TraceData(this, f(intensity))
	} 
	case class TraceData[T](trace:Trace, data:T) {
		def map[D](f:T => D) = TraceData(trace, f(data))
	}
	case class RatioTrace(
			channel1:Channel, 
			channel2:Channel, 
			expRatio:Double, 
			ratios:Array[Double]
	) {
		def map(f:Array[Double] => Array[Double]) =
			RatioTrace(channel1, channel2, expRatio, f(ratios))
	}
	class RatioTraceGroup(
			val nChannels:Int,
			val traces:Seq[RatioTrace]
	) extends RatioIndexHolder {
		def n = nChannels
		
		def byChannels(ui:Int, di:Int) =
			traces(indexByChannels(ui, di))
		
		def map(f:Array[Double] => Array[Double]) =
			new RatioTraceGroup(nChannels, traces.map(_.map(f)))
	}
	
	case class Assay(
			id:String,
			protein:String,
			expectedRt:Double,
			mz:Double,
			z:Int,
			ms1Channels:Seq[Channel],
			ms2Channels:Seq[Channel]
	) {
		override def equals(o:Any) = 
			o match {
				case that: Assay => that.id == this.id
				case that: AssayTrace => that.assay.id == this.id
			    case _ => false
			}
		override def hashCode = id.hashCode
	}
			
	case class AssayTrace(
			assay:Assay,
			times:Array[Double],
			ms1Traces:Seq[Trace],
			ms2Traces:Seq[Trace]
	) {
		override def equals(o:Any) = 
			o match {
				case that: Assay => that.id == this.assay.id
				case that: AssayTrace => that.assay.id == this.assay.id
			    case _ => false
			}
		override def hashCode = assay.id.hashCode
	}
			
	
	def computeRatioTraces(traces:Seq[Trace]) = 
		new RatioTraceGroup(
				traces.length,
				for {
					i1 <- 0 until traces.length
					i2 <- (i1+1) until traces.length
				} yield RatioTrace(
						traces(i1).channel, 
						traces(i2).channel, 
						traces(i1).channel.expIntensity / traces(i2).channel.expIntensity,
						divide(traces(i1), traces(i2))
					)
				)
		
	def divide(t1:Trace, t2:Trace) = {
		val a = new Array[Double](math.min(t1.intensity.length, t2.intensity.length))
		for (i <- 0 until a.length) 
			a(i) = t1.intensity(i) / t2.intensity(i)
		a
	}
	
	def zipTraceData[T](traces:Seq[Trace], data:Seq[T]):Seq[TraceData[T]] =
		for (i <- 0 until math.min(traces.length, data.length)) yield TraceData(traces(i), data(i))
	
	def zipChannelData(channels:Seq[Channel], data:Seq[Array[Double]]):Seq[Trace] =
		for (i <- 0 until math.min(channels.length, data.length)) yield Trace(channels(i), data(i))
	
	/**
	 * Approximates the covariance between two stochastic variables
	 * Xi and Xj using their correlation pij and their variance v
	 */
	def estimateCov(pij:Double, v:Double):Double = 
		return pij * (3.263 + pij * (0.710 + pij*0.027)) + 
				(0.727 + pij*(0.327 - pij*(0.768 - pij*0.331)))/v
	
	
	/**
	 * Calculates the fishers statistic of a number of p-values:
	 * F = -2 * sum( ln(p) )
	 */
	def fishers(pvalues:Seq[Double]):Double = 
		return -2 * pvalues.map(p => math.log(p)).sum
		
	
	
	
	
	def getPValueTable(
			ratioTraces:Seq[RatioTrace],
			upperBound:Double
	):Seq[Double] = {
		val ub = upperBound
		val lb = 1 / upperBound
		for (r <- ratioTraces) yield 
			r.ratios.count(d => 	
							d > 0.0
						&&	d < Double.PositiveInfinity
						&&	d < r.expRatio * ub
	            		&&	d > r.expRatio * lb
	            ) / r.ratios.length.toDouble
	}
	
	
	
	def markovPValues(r:Seq[Double], target:Double, upperBound:Double) = {
		val ub = upperBound
		val lb = 1 / upperBound
		var counts = Array(Array(0, 0), Array(0, 0))
		for (i <- 1 until r.length) {
			val r0 = r(i-1)
			val r1 = r(i)
			var t0 = if (	r0 > 0.0
					&&	r0 < Double.PositiveInfinity
					&&	r0 < target * ub
					&&	r0 > target * lb) 1 else 0
			var t1 = if (	r1 > 0.0
					&&	r1 < Double.PositiveInfinity
					&&	r1 < target * ub
					&&	r1 > target * lb) 1 else 0
			counts(t0)(t1) += 1
		}
		println("    \t   TO")
		println("FROM\tok\tnot")
		println("not\t%d\t%d".format(counts(0)(1), counts(0)(0)))
		println("ok\t%d\t%d".format(counts(1)(1), counts(1)(0)))
	}
	
	
	
	def calculateMarkovDists(
			ratioGroup:RatioTraceGroup,
			n:Int,
			upperBound:Double
	):Seq[MarkovChainDistribution] = 
		for (r <- ratioGroup.traces) yield
	        calculateMarkovDist(r.ratios, r.expRatio, n, upperBound)
		
	
	
	
	def calculateMarkovDist(
			r:Seq[Double], 
			target:Double, 
			n:Int,
			upperBound:Double
	):MarkovChainDistribution = {
		val ub = upperBound
		val lb = 1 / upperBound
		var counts = Array(Array(0, 0), Array(0, 0))
		for (i <- 1 until r.length) {
			val r0 = r(i-1)
			val r1 = r(i)
			var t0 = if (	r0 > 0.0
					&&	r0 < Double.PositiveInfinity
					&&	r0 < target * ub
					&&	r0 > target * lb) 1 else 0
			var t1 = if (	r1 > 0.0
					&&	r1 < Double.PositiveInfinity
					&&	r1 < target * ub
					&&	r1 > target * lb) 1 else 0
			counts(t0)(t1) += 1
		}
		val sums = counts.map(_.sum)
		val p1 = (counts(0)(1)+counts(1)(1)) / sums.sum.toDouble
		val p11 = if (sums(1) == 0) 0.0 else counts(1)(1) / sums(1).toDouble
		val p22 = if (sums(0) == 0) 0.0 else counts(0)(0) / sums(0).toDouble
		return new MarkovChainDistribution(n, p1, p11, p22)
	}
}