package se.lth.immun.diana

import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.Queue

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation
import org.apache.commons.math3.stat.StatUtils
import org.apache.commons.math3.distribution.ChiSquaredDistribution

import se.lth.immun.math.Matrix
import se.lth.immun.markov.MarkovChainDistribution

import DianaSignalProcessor._


object DianaPeakCandidate {
	
	import DianaLib._
	
	case class PC(istart:Int, iapex:Int, iend:Int, trace:Trace)
	
	
	
	
	class Validation(
		val icurve:Int,
		width:Int
	) {
		var ok 		= new Array[Boolean](width)
		var flag 	= false
		for (i <- 0 until width) ok(i) = false
	}
	
	
	
	
	
	object GroupValidation {
		
		def apply(nValids:Int, iwidth:Int) = 
			new GroupValidation(
				for (i <- 0 until nValids) yield new Validation(i, iwidth),
				for {
					i <- 0 until nValids
					j <- (i+1) until nValids
				} yield new Validation(0, iwidth)
			)
		
	}
	
	class GroupValidation(
		val valids:Seq[Validation],
		val rValids:Seq[Validation]
	) extends DianaLib.RatioIndexHolder {
		val n = valids.length
		def rValidByChannels(ui:Int, di:Int) = 
			rValids(indexByChannels(ui, di))
			
	}
	
	
	
	
	
	object GroupEstimation {
		
		def apply(nChrom:Int, niwidth:Int) = {
			val est = new GroupEstimation
			
			est.estimateApex = 0.0
			est.estimates 	= new Array[Array[Double]](nChrom)
	        for (e <- 0 until nChrom) 
	        	est.estimates(e) = new Array[Double](niwidth)
	        
	        est.correctedAreas = new Array(nChrom)
	        est.rawAreas = new Array(nChrom)
	        est
		}
	}
	
	class GroupEstimation {
		var estimates:Array[Array[Double]] = null
		var estimateApex:Double = -1.0
		var iEstimateApex:Int = -1
		var correctedAreas:Array[Double] = null
		var rawAreas:Array[Double] = null
			
		def rawArea = if (rawAreas == null || rawAreas.isEmpty) 0.0 else rawAreas.sum
		def correctedArea = if (correctedAreas == null || correctedAreas.isEmpty) 0.0 else correctedAreas.sum
	}
	
	
	
	
	
	class PCGroup {
		var pcs 	= new ArrayBuffer[PC]
		var pvalueTable:Seq[Double] 	= Nil	// Serialized in ratio index order
		
		override def toString = istart + "-"+ iend
	
		def iend 	= StatUtils.percentile(pcs.map(_.iend.toDouble).toArray, 50).toInt //pcs.map(_.iend).max
		def istart 	= StatUtils.percentile(pcs.map(_.istart.toDouble).toArray, 50).toInt //pcs.map(_.istart).min
		def iwidth 	= iend - istart
		
		var nistart = 0
		var niend = 0
		def niwidth = niend - nistart
		
		def getQuota(ri:Int, gv:GroupValidation) = {
			val ok = gv.rValids(ri).ok
			val count = 
				math.max(
						0, 
						ok.slice(nistart, niend).count(b => b)-1
					)
			count.toDouble / niwidth
		}
		
		
		
		def validate(
				msLevel:Int,
				smooths:Seq[TraceData[SmoothAndBase]],
				state:DianaChromatogramState,
				params:DianaLibParams
		):GroupValidation = {
			//var istart 	= group.istart
			//var iend 	= group.iend 
			//var iwidth	= group.iwidth
			val channels 	= smooths.map(_.trace.channel)
			val nChannels	= channels.length
			val ub 		= state.upperBound
			val lb 		= state.lowerBound
			
			val handled 	= for (i <- 0 until state.ratioGroup.nChannels) yield false
			val gv 			= GroupValidation(nChannels, iwidth)
			
			val q = new Queue[Validation]
			for (pc <- pcs.filter(_.trace.channel.msLevel == msLevel)) {
				val iChannel = channels.indexOf(pc.trace.channel)
				q += gv.valids(iChannel)
				for (ir <- math.max(istart, pc.istart) until math.min(iend, pc.iend))
					gv.valids(iChannel).ok(ir - istart) = true
			}
			
			while (!q.isEmpty) {
				val v = q.dequeue
				for (i <- 0 until nChannels) {
					if (i == v.icurve)
						gv.rValidByChannels(i, i).flag = true
					else if (!gv.rValidByChannels(v.icurve, i).flag) {
						val ui = math.min(i, v.icurve)
		            	val di = math.max(i, v.icurve)
		            	val r		= state.ratioGroup.byChannels(ui, di)
		            	val y		= smooths(i).data.smooth
		            	val bl		= smooths(i).data.base
		            	val bL		= bl.length
		            	val rValid	= gv.rValidByChannels(ui, di)
		            	for (ir <- istart until iend) {
		            		if (		r.ratios(ir) < r.expRatio * ub
		            				&&	r.ratios(ir) > r.expRatio * lb
		            				&&	v.ok(ir - istart)
		            				&& 	y(ir) > bl(math.max(0, math.min(bL-1, ir - state.params.binSize / 2)))) {
		            			gv.valids(i).ok(ir - istart) = true
		            			rValid.ok(ir - istart) = true
		            		} else
		            			rValid.ok(ir - istart) = false
		            	}
						gv.rValidByChannels(v.icurve, i).flag = true
						gv.rValidByChannels(i, v.icurve).flag = true
						if (gv.rValidByChannels(i, i).flag)
							q += gv.valids(i)
					}
				}
			}
			
			for (i <- 0 until nChannels) {
				val counts = new Array[Int](iwidth)
				for (j <- 0 until nChannels) {
					if (i != j) {
						val ui = math.min(i, j)
			            val di = math.max(i, j)
			            val ok = gv.rValidByChannels(ui, di).ok
						for (k <- 0 until iwidth) {
				            if (ok(k))
				            	counts(k) += 1
						}
					}
				}
				gv.valids(i).ok = counts.map(_ >= params.minRatioValidity)
			}
			
			val rvs 	= gv.rValids.map(_.ok)
				
			val rvsSummed 	= rvs.transpose.map(_.count(b => b))
			val vTotal 		= rvsSummed.sum
			val temp 		= 2.0 / rvs.length
			val ratio 		= vTotal.toDouble / (rvs.length * iwidth)
			
			nistart 	= 0
			niend 		= iwidth
			while (rvsSummed(nistart) * temp < ratio) nistart += 1
			while (rvsSummed(niend-1) * temp < ratio) niend -= 1
			
			return gv
		}
		
		
		
		def estimateAndIntegrate(
				state:DianaChromatogramState, 
				y:Seq[Array[Double]],
				gv:GroupValidation
		):GroupEstimation = {
	        def stable(rt:DianaLib.RatioTrace):Boolean = true//0.5 * r.mean > r.stdDev
	
	        val L		= state.ratioGroup.nChannels
	        val t0		= istart
	        val est 	= GroupEstimation(L, niwidth)
	        
	
	        for (t <- nistart until niend) {
	            for (b <- 0 until L) {
	                try {
	                	if (gv.valids(b).ok(t))
		                    est.estimates(b)(t - nistart) = y(b)(t0 + t)
		                else {
		                    var localEstimates = List[Double]()
		                    for (e <- 0 until L) {
		                        if (
		                        			e != b 
		                        		&& 	gv.valids(e).ok(t) 
		                        		&& 	stable(state.ratioGroup.byChannels(b, e))
		                        ) {
		                        	var multiplier = state.ratioGroup.byChannels(b, e).expRatio
		                        	localEstimates = localEstimates :+ (y(e)(t0 + t) * multiplier)
		                        }
		                    }
		                    est.estimates(b)(t - nistart) = 
		                    	if (localEstimates.isEmpty)	 y(b)(t0 + t) 
		                    	else 
		                    		math.min(StatUtils.mean(localEstimates.toArray), y(b)(t0 + t))
		                }
	                } catch {
	                	case e:Exception => e.printStackTrace()
	                }
	            }
	        }
	        
	        
	        for (ifrag <- 0 until L) {
	            est.correctedAreas(ifrag) = 0
	            var vals = y(ifrag)
	            for (et <- 0 until niwidth) {
	            	val x = 
		                if (vals(t0 + nistart + et) <= est.estimates(ifrag)(et))
		                	vals(t0 + nistart + et)
						else 
							est.estimates(ifrag)(et)
					est.correctedAreas(ifrag) += x
					if (x > est.estimateApex) {
						est.estimateApex = x
						est.iEstimateApex = t0 + nistart + et
					}
	            }
	        }
	        
	        est
		}
	}
	
	
	
	def findPCs(
			y:Seq[Double],
			dy:Seq[Double],
			ddy:Seq[Double],
			baseline:Seq[Double],
			binSize:Int,
			trace:Trace
	):Seq[PC] = {
		val ret = new ArrayBuffer[PC]
		
		var istart 	= -1
		var iapex 	= -1
		for (i <- 0 until dy.length-1) {
			if (math.signum(dy(i)) != math.signum(dy(i+1))) {
				if (ddy(i) > 0) { // local min
					if (iapex >= 0 && istart >= 0) {
						val max 	= y(iapex)
						val starty 	= y(istart)
						val endy	= y(i+1)
						val bl = baseline(math.max(0, math.min(baseline.length-1, iapex - binSize / 2)))
						if (max > 2*math.max(starty, endy) && max > 2 * bl)
							ret += new PC(istart, iapex, i+1, trace)
					}
						
					istart = i+1
					iapex = -1
				} else { // local max
					iapex = i+1
				}
			}
		}
		
		return ret
	}
	
	
	
	def groupPCs(
			pcs:Seq[Seq[PC]],
			findEdges:Seq[PC] => Array[List[Int]]
	):Seq[PCGroup] = {
			
		val sorted 	= pcs.flatten.sortBy(_.iapex)
		val edges	= findEdges(sorted)
		
		val L				= sorted.length
        val groupRefList 	= new Array[PCGroup](L)
        val groupList 		= new ArrayBuffer[PCGroup]

        for (i <- 0 until L) {
            for (e <- edges(i))
                if (groupRefList(e) != null) groupRefList(i) = groupRefList(e)
            
            if (groupRefList(i) == null) {
                val pc = new PCGroup
                groupList += pc
                groupRefList(i) = pc
            }
            for (e <- edges(i))
                if (groupRefList(e) == null) groupRefList(e) = groupRefList(i)
        }

        for (i <- 0 until L)
            groupRefList(i).pcs += sorted(i)
        
        return groupList
	}
}
