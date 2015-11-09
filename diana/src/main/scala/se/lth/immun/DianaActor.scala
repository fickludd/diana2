package se.lth.immun

import akka.actor._
import akka.routing.SmallestMailboxRoutingLogic
import java.net.InetSocketAddress
import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.Queue
import scala.collection.mutable.HashSet
import scala.collection.mutable.HashMap
import scala.collection.JavaConversions._
import se.lth.immun.diana.DianaLib.Assay
import se.lth.immun.diana.DianaLib.AssayTrace
import se.lth.immun.diana.DianaLib
import se.lth.immun.diana.DianaAnalysis.Result
import se.lth.immun.protocol.MSDataProtocol
import se.lth.immun.protocol.MSDataProtocolActors
import akka.routing.SmallestMailboxPool

object DianaActor {
	case class Analyze(assays:Seq[Assay], address:InetSocketAddress)
	case class Done(results:Seq[DianaAnalysisActor.AnalysisComplete])
	
	case class ProcessStats(
			assayTodo:Int,
			tracingPending:Int,
			traceTodo:Int,
			analysisPending:Int,
			result:Int,
			failed:Int
		) {
		override def toString = 
			"%5d %3d %3d %3d %5d %5d".format(assayTodo, tracingPending, traceTodo, analysisPending, result, failed)
			
		def done = result + failed
		def total =	assayTodo + tracingPending + traceTodo + analysisPending + result + failed
	}
}

class DianaActor(params:DianaParams) extends Actor {

	import DianaActor._
	import DianaAnalysisActor._
	import MSDataProtocol._
	import MSDataProtocolActors._
	
	var t0 = 0L
	var id = 0
	var customer:ActorRef = _
	var pantherConnection:ActorRef = _
	val analysisRouter:ActorRef = 
		context.actorOf(
			SmallestMailboxPool(params.concurrency)
				.props(DianaAnalysisActor.props(params.adv))
		)
		
	val assayTodo = new Queue[Assay]
	val tracingPending = new HashMap[Int, Assay]
	val traceTodo = new Queue[AssayTrace]
	val analysisPending = new HashMap[AssayTrace, Long]
	val results = new ArrayBuffer[AnalysisComplete]
	val failed = new ArrayBuffer[Assay]
	
	var headerPrinted = false
	
	def reportError(str:String) = {
		println(str)
	}
	
	def receive = {
		case x =>
			//println("DIANA-ACTOR got "+x.toString.take(20))
			x match {
				case Analyze(assays, address) =>
					t0 = System.currentTimeMillis
					assayTodo ++= assays
					context.actorOf(ClientInitiator.props(address, self))
					customer = sender
					
				case MSDataProtocolConnected(remote, local) =>
					pantherConnection = sender
					processAssays
					
				case MSDataReply(msg, nBytes, checkSum, timeTaken, remote) =>
					//customer ! processStats + " "
					//customer ! "DIANA| parsed %d bytes in %d ms. CHECKSUM=%d\n".format(nBytes, timeTaken, checkSum)
					//customer ! "First frag prec mz:"+msg.getTraces.getFragmentList.head.getFragment.getPrecursor
					
					val id = msg.getId
					if (tracingPending.contains(id)) {
						traceTodo += makeAssayTrace(tracingPending(id), msg)
						tracingPending -= id
					} else
						reportError("Got MSDataReply with id that was not pending. Mysterious!")
					
					process
					
				case AnalysisComplete(at, atResults) =>
					//customer ! processStats + " "
					
					if (!headerPrinted) {
						customer ! "assayTodo tracePending traceTodo analysisPending succeeded failed %s nChannels assayId\n"
										.format(atResults.timings.map(_._1).mkString(" "))
						headerPrinted = true
					}
					
					analysisPending.get(at) match {
						case Some(localT0) =>
							analysisPending -= at
							results += AnalysisComplete(at, atResults)
							
							val stats = processStats
							if (stats.done % params.verboseFreq == 0)
								customer ! "%s %s %5ds %d %s\n".format(
										stats, 
										atResults.timings.map(t => "%4dms".format(t._2)).mkString(" "), 
										(System.currentTimeMillis - t0) / 1000,
										at.assay.ms1Channels.length + at.assay.ms2Channels.length,
										at.assay.id
									)
							
						case None =>
							failed += at.assay
							reportError("Got analysis results for assay trace that's not pending. Peculiar!")
					}
					
					process
					
				case AnalysisError(at, e) =>
					reportError("Something went wrong with "+at)
					reportError(e.getMessage)
					
					analysisPending -= at
					failed += at.assay
					
					processAssays
			}
	}
	
	
	def processStats = 
		ProcessStats(assayTodo.size, tracingPending.size, traceTodo.size, analysisPending.size, results.size, failed.size)
	
	
	
	def process = {
		processAssays
		processTraces
	}
		
		
		
	def processAssays = {
		if (traceTodo.size < params.concurrency && tracingPending.isEmpty && assayTodo.nonEmpty) {
			val assay = assayTodo.dequeue
			val req = makeRequest(assay)
			tracingPending += req.getId -> assay
			//println("DIANA-ACTOR sends for tracing: "+req.toString.take(11))
			pantherConnection ! req
		}
	}
	
	
	
	def processTraces = {
		val nPending = analysisPending.size
		for (trace <- dequeue(3 * params.concurrency - nPending, traceTodo)) {
			analysisPending += trace -> System.currentTimeMillis
			//println("DIANA-ACTOR sends for analysis: "+trace.toString.take(10))
			analysisRouter ! AnalyzeAssayTrace(trace)
		}
		if (analysisPending.isEmpty && tracingPending.isEmpty && assayTodo.isEmpty && traceTodo.isEmpty) {
			//println("DIANA-ACTOR sends for done: "+results.length)
			customer ! Done(results)
		}
	}
	
	
	
	
	def makeRequest(assay:Assay) = {
		val req = GetTracesFor.newBuilder
		def makeBounds(mz:Double) = {
			 val diff = mz * params.extractWidthPPM / 1e6
			 Bounds.newBuilder.setLmz(mz - diff).setHmz(mz + diff)
		}
		
		for (channel <- assay.ms1Channels)
			req.addPrecursor(makeBounds(channel.mz))
			
		val precursorMz = assay.mz
		val precursorBounds = makeBounds(precursorMz)
		for (channel <- assay.ms2Channels)
			req.addFragment(
				FragmentBounds.newBuilder
					.setPrecursor(precursorBounds)
					.setFragment(makeBounds(channel.mz))
			)
		
		req.setTransferMode(TransferMode.DOUBLE)
		id += 1
		MasterRequest.newBuilder.setGetTracesFor(req).setId(id).build
	}
	
	
	
	def makeAssayTrace(assay:Assay, msg:MasterReply) = {
		def within(mz:Double, bounds:Bounds) =
			bounds.getLmz < mz && mz < bounds.getHmz
		
		val precTraces =
			for (precTrace <- msg.getTraces.getPrecursorList) yield
				assay.ms1Channels.find(ch => within(ch.mz, precTrace.getPrecursor)) match {
					case Some(channel) =>
						val intensity = precTrace.getTrace.getIntensityList().map(_.toDouble)
						DianaLib.Trace(channel, intensity.toArray)
					case None =>
						throw new Exception("MSDataReply precursor trace mz doesn't match assay!")
				}
		
		val fragTraces =
			for (fragTrace <- msg.getTraces.getFragmentList) yield
				assay.ms2Channels.find(ch => within(ch.mz, fragTrace.getFragment.getFragment)) match {
					case Some(channel) =>
						val intensity = fragTrace.getTrace.getIntensityList().map(_.toDouble)
						DianaLib.Trace(channel, intensity.toArray)
					case None =>
						throw new Exception("MSDataReply fragment trace mz doesn't match assay!")
				}
		
		AssayTrace(
				assay, 
				msg.getTraces.getFragment(0).getTrace.getTimeList.map(_.toDouble).toArray, 
				precTraces, fragTraces)
	}
	
	
	
	def dequeue[T](n:Int, q:Queue[T]) = {
		val k = Math.min(q.size, n)
		for (i <- 0 until k) yield q.dequeue
	}
}