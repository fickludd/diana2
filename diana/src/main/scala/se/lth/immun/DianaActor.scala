package se.lth.immun

import akka.actor._
import akka.routing.SmallestMailboxRoutingLogic
import java.net.InetSocketAddress
import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.Queue
import scala.collection.mutable.HashSet
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
}

class DianaActor(params:DianaParams) extends Actor {

	import DianaActor._
	import DianaAnalysisActor._
	import MSDataProtocol._
	import MSDataProtocolActors._
	
	var id = 0
	var customer:ActorRef = _
	var pantherConnection:ActorRef = _
	val analysisRouter:ActorRef = 
		context.actorOf(
			SmallestMailboxPool(params.concurrency)
				.props(DianaAnalysisActor.props(params.adv))
		)
		
	val assayTodo = new Queue[Assay]
	val tracingPending = new HashSet[(Int, Assay)]
	val traceTodo = new Queue[AssayTrace]
	val analysisPending = new HashSet[AssayTrace]
	val results = new ArrayBuffer[AnalysisComplete]
	
	def reportError(str:String) = {
		println(str)
	}
	
	def receive = {
		case x =>
			//println("DIANA-ACTOR got "+x.toString.take(100))
			x match {
				case Analyze(assays, address) =>
					assayTodo ++= assays
					context.actorOf(ClientInitiator.props(address, self))
					customer = sender
					
				case MSDataProtocolConnected(remote, local) =>
					pantherConnection = sender
					processAssays
					
				case MSDataReply(msg, nBytes, checkSum, timeTaken, remote) =>
					println("DIANA| parsed %d bytes in %d ms. CHECKSUM=%d".format(nBytes, timeTaken, checkSum))
					println("First frag prec mz:"+msg.getTraces.getFragmentList.head.getFragment.getPrecursor)
					
					tracingPending.find(_._1 == msg.getId) match {
						case Some((id, assay)) =>
							traceTodo += makeAssayTrace(assay, msg)
						case None =>
							reportError("Get MSDataReply with id that was not pending. Mysterious!")
					}
					
					processAssays
					processTraces
					
				case AnalysisComplete(at, atResults) =>
					println("DIANA| completed analysis of "+at)
					
					if (!analysisPending.contains(at))
						reportError("Got analysis results for assay trace that's not pending. Peculiar!")
						
					analysisPending -= at
					results += AnalysisComplete(at, atResults)
					
					processTraces
					
				case AnalysisError(at, e) =>
					reportError("Something went wrong with "+at)
					reportError(e.getMessage)
					
					analysisPending -= at
					
					processTraces
			}
	}
	
	
	
	def processAssays = {
		val nPending = tracingPending.size
		for (assay <- dequeue(params.concurrency - nPending, assayTodo)) {
			val req = makeRequest(assay)
			tracingPending += req.getId -> assay
			//println("DIANA-ACTOR sends: "+req)
			pantherConnection ! req
		}
	}
	
	
	
	def processTraces = {
		val nPending = analysisPending.size
		for (trace <- dequeue(3 * params.concurrency - nPending, traceTodo)) {
			analysisPending += trace
			analysisRouter ! AnalyzeAssayTrace(trace)
		}
		if (analysisPending.size == 0) {
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