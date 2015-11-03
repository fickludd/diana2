package se.lth.immun

import se.lth.immun.diana.DianaLib._
import se.lth.immun.diana.DianaAnalysis

import akka.actor._
import scala.concurrent.duration.Duration
import java.util.concurrent.TimeUnit
import java.net.InetSocketAddress

class DiaAssayScorer(params:DianaParams) {

	val actorSystem 	= ActorSystem("actor-system")
	val actorInbox 		= Inbox.create(actorSystem)
	
	val dianaActor		= actorSystem.actorOf(Props(new DianaActor(params)))
	
	
	def await[T](pf:PartialFunction[AnyRef, T]):T = {
		var res:Option[T] = None
		while (res.isEmpty)
			try {
				actorInbox.receive(Duration.create(1, TimeUnit.SECONDS )) match {
					case str:String => if (params.verbose) print(str)
					case a:AnyRef => 
						if (pf.isDefinedAt(a)) 
							res = Some(pf(a))
				}
			} catch {
				case e:java.util.concurrent.TimeoutException => {}
			}
		res.get
	} 
	def awaitAnalysis = await({ case DianaActor.Done(results) => results })
	
	
	def analyze(assays:Seq[Assay]):Seq[DianaAnalysisActor.AnalysisComplete] = {
		val ipRE = """(\d+\.\d+\.\d+\.\d+):(\d+)""".r
		val localhostRE = """localhost:(\d+)""".r
		val address = 
			params.ip.value match {
				case ipRE(ip, port) => new InetSocketAddress(ip, port.toInt)
				case localhostRE(port) => new InetSocketAddress("localhost", port.toInt)
		}
		
		actorInbox.send(dianaActor, DianaActor.Analyze(assays, address))
		
		val results = awaitAnalysis
		actorSystem.awaitTermination
		results
	}
}