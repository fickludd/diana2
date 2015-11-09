package se.lth.immun

import akka.actor._

import se.lth.immun.diana.DianaLibParams
import se.lth.immun.diana.DianaLib.AssayTrace
import se.lth.immun.diana.DianaAnalysis
import se.lth.immun.diana.DianaAnalysis.AssayResults

object DianaAnalysisActor {
	case class AnalyzeAssayTrace(at:AssayTrace)
	case class AnalysisComplete(at:AssayTrace, result:AssayResults)
	case class AnalysisError(at:AssayTrace, e:Exception)
	
	def props(params:DianaLibParams) = 
		Props(new DianaAnalysisActor(params))
}

class DianaAnalysisActor(params:DianaLibParams) extends Actor {

	import DianaAnalysisActor._
	
	def receive = {
		case AnalyzeAssayTrace(at) =>
			try {
				sender ! AnalysisComplete(at, DianaAnalysis.run(at, params))
			} catch {
				case e:Exception =>
					sender ! AnalysisError(at, e)
			}
	}
}