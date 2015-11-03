package se.lth.immun

import java.io.File
import java.io.FileReader
import java.io.BufferedReader
import se.lth.immun.xml.XmlReader

import se.lth.immun.traml.ghost.GhostTraML
import se.lth.immun.diana.DianaLib

object TraML {
	
	import DianaLib._

	def parse(f:File) = {
		val r = new XmlReader(new BufferedReader(new FileReader(f)))
		val gt = GhostTraML.fromFile(r)
		
		for ((pk, transitions) <- gt.transitionGroups.toSeq) yield {
			val firstTrans = transitions.head
			Assay(
				pk.toString,
				gt.peptides(pk.pepCompId).proteins.mkString(";"),
				firstTrans.rt.map(_.t).getOrElse(0.0),
				pk.mz,
				firstTrans.q1z,
				gt.includeGroups.getOrElse(pk, Nil).map(gTarget =>
					Channel(gTarget.q1, gTarget.q1z, gTarget.id, 1, gTarget.intensity.getOrElse(
						throw new Exception("[%s] Couldn't parse target intensity".format(gTarget.id))
					))
				),
				transitions.map(gTrans =>
					Channel(gTrans.q3, gTrans.q3z, gTrans.id, 2, gTrans.intensity.getOrElse(
						throw new Exception("[%s] Couldn't parse transition intensity".format(gTrans.id))
					))
				)
			)
		}
	}
}