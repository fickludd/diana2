package se.lth.immun.diana

import se.jt.Params

class DianaLibParams extends Params {

	import Params._
	
	val ratioUppedBound = 1.5	## "The maximal expected ratio deviation to accept"
	val binSize = 20			## "Bin size used for signal smoothing"
	val singleFragmentMinArea = 25.0	## "Minimum area required of a fragment in a peak candidate with only one fragment"
	val minRatioValidity = 2	## "The minimum number of supporting fragments to propagate validation"
	val pCutoff = 0.99 			## "p-value cutoff for peak candidates (default: 0.99)"
	
	lazy val getSignalProcessor = DianaSignalProcessor.getDefault
		
}