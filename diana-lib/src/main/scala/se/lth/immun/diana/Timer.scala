package se.lth.immun.diana

class Timer {

	var t0 = System.currentTimeMillis
	var tMid = t0
	
	def click = {
		val temp = tMid
		tMid = System.currentTimeMillis
		tMid - temp
	}
	
	def total = {
		System.currentTimeMillis - t0
	}
	
	def reset = {
		val tot = total
		t0 = System.currentTimeMillis
		tMid = t0
	}
	
}