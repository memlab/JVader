package edu.upenn.psych.memory.jvader

import Math.{floor, ceil, max, min}

import scala.collection.mutable.ListBuffer

/** * 
 * Voice audio detection algorithm similar to the one described in Rabiner and Sambur (1975).
 * We've refined things to detect multiple utterances.
 * We're assuming the SNR ratio is pretty good, which should be true in the testing room.
 * 
 * [Port of pywr to the JVM. Documentation is from the original.]
 * 
 * @author original Computational Memory Lab authors, port by Yuvi Masory
 */
object VoicedEndpoints {

	/**
	 * The length of the smallest utterance to keep, in milliseconds.
	 * Helps eliminate short-lived popping and lip smacking.
	 */
	val smallestSegment = 100
	
	val IzctConstant = 0
	
	/**
	 * Detect human speech in provided audio samples.
	 * 
	 * @param input - Audio samples
	 * @samplingRate - The audio sampling rate of <tt>input</tt>
	 * 
	 * @return A <tt>List</tt> of tuples of the form (voiceStarts, voiceEnds), in order, in milliseconds
	 */
	def voicedEndpoints(input: Array[Double], samplingRate: Int): List[(Int, Int)] = {
		//~10ms window size
		val windowSize: Int = floor(samplingRate / 100.0).toInt

		//~10ms window step
		val windowStep: Int = floor(samplingRate / 100.0).toInt		
		
		//smallest segment in samples
		val smallestSamples: Double = ceil((smallestSegment * samplingRate) / 1000.0)
		
		//smallest segment size in windows
		val smallestWindows: Double = ceil(smallestSamples / windowSize)

		//move the window across the input, count the number of zero crossings and the energy of each segment
		val numWindows: Int = (floor((input.length - windowSize) / windowStep)).toInt + 1
		val zeroXings = new Array[Double](numWindows)
		val energy = new Array[Double](numWindows)
		for (i <- 0 until numWindows) {
			val window = input.slice(i * windowStep, i * windowStep + windowSize)
			val correctedWindow = window.filter(_ != 0)
			for (j <- 0 until correctedWindow.length)				
				correctedWindow(j) = correctedWindow(j).abs
			energy(i) = correctedWindow.sum
		}

		//assume first ~100ms is silence
		val numSilentWindows = (math.floor(((samplingRate / 10.0) - windowSize) / windowStep) + 1).toInt
		val numZcWindows = 0
		
		//zero-crossing threshold. if the number of zero crossings in a window is more than 2 standard deviations
		//away from the number in a window of silence, that window is likely not background noise
		val izct = IzctConstant
		
		//peak energy
		val imx = energy.max
		
		//mean energy during silence
		val imn = mean(energy.slice(0, numSilentWindows))
		
		//lower energy threshold.
		//the smaller of 3 percent of the peak energy above silence energy, or four times the silence energy
		val itl = min(0.03 * (imx - imn) + imn, 4 * imn)
		
		//upper energy threshold
		val itu = 5 * itl
		
		//look for a window with energy below the lower threshold
		var endpoints = new ListBuffer[(Int, Int)]
		var p1 = -1
		var n1 = -1
		var n2 = -1
		for (i <- 0 until numWindows) {
			//if we don't have a potential starting point and we found a window with energy greater than the
			//lower threshold, mark it as a potential starting point
			if (p1 < 0 && energy(i) >= itl)
				p1 = i
			
			//if we have a potential starting point and the current energy is less than the lower threshold,
			//unmark the currently marked starting point
			if (p1 >= 0 && energy(i) < itl)
				p1 = -1

			/* 
			 * if we haven't found the real starting point, but we have a potential starting point and we either
			 * found a window with energy greater than the upper threshold or the lower energy threshold has
			 * been maintained for more than twice the smallest segment, mark the "potential" starting point as a real
			 * starting point. the criterion of maintaining the lower energy threshold for a period of time without
			 * crossing the upper energy threshold is something that did not appear in Rabiner and Sambur (1975),
			 * but has been found to be useful with our sound files
			 */
			if (n1 < 0 && p1 >= 0 && (energy(i) >= itu || i - p1 >= 2 * smallestWindows)) {
				n1 = p1
				
				//as a further refinement, look at the previous numZcWindows windows. if three or more exceed
				//the zc threshold, set the starting point to the window where the zc threshold was first surpassed
				var numSurpassed = 0
				var firstWindow = -1
				for (wb <- (n1 - 1).until(max(n1 - 1 - numZcWindows, numSilentWindows - 1), -1)) {
					if (zeroXings(wb) > izct) {
						numSurpassed += 1
						firstWindow = wb
					}
				}
				if (numSurpassed >= 3)
					n1 = firstWindow
			}
			
			//if we have a real starting point, but no ending point, and the current window is below the lower
			//energy threshold, mark the previous window as the ending point.
			if (n1 > 0 && n2 < 0 && energy(i) < itl) {
				n2 = i - 1
				
				//as a further refinement, look at the next numZcWindows windows.
				//if three or more exceed the zc threshold, set the ending point to the window where the zc threshold
				//was last surpassed
				var numSurpassed = 0
				var lastWindow = -1
				for (wb <- (n2 + 1) until min(n2 + 1 + numZcWindows, numWindows)) {
					if (zeroXings(wb) > izct) {
						numSurpassed += 1
						lastWindow = wb
					}
				}
				if (numSurpassed >= 3)
					n2 = lastWindow
			}
			
			//if we have a real starting point, but no ending point, and this is the last frame, mark this as the ending
			if (n1 > 0 && n2 < 0 && i == numWindows - 1)
				n2 = i
				
			//if we found both a starting point and an ending point, add them to the list and start looking for the
			//next set
			if(n1 > 0 && n2 > 0) {
				if ((n2 - n1) >= smallestWindows) {
					val pts: (Double, Double) = 
						(ceil(n1 * windowStep + windowSize / 2.0).toInt, 
								ceil((n2 + 1) * windowStep - windowSize / 2.0).toInt)
				}
				p1 = -1
				n1 = -1
				n2 = -1
			}
		}
		
		endpoints.toList
	}
	
	private def mean(arr: Array[Double]): Double = {
		var sum = 0.0
		for (el <- arr)
			sum = sum + el
		sum/arr.length
	}
}