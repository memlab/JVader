package edu.upenn.psych.memory.jvader

import Math.{floor, ceil}

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
	
	/**
	 * Detect human speech in provided audio samples.
	 * 
	 * @param input - Audio samples
	 * @samplingRate - The audio sampling rate of <tt>input</tt>
	 * 
	 * @return A <tt>List</tt> of tuples of the form (voiceStarts, voiceEnds), in order, in milliseconds
	 */
	def voicedEndpoints(input: Array[Double], samplingRate: Int) = List[Tuple2[Double, Double]] {
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

		
		null
	}
}