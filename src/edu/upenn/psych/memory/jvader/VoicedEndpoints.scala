package edu.upenn.psych.memory.jvader

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
	def voicedEndpoints(input: Array[Double], samplingRate: Int) = List[Tuple2[Double, Double]]{
		val windowSize: Int = Math.floor(samplingRate / 100.0).toInt
		println("windowSize: " + windowSize)
		val windowStep: Int = Math.floor(samplingRate / 100.0).toInt
		println("windowStep: " + windowStep)
		val numWindows = (Math.floor((input.length - windowSize) / windowStep)).toInt + 1
		println("numWindows: " + numWindows)
		
		val smallestSamples: Double = Math.ceil((smallestSegment * samplingRate) / 100.0)
		println("smallestSamples: " + smallestSamples)
		val smallestWindows: Double = Math.ceil(smallestSamples / windowSize)
		println("smallestWindows: " + smallestWindows)
		
		null
	}
}