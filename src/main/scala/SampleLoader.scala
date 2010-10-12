package edu.upenn.psych.memory.jvader

import java.io.File

import javax.sound.sampled._
import de.dfki.lt.signalproc.util._

class SampleLoader(val audioFile: File) {
	val sampleRate = AudioSystem.getAudioInputStream(audioFile).getFormat().getSampleRate.toInt
	
	/**
	 * Extracts all samples from the audio file.
	 * 
	 * @return The audio samples
	 */
	def loadSamples(): Array[Double] = {
		loadSamples(-1)
	}
	
	/**
	 * Extracts the provided number of samples from the audio file, starting with the first sample.
	 * 
	 * @param numSamples - The number of samples to return, or -1 for all samples
	 * @return The audio samples
	 */
	def loadSamples(numSamples: Int): Array[Double] = {
		val ais = AudioSystem.getAudioInputStream(audioFile)
		val adds = new AudioDoubleDataSource(ais)
		val available = adds.getDataLength.toInt
		require(numSamples <= available)
		val samplesToExtract = if (numSamples >= 0) numSamples else available
		val normSamples = new Array[Double](samplesToExtract)
		adds.getData(normSamples)
		val samples = new Array[Double](normSamples.length)
		val deNormConstant = ais.getFormat().getFrameSize * 8 - 1
		for (i <- 0 until normSamples.length)
			samples(i) = normSamples(i) * deNormConstant //de-normalize the data
		
		samples
	}
}

object SampleLoader {
	def apply(audioFile: File) = new SampleLoader(audioFile)
}