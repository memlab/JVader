package edu.upenn.psych.memory.jvader

import java.io.File

import javax.sound.sampled._
import de.dfki.lt.signalproc.util._

class SampleLoader(val audioFile: File) {
	private val ais = AudioSystem.getAudioInputStream(audioFile);
	val sampleRate = ais.getFormat().getSampleRate.toInt
	
	def loadSamples(): Array[Double] = {		
		val adds = new AudioDoubleDataSource(ais)
		val nSamples = new Array[Double](adds.getDataLength.toInt)
		adds.getData(nSamples)
		val samples = new Array[Double](nSamples.length)
		val deNormConstant = ais.getFormat().getFrameSize * 8 - 1
		for (i <- 0 until nSamples.length)
			samples(i) = nSamples(i) * deNormConstant //de-normalize the data
		
		samples
	}	
}

object LoadSpeech {
	def apply(audioFile: File) = new SampleLoader(audioFile)
}