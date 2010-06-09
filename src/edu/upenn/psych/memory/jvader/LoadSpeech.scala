package edu.upenn.psych.memory.jvader

import java.io.File

import javax.sound._
import javax.sound.sampled._

import de.dfki.lt.signalproc.filter._
import de.dfki.lt.signalproc.util._

class LoadSpeech(audioFile: File) {
	
	def loadSpeech(): Array[Double] = {
		val ais = AudioSystem.getAudioInputStream(audioFile);
		val adds = new AudioDoubleDataSource(ais)
		val samples = new Array[Double](adds.getDataLength.asInstanceOf[Int])
		samples
	}	
}

object LoadSpeech {
	def apply(audioFile: File) = new LoadSpeech(audioFile)
}