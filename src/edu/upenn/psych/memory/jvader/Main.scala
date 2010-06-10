package edu.upenn.psych.memory.jvader

import java.io.File

import scala.io._

/**
 * Entry point of application.
 * 
 * @author Yuvi Masory
 */
object Main {
	
	private val usage = "Usage: scala Main file1.wav [file2.wav file3.wav ...]\n\n" +
						"For each file the filename is printed followed by ordered pairs corresponding to predicted word endpoints."
				   
	def main(args: Array[String]) {
		if (args.length < 1) { 
			Console.err.println(usage)
			System.exit(1)
		}
		
		def filterFiles(name: String): List[File] = {
			val file = new File(name)
			if (file.exists)
				List(file)
			else {
				Console.err.println(file + " does not exist. Skipping.")
				Nil
			}
		}
		
		val files = args.flatMap(filterFiles)
		
		for (file <- files) {
		}
		
		val lines = Source.fromFile(new File("/home/yuvi/workspace/JVader/dev/pywr/input_to_voicedEndpoints.txt")).getLines("\n")
		val input: Array[Double] = 
			(for (line <- lines) yield line.toDouble).toArray
		VoicedEndpoints.voicedEndpoints(input, 11025)
	}
}