package edu.upenn.psych.memory.jvader

object Main {
	
	val usage = "Usage: scala Main file1.wav [file2.wav file3.wav ...]\n\n" +
				"For each file the filename is printed followed by ordered pairs corresponding to predicted word endpoints."
				   
	def main(args: Array[String]) {
		if (args.length < 1)	usageAndExit()
		val files = 
			for {arg <- args} 
				yield new java.io.File(arg)
		for {file <- files} {
			if (file.exists == false) {
				Console.err.println(file + " does not exist!")
				System.exit(1)
			}
		}
	}
	
	def usageAndExit() {
		Console.err.println(usage)
		System.exit(1)
	}
}