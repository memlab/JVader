package edu.upenn.psych.memory.jvader

import java.io.File

trait VoiceEndpointer {

	
	/**
	 * Identifies segments of human speech in an audio file.
	 *
	 * @param audioFile - A mono audio file in wav, aiff, or au formats 
	 *
	 * @return A <tt>List</tt> of tuples of the form (voiceStarts, voiceEnds), in order, in frames
	 */
	def findEndpoints(audioFile: File): List[(Int, Int)]
}