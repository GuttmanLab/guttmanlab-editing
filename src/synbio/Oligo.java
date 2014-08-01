package synbio;

import primer.PrimerPair;
import broad.core.sequence.Sequence;

public class Oligo {
	/**
	 * Check whether a primer pair is compatible with a probe
	 * @param primerPair Primer pair
	 * @param probeSequenceExcludingPrimers Probe sequence without primers
	 * @return True if there are no matches of the primers to the probe
	 */
	public static boolean primerPairCompatibleWithProbe(PrimerPair primerPair, String probeSequenceExcludingPrimers) {
		return primerPairCompatibleWithProbe(primerPair.getLeftPrimer(), primerPair.getRightPrimer(), probeSequenceExcludingPrimers);
	}
	
	/**
	 * Check whether a primer pair is compatible with a probe
	 * @param leftPrimer Left primer
	 * @param rightPrimer Right primer
	 * @param probeSequenceExcludingPrimers Probe sequence without primers
	 * @return True if there are no matches of the primers to the probe
	 */
	public static boolean primerPairCompatibleWithProbe(String leftPrimer, String rightPrimer, String probeSequenceExcludingPrimers) {
		String leftPrimer3primeEnd = leftPrimer.substring(leftPrimer.length() - 8);
		String leftPrimer3primeEndRC = Sequence.reverseSequence(leftPrimer3primeEnd);
		String rightPrimer3primeEnd = rightPrimer.substring(rightPrimer.length() - 8);
		String rightPrimer3primeEndRC = Sequence.reverseSequence(rightPrimer3primeEnd);
		// Check that the two primer ends do not appear in the oligo
		if(probeSequenceExcludingPrimers.contains(leftPrimer3primeEnd)) {
			return false;
		}
		if(probeSequenceExcludingPrimers.contains(rightPrimer3primeEndRC)) {
			return false;
		}
		// Check that the RCed primer ends do not appear at all in the oligo
		if(probeSequenceExcludingPrimers.contains(leftPrimer3primeEndRC)) {
			return false;
		}
		if(probeSequenceExcludingPrimers.contains(rightPrimer3primeEnd)) {
			return false;
		}
		return true;
	}
	
	/**
	 * Check whether a primer pair is compatible with a full oligo sequence
	 * @param primerPair Primer pair
	 * @param oligoSequenceIncludingPrimers Full oligo sequence including this primer pair at the ends
	 * @return True if there are no nonspecific matches of the primers to the oligo
	 */
	public static boolean primerPairCompatibleWithFullOligo(PrimerPair primerPair, String oligoSequenceIncludingPrimers) {
		return primerPairCompatibleWithFullOligo(primerPair.getLeftPrimer(), primerPair.getRightPrimer(), oligoSequenceIncludingPrimers);
	}

	/**
	 * Check whether a primer pair is compatible with a full oligo sequence
	 * @param leftPrimer Left primer
	 * @param rightPrimer Right primer
	 * @param oligoSequenceIncludingPrimers Full oligo sequence including this primer pair at the ends
	 * @return True if there are no nonspecific matches of the primers to the oligo
	 */
	public static boolean primerPairCompatibleWithFullOligo(String leftPrimer, String rightPrimer, String oligoSequenceIncludingPrimers) {
		String leftPrimer3primeEnd = leftPrimer.substring(leftPrimer.length() - 8);
		String leftPrimer3primeEndRC = Sequence.reverseSequence(leftPrimer3primeEnd);
		String rightPrimer3primeEnd = rightPrimer.substring(rightPrimer.length() - 8);
		String rightPrimer3primeEndRC = Sequence.reverseSequence(rightPrimer3primeEnd);
		// Check the each oligo contains the 3' ends of the primers once
		int firstOccurrenceLeftPrimer3primeEnd = oligoSequenceIncludingPrimers.indexOf(leftPrimer3primeEnd);
		int lastOccurrenceLeftPrimer3primeEnd = oligoSequenceIncludingPrimers.lastIndexOf(leftPrimer3primeEnd);
		int firstOccurrenceRightPrimer3primeEndRC = oligoSequenceIncludingPrimers.indexOf(rightPrimer3primeEndRC);
		int lastOccurrenceRightPrimer3primeEndRC = oligoSequenceIncludingPrimers.lastIndexOf(rightPrimer3primeEndRC);
		// Check that the two primer ends appear in the oligo
		if(firstOccurrenceLeftPrimer3primeEnd == -1) {
			return false;
		}
		if(firstOccurrenceRightPrimer3primeEndRC == -1) {
			return false;
		}
		// Check that the two primer ends appear at most once in the oligo
		if(firstOccurrenceLeftPrimer3primeEnd != lastOccurrenceLeftPrimer3primeEnd) {
			return false;
		}
		if(firstOccurrenceRightPrimer3primeEndRC != lastOccurrenceRightPrimer3primeEndRC) {
			return false;
		}
		// Check that the RCed primer ends do not appear at all in the oligo
		if(oligoSequenceIncludingPrimers.contains(leftPrimer3primeEndRC)) {
			return false;
		}
		if(oligoSequenceIncludingPrimers.contains(rightPrimer3primeEnd)) {
			return false;
		}
		return true;
	}
	

}
