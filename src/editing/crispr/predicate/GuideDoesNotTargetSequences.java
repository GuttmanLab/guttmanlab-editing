package editing.crispr.predicate;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;

import nextgen.core.utils.AlignmentUtils;

import org.apache.log4j.Logger;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;

import editing.crispr.GuideRNA;

/**
 * Check whether a guide RNA potentially targets some sequence in a set of interest
 * @author prussell
 *
 */
public class GuideDoesNotTargetSequences implements GuideRNAPredicate {
	
	private Collection<Sequence> sequences;
	private int maxMatchesFirst16nt;
	private int maxMatchesLast4nt;
	private Collection<GuideRNA> guidesInSequences;
	public static Logger logger = Logger.getLogger(GuideDoesNotTargetSequences.class.getName());
	
	/**
	 * @param sequenceFasta Sequences that should not be targeted
	 * @param maxMatchesFirst16 Max tolerable number of matches of the first 16nt to a region of another sequence that is followed by a PAM
	 * @param maxMatchesLast4 Max tolerable number of matches in the last 4nt
	 * @param includeNAG Include NAGs as possible PAM sequences in the other sequences
	 */
	public GuideDoesNotTargetSequences(String sequenceFasta, int maxMatchesFirst16, int maxMatchesLast4, boolean includeNAG) throws IOException {
		sequences = FastaSequenceIO.loadSequences(new File(sequenceFasta));
		maxMatchesFirst16nt = maxMatchesFirst16;
		maxMatchesLast4nt = maxMatchesLast4;
		guidesInSequences = new ArrayList<GuideRNA>();
		for(Sequence seq : sequences) {
			guidesInSequences.addAll(GuideRNA.findAll(seq, 0, seq.getLength(), null, includeNAG));
		}
	}
	
	/**
	 * @param g1 Guide RNA 1
	 * @param g2 Guide RNA 2
	 * @return True if the guides have a match over the first 16nt or the last 4nt
	 */
	private boolean match(GuideRNA g1, GuideRNA g2) {
		String seq1 = g1.getSequenceString();
		String seq2 = g2.getSequenceString();
		String seq1prefix = seq1.substring(0, 16);
		String seq2prefix = seq2.substring(0, 16);
		boolean matchFirst16 = AlignmentUtils.hammingDistanceAtMost(seq1prefix, seq2prefix, 16 - maxMatchesFirst16nt - 1, true);
		String seq1suffix = seq1.substring(16);
		String seq2suffix = seq2.substring(16);
		boolean matchLast4 = AlignmentUtils.hammingDistanceAtMost(seq1suffix, seq2suffix, 4 - maxMatchesLast4nt - 1, true);
		if(matchFirst16 && matchLast4) {
			//logger.debug("Guides match over first 16nt and last 4nt: " + g1.getSequenceString() + " " + g2.getSequenceString());
			return true;
		}
		return false;
	}
	
	/**
	 * @param g Guide RNA
	 * @param others Other sequences to check for possible targeting
	 * @return True if some other sequence contains a PAM sequence and guide RNA that matches the provided guide over first 16nt or last 4nt
	 */
	private boolean match(GuideRNA g, Collection<GuideRNA> others) {
		for(GuideRNA guide : others) {
			if(match(g, guide)) {
				return true;
			}
		}
		return false;
	}
	
	@Override
	public boolean evaluate(GuideRNA g) {
		return !match(g, guidesInSequences);
	}

	@Override
	public String getPredicateName() {
		return "guide_does_not_target_sequences";
	}

	@Override
	public String getShortFailureMessage(GuideRNA g) throws IOException, InterruptedException {
		return "guide_targets_some_sequence_in_set";
	}

}
