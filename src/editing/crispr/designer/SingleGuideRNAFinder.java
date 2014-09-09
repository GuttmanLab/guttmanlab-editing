package editing.crispr.designer;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;


import nextgen.core.annotation.Gene;

import org.apache.log4j.Logger;

import editing.crispr.GuideRNA;
import editing.crispr.predicate.GuideDoesNotTargetSequences;
import editing.crispr.predicate.GuideSufficientEfficacy;
import editing.crispr.score.GuideEfficacyScore;
import editing.crispr.score.GuideOffTargetScore;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;

/**
 * Find singleton guide RNAs contained within intervals of chromosomes
 * Ability to add filters
 * @author prussell
 *
 */
public class SingleGuideRNAFinder {

	public static Logger logger = Logger.getLogger(SingleGuideRNAFinder.class.getName());
	private Map<String, Sequence> chrsByName;
	
	// Off target filter
	private static boolean ENFORCE_MIN_OFF_TARGET_SCORE = false;
	public static int MIN_OFF_TARGET_SCORE = 30;
	private static File OFF_TARGET_BITS = null;

	// Guide efficacy filter
	private static boolean ENFORCE_MAX_GUIDE_EFFICACY_SCORE = false;
	public static double MAX_GUIDE_EFFICACY_SCORE = 0.6;

	// Filter for potential to target another set of sequences
	private static boolean ENFORCE_DOES_NOT_TARGET_OTHER_SEQS = false;
	private static GuideDoesNotTargetSequences GUIDE_DOES_NOT_TARGET_SEQS;
	
	/**
	 * @param genomeFasta Genome fasta file
	 * @throws IOException
	 */
	public SingleGuideRNAFinder(String genomeFasta) throws IOException {
		logger.info("");
		logger.info("Instantiating guide RNA finder for sequences in " + genomeFasta + "...");
		chrsByName = FastaSequenceIO.getChrSequencesFromFasta(genomeFasta);
	}
	
	/**
	 * Add filter to ensure that guide RNAs do not target any sequence in another set of sequences
	 * @param sequenceFasta Sequences that should not be targeted
	 * @param maxMatches20 Max tolerable number of matches of the 20nt sequence to a region of another sequence that is followed by a PAM
	 * @param maxMatchesLast10 Max tolerable number of matches in the last 10nt
	 * @param includeNAG Include NAGs as possible PAM sequences in the other sequences
	 * @throws IOException
	 */
	public void addGuideDoesNotTargetOtherSeqs(String sequenceFasta, int maxMatches20, int maxMatchesLast10, boolean includeNAG) throws IOException {
		ENFORCE_DOES_NOT_TARGET_OTHER_SEQS = true;
		GUIDE_DOES_NOT_TARGET_SEQS = new GuideDoesNotTargetSequences(sequenceFasta, maxMatches20, maxMatchesLast10, includeNAG);
	}
	
	/**
	 * Add off target filter with default min off target score
	 * @param offTargetBitsFile Bed file containing off target sites, with the 23-base off-target sequence in the name column
	 */
	public void addOffTargetFilter(String offTargetBitsFile) {
		addOffTargetFilter(offTargetBitsFile, MIN_OFF_TARGET_SCORE);
	}
	
	/**
	 * Add off target filter
	 * @param offTargetBitsFile Bed file containing off target sites, with the 23-base off-target sequence in the name column
	 * @param minScore Min off target score
	 */
	public void addOffTargetFilter(String offTargetBitsFile, int minScore) {
		logger.info("");
		logger.info("Adding off target filter with off target file " + offTargetBitsFile + " and min score " + minScore + ".");
		ENFORCE_MIN_OFF_TARGET_SCORE = true;
		MIN_OFF_TARGET_SCORE = minScore;
		OFF_TARGET_BITS = new File(offTargetBitsFile);
	}
	
	/**
	 * Add guide efficacy filter with default max score
	 */
	public void addEfficacyFilter() {
		addEfficacyFilter(MAX_GUIDE_EFFICACY_SCORE);
	}
	
	/**
	 * Add guide efficacy filter
	 * @param maxScore Max efficacy score
	 */
	public void addEfficacyFilter(double maxScore) {
		logger.info("");
		logger.info("Adding efficacy filter with max score " + maxScore + ".");
		ENFORCE_MAX_GUIDE_EFFICACY_SCORE = true;
		MAX_GUIDE_EFFICACY_SCORE = maxScore;
	}
	
	/**
	 * Filter the collection of guides in place
	 * @param guides Guide RNA collection to filter
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private void applyFilters(Collection<GuideRNA> guides) throws IOException, InterruptedException {
		// Apply guide efficacy filter
		if(ENFORCE_MAX_GUIDE_EFFICACY_SCORE) {
			GuideSufficientEfficacy ge = new GuideSufficientEfficacy(new GuideEfficacyScore(guides), MAX_GUIDE_EFFICACY_SCORE);
			Iterator<GuideRNA> iter = guides.iterator();
			while(iter.hasNext()) {
				GuideRNA guide = iter.next();
				if(!ge.evaluate(guide)) {
					//logger.debug("Guide " + guide.toString() + " removed by efficacy filter.");
					iter.remove();
				} else {
					//logger.debug("Guide " + guide.toString() + " passes efficacy filter.");
				}
			}
		}
		// Apply guide off target filter
		if (ENFORCE_MIN_OFF_TARGET_SCORE) {
			GuideOffTargetScore scorer = new GuideOffTargetScore(OFF_TARGET_BITS);
			Iterator<GuideRNA> iter = guides.iterator();
			while(iter.hasNext()) {
				GuideRNA guide = iter.next();
				if(scorer.getScore(guide) < MIN_OFF_TARGET_SCORE) {
					logger.debug("Guide " + guide.toString() + " removed by off target filter.");
					iter.remove();
				} else {
					logger.debug("Guide " + guide.toString() + " passes off target filter.");
				}
			}
		}
		// Apply match filter
		if(ENFORCE_DOES_NOT_TARGET_OTHER_SEQS) {
			Iterator<GuideRNA> iter = guides.iterator();
			while(iter.hasNext()) {
				GuideRNA guide = iter.next();
				if(!GUIDE_DOES_NOT_TARGET_SEQS.evaluate(guide)) {
					//logger.debug("Guide " + guide.toString() + " removed by other seqs filter.");
					iter.remove();
				} else {
					//logger.debug("Guide " + guide.toString() + " passes other seqs filter.");
				}
			}
		}
	}
	
	/**
	 * Get filtered set of guide RNAs within the region
	 * Get a certain number of valid guides
	 * No guarantee about which ones will be returned
	 * @param regionChr Region chromosome
	 * @param regionStart Region start
	 * @param regionEnd Region end
	 * @param numToGet Number of guides to get
	 * @return All valid guide RNAs contained within the region
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public Collection<GuideRNA> getFilteredGuideRNAs(String regionChr, int regionStart, int regionEnd, int numToGet) throws IOException, InterruptedException {
		Gene region = new Gene(regionChr, regionStart, regionEnd);
		region.setName(region.toUCSC());
		Collection<GuideRNA> guides = GuideRNA.findAll(chrsByName.get(region.getChr()), region.getStart(), region.getEnd(), region);
		// Use a hash set so order will be somewhat random
		HashSet<GuideRNA> guidesHash = new HashSet<GuideRNA>();
		guidesHash.addAll(guides);
		Collection<GuideRNA> rtrn = new ArrayList<GuideRNA>();
		for(GuideRNA guide : guidesHash) {
			Collection<GuideRNA> thisGuide = new ArrayList<GuideRNA>();
			thisGuide.add(guide);
			try {
				applyFilters(thisGuide);
			} catch(IllegalStateException e) {
				logger.warn("CAUGHT EXCEPTION");
				e.printStackTrace();
			}
			if(!thisGuide.isEmpty()) {
				rtrn.addAll(thisGuide);
			}
			if(rtrn.size() >= numToGet) {
				return rtrn;
			}
		}
		logger.warn("Only found " + rtrn.size() + " valid guide RNAs for region " + region.toUCSC());
		return rtrn;
	}

	
	/**
	 * Get filtered set of guide RNAs within the region
	 * @param regionChr Region chromosome
	 * @param regionStart Region start
	 * @param regionEnd Region end
	 * @return All valid guide RNAs contained within the region
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public Collection<GuideRNA> getFilteredGuideRNAs(String regionChr, int regionStart, int regionEnd) throws IOException, InterruptedException {
		//logger.debug("");
		//logger.debug("Getting filtered guide RNAs for region " + regionChr + ":" + regionStart + "-" + regionEnd);
		Gene region = new Gene(regionChr, regionStart, regionEnd);
		region.setName(region.toUCSC());
		Collection<GuideRNA> guides = GuideRNA.findAll(chrsByName.get(region.getChr()), region.getStart(), region.getEnd(), region);
		//logger.debug("There are " + guides.size() + " possible guide RNAs.");
		try {
			applyFilters(guides);
		} catch(IllegalStateException e) {
			logger.warn("CAUGHT EXCEPTION, RETURNING EMPTY SET FOR REGION " + region.toUCSC());
			e.printStackTrace();
		}
		//logger.debug("There are " + guides.size() + " guide RNAs that pass filters.");
		if(guides.isEmpty()) {
			logger.debug("NO VALID GUIDE RNAS FOR REGION " + region.toUCSC());
		}
		return guides;
	}

	/**
	 * @return The reference sequences by name
	 */
	public Map<String, Sequence> getReferenceSequences() {
		return chrsByName;
	}
	

}
