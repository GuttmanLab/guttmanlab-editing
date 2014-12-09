package editing.crispr;

import java.util.ArrayList;
import java.util.Collection;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;
import nextgen.core.general.TabbedReader;
import broad.core.error.ParseException;
import broad.core.sequence.Sequence;

/**
 * A guide RNA with a genomic position and a target gene for CRISPR editing
 * @author prussell
 */
public class GuideRNA extends BasicAnnotation {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	public static Logger logger = Logger.getLogger(GuideRNA.class.getName());
	private Sequence sequenceWithoutPAM, sequenceWithPAM;
	private Gene target = null;
	public static boolean ALLOW_NAG_AS_PAM = false;

	
	/**
	 * @param targetGene Target gene
	 * @param chromosome Chromosome
	 * @param start Start position of guide sequence without PAM
	 * @param end Position after last position of guide sequence without PAM
	 * @param orientation Strand
	 */
	public GuideRNA(Gene targetGene, Sequence chromosome, int start, int end, Strand orientation) {
		super(chromosome.getId(), start, end, orientation);
		validateStartEnd(start, end);
		validateStrand(orientation);
		setNameFromTarget();
		target = targetGene;
		
		int startWithoutPAM;
		int startWithPAM;
		int endWithoutPAM;
		int endWithPAM;
		Strand strand;
		
		// Validate
		strand = orientation;
		if(strand.equals(Strand.POSITIVE)) {
			startWithoutPAM = getStart();
			startWithPAM = startWithoutPAM;
			endWithoutPAM = getEnd();
			endWithPAM = endWithoutPAM + 3;
		} else {
			startWithoutPAM = getStart();
			startWithPAM = getStart() - 3;
			endWithoutPAM = getEnd();
			endWithPAM = getEnd();
		}
		Sequence shortSeq = chromosome.getSubSequence("", startWithoutPAM, endWithoutPAM);
		shortSeq.setSequenceBases(shortSeq.getSequenceBases().toUpperCase());
		Sequence longSeq = chromosome.getSubSequence("", startWithPAM, endWithPAM);
		longSeq.setSequenceBases(longSeq.getSequenceBases().toUpperCase());
		sequenceWithoutPAM = strand.equals(Strand.POSITIVE) ? shortSeq : Sequence.reverseSequence(shortSeq);
		sequenceWithPAM = strand.equals(Strand.POSITIVE) ? longSeq : Sequence.reverseSequence(longSeq);
		validateSequenceWithPAM(sequenceWithPAM);
	}

	/**
	 * @param a The annotation of the guide RNA itself without PAM sequence
	 * @param sequenceIncludingPAM Sequence of the guide RNA with PAM
	 */
	public GuideRNA(Annotation a, String sequenceIncludingPAM) {
		super(a);
		validateStartEnd(getStart(), getEnd());
		validateStrand(getOrientation());
		setNameFromTarget();
		sequenceWithoutPAM = new Sequence(getName());
		sequenceWithoutPAM.setSequenceBases(sequenceIncludingPAM.substring(0,sequenceIncludingPAM.length() - 3));
		sequenceWithoutPAM.setForwardStrand(getOrientation() == Strand.POSITIVE);
		
		sequenceWithPAM = new Sequence(getName());
		sequenceWithPAM.setSequenceBases(sequenceIncludingPAM);
		sequenceWithPAM.setForwardStrand(getOrientation() == Strand.POSITIVE);
		
		validateSequenceWithPAM(sequenceWithPAM);
		
	}
	
	
	private void setNameFromTarget() {
		setName(getReferenceName() + ":" + getStart() + "-" + getEnd() + ":" + getOrientation().toString());
	}
	
	
	public boolean isPlusStrand() {
		return getStrand().equals(Strand.POSITIVE);
	}
	
	public boolean isMinusStrand() {
		return getStrand().equals(Strand.NEGATIVE);
	}
	
	public String getSequenceString() {
		return sequenceWithoutPAM.getSequenceBases();
	}
	
	public Sequence getSequence() {
		return sequenceWithoutPAM;
	}
	
	public Sequence getSequenceWithPAM() {
		return sequenceWithPAM;
	}
	
	public String getSequenceStringWithPAM() {
		return getSequenceWithPAM().getSequenceBases();
	}
	
	public Gene getTargetGene() {
		return target;
	}

	/**
	 * Find all 20nt guide RNAs on either strand within the window
	 * @param chr Chromosome
	 * @param start Start position of window to look in
	 * @param end Position after last position of window
	 * @return All guide RNAs followed by PAM sequence whose sequence without PAM is fully contained in the window
	 */
	public static Collection<GuideRNA> findAll(Sequence chr, int start, int end, Gene targetGene, boolean allowNAG) {
		return findAll(chr, start, end, 20, 20, targetGene, allowNAG);
	}
	
	/**
	 * Find all guide RNAs on either strand within the window
	 * @param chr Chromosome
	 * @param start Start position of window to look in
	 * @param end Position after last position of window
	 * @param minLen Min length of guide sequence without PAM
	 * @param maxLen Max length of guide sequence without PAM
	 * @return All guide RNAs followed by PAM sequence whose sequence without PAM is fully contained in the window
	 */
	public static Collection<GuideRNA> findAll(Sequence chr, int start, int end, int minLen, int maxLen, Gene targetGene, boolean allowNAG) {
		//logger.debug("");
		ALLOW_NAG_AS_PAM = allowNAG;
		Collection<GuideRNA> rtrn = new ArrayList<GuideRNA>();
		Collection<Annotation> pams = findAllNGGs(chr, start, end);
		if(ALLOW_NAG_AS_PAM) {
			pams.addAll(findAllNAGs(chr, start, end));
		}
		for(Annotation pam : pams) {
			for(int len = minLen; len <= maxLen; len++) {
				GuideRNA g = adjacentGuideRNA(chr, start, end, len, pam, targetGene);
				if(g != null) {
					//logger.debug("Added " + g.toString());
					rtrn.add(g);
				}
			}
		}
		return rtrn;
	}
	
	public String toString() {
		return getReferenceName() + ":" + getStart() + "-" + getEnd() + ":" + getStrand().toString() + ":" + sequenceWithoutPAM.getSequenceBases();
	}
	
	/**
	 * Get the correctly oriented guide RNA that is fully contained in the window and ends with the PAM specified
	 * @param windowChr Window chromosome
	 * @param windowStart Window start
	 * @param windowEnd Position after last position of window
	 * @param length Length of guide RNA sequence without PAM
	 * @param pam Oriented PAM location
	 * @param targetGene Target gene
	 * @return The guide RNA or null if not fully contained in window
	 */
	private static GuideRNA adjacentGuideRNA(Sequence windowChr, int windowStart, int windowEnd, int length, Annotation pam, Gene targetGene) {
		validatePAM(windowChr, pam);
		//logger.debug("Getting guide RNA adjacent to " + pam.toUCSC() + ":" + pam.getOrientation().toString());
		if(pam.getOrientation().equals(Strand.POSITIVE)) {
			if(pam.getStart() - windowStart < length) {
				// Guide RNA cannot be fully contained in window
				logger.debug("Guide RNA neighboring " + pam.toUCSC() +" not fully contained in " + windowChr.getId() + ":" + windowStart + "-" + windowEnd);
				return null;
			}
			GuideRNA rtrn = new GuideRNA(targetGene, windowChr, pam.getStart() - length, pam.getStart(), Strand.POSITIVE);
			//logger.debug("Guide RNA neighboring " + pam.toUCSC() + " is " + rtrn.toString());
			return rtrn;
		}
		if(pam.getOrientation().equals(Strand.NEGATIVE)) {
			if(pam.getEnd() + length > windowEnd) {
				// Guide RNA cannot be fully contained in window
				logger.debug("Guide RNA neighboring " + pam.toUCSC() +" not fully contained in " + windowChr.getId() + ":" + windowStart + "-" + windowEnd);
				return null;
			}
			GuideRNA rtrn = new GuideRNA(targetGene, windowChr, pam.getEnd(), pam.getEnd() + length, Strand.NEGATIVE);
			//logger.debug("Guide RNA neighboring " + pam.toUCSC() + " is " + rtrn.toString());
			return rtrn;
		}
		throw new IllegalArgumentException("Strand must be known");
	}
	
	private static void validatePAM(Sequence chr, Annotation pam) {
		if(pam.numBlocks() != 1) {
			throw new IllegalArgumentException("Must have one block");
		}
		Sequence seq = chr.getSubsequence(pam);
		String substr = seq.getSequenceBases().substring(1).toUpperCase();
		if(!ALLOW_NAG_AS_PAM && !substr.equals("GG")) {
			throw new IllegalArgumentException("Sequence must be NGG. Is " + seq.getSequenceBases());
		}
		if(ALLOW_NAG_AS_PAM && !(substr.equals("GG") || substr.equals("AG"))) {
			throw new IllegalArgumentException("Sequence must be NGG or NAG. Is " + seq.getSequenceBases());
		}
	}
	
	/**
	 * Get all NAG sequences that are fully contained within the window on either strand
	 * @param chr Chromosome
	 * @param start First position of window
	 * @param end Position after last position of window
	 * @return All NAGs as annotations with strand
	 */
	private static Collection<Annotation> findAllNAGs(Sequence chr, int start, int end) {
		String seq = chr.getSubSequence("", start, end).getSequenceBases();
		Collection<Annotation> rtrn = new TreeSet<Annotation>();
		
		// Find all NAGs
		int i1 = 0;
		while(i1 < seq.length()) {
			int pos = seq.indexOf("AG", i1);
			if(pos == -1) {
				break;
			}
			if(pos == 0) {
				// Won't be fully contained
				i1++;
				continue;
			}
			if(pos + 2 > seq.length()) {
				// Won't be fully contained
				break;
			}
			// Create annotation
			Annotation pam = new BasicAnnotation(chr.getId(), start + pos - 1, start + pos + 2, Strand.POSITIVE);
			rtrn.add(pam);
			i1 = pos + 1;
		}
		
		int i2 = 0;
		while(i2 < seq.length()) {
			int pos = seq.indexOf("ag", i2);
			if(pos == -1) {
				break;
			}
			if(pos == 0) {
				// Won't be fully contained
				i2++;
				continue;
			}
			if(pos + 2 > seq.length()) {
				// Won't be fully contained
				break;
			}
			// Create annotation
			Annotation pam = new BasicAnnotation(chr.getId(), start + pos - 1, start + pos + 2, Strand.POSITIVE);
			rtrn.add(pam);
			i2 = pos + 1;
		}
		
		int i3 = 0;
		while(i3 < seq.length()) {
			int pos = seq.indexOf("Ag", i3);
			if(pos == -1) {
				break;
			}
			if(pos == 0) {
				// Won't be fully contained
				i3++;
				continue;
			}
			if(pos + 2 > seq.length()) {
				// Won't be fully contained
				break;
			}
			// Create annotation
			Annotation pam = new BasicAnnotation(chr.getId(), start + pos - 1, start + pos + 2, Strand.POSITIVE);
			rtrn.add(pam);
			i3 = pos + 1;
		}
		
		int i4 = 0;
		while(i4 < seq.length()) {
			int pos = seq.indexOf("aG", i4);
			if(pos == -1) {
				break;
			}
			if(pos == 0) {
				// Won't be fully contained
				i4++;
				continue;
			}
			if(pos + 2 > seq.length()) {
				// Won't be fully contained
				break;
			}
			// Create annotation
			Annotation pam = new BasicAnnotation(chr.getId(), start + pos - 1, start + pos + 2, Strand.POSITIVE);
			rtrn.add(pam);
			i4 = pos + 1;
		}
		
		// Find all CCNs
		int j1 = 0;
		while(j1 < seq.length()) {
			int pos = seq.indexOf("CT", j1);
			if(pos == -1) break;
			if(pos + 3 > seq.length()) {
				// Won't be fully contained
				break;
			}
			// Create annotation
			Annotation ccn = new BasicAnnotation(chr.getId(), start + pos, start + pos + 3, Strand.NEGATIVE);
			rtrn.add(ccn);
			j1 = pos + 1;
		}
		
		int j2 = 0;
		while(j2 < seq.length()) {
			int pos = seq.indexOf("cT", j2);
			if(pos == -1) break;
			if(pos + 3 > seq.length()) {
				// Won't be fully contained
				break;
			}
			// Create annotation
			Annotation ccn = new BasicAnnotation(chr.getId(), start + pos, start + pos + 3, Strand.NEGATIVE);
			rtrn.add(ccn);
			j2 = pos + 1;
		}
		
		int j3 = 0;
		while(j3 < seq.length()) {
			int pos = seq.indexOf("ct", j3);
			if(pos == -1) break;
			if(pos + 3 > seq.length()) {
				// Won't be fully contained
				break;
			}
			// Create annotation
			Annotation ccn = new BasicAnnotation(chr.getId(), start + pos, start + pos + 3, Strand.NEGATIVE);
			rtrn.add(ccn);
			j3 = pos + 1;
		}
		
		int j4 = 0;
		while(j4 < seq.length()) {
			int pos = seq.indexOf("Ct", j4);
			if(pos == -1) break;
			if(pos + 3 > seq.length()) {
				// Won't be fully contained
				break;
			}
			// Create annotation
			Annotation ccn = new BasicAnnotation(chr.getId(), start + pos, start + pos + 3, Strand.NEGATIVE);
			rtrn.add(ccn);
			j4 = pos + 1;
		}
		
		return rtrn;
		
	}

	
	/**
	 * Get all NGG sequences that are fully contained within the window on either strand
	 * @param chr Chromosome
	 * @param start First position of window
	 * @param end Position after last position of window
	 * @return All NGGs as annotations with strand
	 */
	private static Collection<Annotation> findAllNGGs(Sequence chr, int start, int end) {
		String seq = chr.getSubSequence("", start, end).getSequenceBases();
		Collection<Annotation> rtrn = new TreeSet<Annotation>();
		
		// Find all NGGs
		int i1 = 0;
		while(i1 < seq.length()) {
			int pos = seq.indexOf("GG", i1);
			if(pos == -1) {
				break;
			}
			if(pos == 0) {
				// Won't be fully contained
				i1++;
				continue;
			}
			if(pos + 2 > seq.length()) {
				// Won't be fully contained
				break;
			}
			// Create annotation
			Annotation ngg = new BasicAnnotation(chr.getId(), start + pos - 1, start + pos + 2, Strand.POSITIVE);
			//logger.debug("NGG in " + chr.getId() + ":" + start + "-" + end + "\t" + ngg.toUCSC() + ":" + ngg.getOrientation().toString());
			rtrn.add(ngg);
			i1 = pos + 1;
		}
		
		int i2 = 0;
		while(i2 < seq.length()) {
			int pos = seq.indexOf("gg", i2);
			if(pos == -1) {
				break;
			}
			if(pos == 0) {
				// Won't be fully contained
				i2++;
				continue;
			}
			if(pos + 2 > seq.length()) {
				// Won't be fully contained
				break;
			}
			// Create annotation
			Annotation ngg = new BasicAnnotation(chr.getId(), start + pos - 1, start + pos + 2, Strand.POSITIVE);
			//logger.debug("NGG in " + chr.getId() + ":" + start + "-" + end + "\t" + ngg.toUCSC() + ":" + ngg.getOrientation().toString());
			rtrn.add(ngg);
			i2 = pos + 1;
		}
		
		int i3 = 0;
		while(i3 < seq.length()) {
			int pos = seq.indexOf("Gg", i3);
			if(pos == -1) {
				break;
			}
			if(pos == 0) {
				// Won't be fully contained
				i3++;
				continue;
			}
			if(pos + 2 > seq.length()) {
				// Won't be fully contained
				break;
			}
			// Create annotation
			Annotation ngg = new BasicAnnotation(chr.getId(), start + pos - 1, start + pos + 2, Strand.POSITIVE);
			//logger.debug("NGG in " + chr.getId() + ":" + start + "-" + end + "\t" + ngg.toUCSC() + ":" + ngg.getOrientation().toString());
			rtrn.add(ngg);
			i3 = pos + 1;
		}
		
		int i4 = 0;
		while(i4 < seq.length()) {
			int pos = seq.indexOf("gG", i4);
			if(pos == -1) {
				break;
			}
			if(pos == 0) {
				// Won't be fully contained
				i4++;
				continue;
			}
			if(pos + 2 > seq.length()) {
				// Won't be fully contained
				break;
			}
			// Create annotation
			Annotation ngg = new BasicAnnotation(chr.getId(), start + pos - 1, start + pos + 2, Strand.POSITIVE);
			//logger.debug("NGG in " + chr.getId() + ":" + start + "-" + end + "\t" + ngg.toUCSC() + ":" + ngg.getOrientation().toString());
			rtrn.add(ngg);
			i4 = pos + 1;
		}
		
		// Find all CCNs
		int j1 = 0;
		while(j1 < seq.length()) {
			int pos = seq.indexOf("CC", j1);
			if(pos == -1) break;
			if(pos + 3 > seq.length()) {
				// Won't be fully contained
				break;
			}
			// Create annotation
			Annotation ccn = new BasicAnnotation(chr.getId(), start + pos, start + pos + 3, Strand.NEGATIVE);
			//logger.debug("NGG in " + chr.getId() + ":" + start + "-" + end + "\t" + ccn.toUCSC() + ":" + ccn.getOrientation().toString());
			rtrn.add(ccn);
			j1 = pos + 1;
		}
		
		int j2 = 0;
		while(j2 < seq.length()) {
			int pos = seq.indexOf("cC", j2);
			if(pos == -1) break;
			if(pos + 3 > seq.length()) {
				// Won't be fully contained
				break;
			}
			// Create annotation
			Annotation ccn = new BasicAnnotation(chr.getId(), start + pos, start + pos + 3, Strand.NEGATIVE);
			//logger.debug("NGG in " + chr.getId() + ":" + start + "-" + end + "\t" + ccn.toUCSC() + ":" + ccn.getOrientation().toString());
			rtrn.add(ccn);
			j2 = pos + 1;
		}
		
		int j3 = 0;
		while(j3 < seq.length()) {
			int pos = seq.indexOf("cc", j3);
			if(pos == -1) break;
			if(pos + 3 > seq.length()) {
				// Won't be fully contained
				break;
			}
			// Create annotation
			Annotation ccn = new BasicAnnotation(chr.getId(), start + pos, start + pos + 3, Strand.NEGATIVE);
			//logger.debug("NGG in " + chr.getId() + ":" + start + "-" + end + "\t" + ccn.toUCSC() + ":" + ccn.getOrientation().toString());
			rtrn.add(ccn);
			j3 = pos + 1;
		}
		
		int j4 = 0;
		while(j4 < seq.length()) {
			int pos = seq.indexOf("Cc", j4);
			if(pos == -1) break;
			if(pos + 3 > seq.length()) {
				// Won't be fully contained
				break;
			}
			// Create annotation
			Annotation ccn = new BasicAnnotation(chr.getId(), start + pos, start + pos + 3, Strand.NEGATIVE);
			//logger.debug("NGG in " + chr.getId() + ":" + start + "-" + end + "\t" + ccn.toUCSC() + ":" + ccn.getOrientation().toString());
			rtrn.add(ccn);
			j4 = pos + 1;
		}
		
		return rtrn;
		
	}
	
	private static void validateStartEnd(int start, int end) {
		if(!(end > start)) {
			throw new IllegalArgumentException("End must be > start");
		}
	}
	
	private static void validateStrand(Strand orientation) {
		if(orientation.equals(Strand.UNKNOWN)) {
			throw new IllegalArgumentException("Strand cannot be unknown");
		}
	}
	
	private static void validateSequenceWithPAM(Sequence sequence) {
		String substr = sequence.getSequenceBases().substring(sequence.getLength() - 2, sequence.getLength());
		if(!ALLOW_NAG_AS_PAM && !substr.equals("GG")) {
			throw new IllegalArgumentException("Sequence must end in GG: " + sequence.getSequenceBases());
		}
		if(ALLOW_NAG_AS_PAM && !(substr.equals("GG") || substr.equals("AG"))) {
			throw new IllegalArgumentException("Sequence must end in GG or AG: " + sequence.getSequenceBases());
		}
	}

	public int hashCode() {
		String s = toBED() + getName() + sequenceWithoutPAM.getSequenceBases();
		if (target != null) s = s + target.toBED();
		return s.hashCode();
	}
	
	public boolean equals(Object o) {
		if(!o.getClass().equals(getClass())) return false;
		GuideRNA g = (GuideRNA)o;
		return hashCode() == g.hashCode();
	}

	
	public String toBedWithSequence() {
		return toBED() + "\t" + sequenceWithoutPAM.getSequenceBases();
	}

	/**
	 * @author engreitz
	 * For files stored in BED12 plus the sequence including PAM in column 13
	 */
	public static class Factory implements TabbedReader.Factory<GuideRNA> {
		@Override
		public GuideRNA create(String[] rawFields) throws ParseException {
			Annotation a = new BasicAnnotation.Factory().create(rawFields);
			return new GuideRNA(a, rawFields[12]);
		}
	}
	
	/**
	 * @author engreitz
	 * For files stored in BED6 with sequence including PAM in the name column
	 */
	public static class FactoryBED6 implements TabbedReader.Factory<GuideRNA> {
		@Override
		public GuideRNA create(String[] rawFields) throws ParseException {
			Annotation a = new BasicAnnotation.Factory().create(rawFields);
			return new GuideRNA(a.trim(0,3), a.getName());
		}
	}
}
