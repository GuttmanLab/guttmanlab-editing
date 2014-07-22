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

	public static Logger logger = Logger.getLogger(GuideRNA.class.getName());
	private Sequence sequence20, sequence23;
	private Gene target = null;

	
	/**
	 * @param targetGene Target gene
	 * @param chromosome Chromosome
	 * @param start Start position of 20mer
	 * @param end Position after last position of 20mer
	 * @param allowNAG Allow NAG as valid PAM sequence
	 * @param orientation Strand
	 */
	public GuideRNA(Gene targetGene, Sequence chromosome, int start, int end, Strand orientation, boolean allowNAG) {
		super(chromosome.getId(), start, end, orientation);
		validateStartEnd(start, end);
		validateStrand(orientation);
		setNameFromTarget();
		target = targetGene;
		
		int start20;
		int start23;
		int end20;
		int end23;
		Strand strand;
		
		// Validate
		strand = orientation;
		if(strand.equals(Strand.POSITIVE)) {
			start20 = getStart();
			start23 = start20;
			end20 = getEnd();
			end23 = end20 + 3;
		} else {
			start20 = getStart();
			start23 = getStart() - 3;
			end20 = getEnd();
			end23 = getEnd();
		}
		Sequence shortSeq = chromosome.getSubSequence("", start20, end20);
		shortSeq.setSequenceBases(shortSeq.getSequenceBases().toUpperCase());
		Sequence longSeq = chromosome.getSubSequence("", start23, end23);
		longSeq.setSequenceBases(longSeq.getSequenceBases().toUpperCase());
		sequence20 = strand.equals(Strand.POSITIVE) ? shortSeq : Sequence.reverseSequence(shortSeq);
		sequence23 = strand.equals(Strand.POSITIVE) ? longSeq : Sequence.reverseSequence(longSeq);
		//logger.debug(name + "\tsequence20\t" + sequence20.getSequenceBases() + "\tsequence23\t" + sequence23.getSequenceBases());
		validateSequence20(sequence20);
		validateSequence23(sequence23, allowNAG);
	}

	
	public GuideRNA(Annotation a, String sequence) {
		super(a);
		validateStartEnd(getStart(), getEnd());
		validateStrand(getOrientation());
		setNameFromTarget();
		sequence20 = new Sequence(getName());
		sequence20.setSequenceBases(sequence.substring(0,20));
		sequence20.setForwardStrand(getOrientation() == Strand.POSITIVE);
		
		sequence23 = new Sequence(getName());
		sequence23.setSequenceBases(sequence);
		sequence23.setForwardStrand(getOrientation() == Strand.POSITIVE);
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
		return sequence20.getSequenceBases();
	}
	
	public Sequence getSequence() {
		return sequence20;
	}
	
	public Sequence getSequenceWithPAM() {
		return sequence23;
	}
	
	
	public Gene getTargetGene() {
		return target;
	}

	/**
	 * Find all guide RNAs on either strand within the window
	 * @param chr Chromosome
	 * @param start Start position of window to look in
	 * @param end Position after last position of window
	 * @return All guide RNAs followed by PAM sequence whose 20nt sequence is fully contained in the window
	 */
	public static Collection<GuideRNA> findAll(Sequence chr, int start, int end, Gene targetGene) {
		return findAll(chr, start, end, targetGene, false);
	}
	
	/**
	 * Find all guide RNAs on either strand within the window
	 * @param chr Chromosome
	 * @param start Start position of window to look in
	 * @param end Position after last position of window
	 * @param includeNAG Include NAG as a valid PAM sequence
	 * @return All guide RNAs followed by PAM sequence whose 20nt sequence is fully contained in the window
	 */
	public static Collection<GuideRNA> findAll(Sequence chr, int start, int end, Gene targetGene, boolean includeNAG) {
		//logger.debug("");
		Collection<GuideRNA> rtrn = new ArrayList<GuideRNA>();
		Collection<Annotation> pams = findAllNGGs(chr, start, end);
		if(includeNAG) {
			pams.addAll(findAllNAGs(chr, start, end));
		}
		for(Annotation pam : pams) {
			GuideRNA g = adjacentGuideRNA(chr, start, end, pam, targetGene, includeNAG);
			if(g != null) {
				//logger.debug("Added " + g.toString());
				rtrn.add(g);
			}
		}
		return rtrn;
	}
	
	public String toString() {
		return getReferenceName() + ":" + getStart() + "-" + getEnd() + ":" + getStrand().toString() + ":" + sequence20.getSequenceBases();
	}
	
	/**
	 * Get the correctly oriented guide RNA that is fully contained in the window and ends with the PAM specified
	 * @param windowChr Window chromosome
	 * @param windowStart Window start
	 * @param windowEnd Position after last position of window
	 * @param pam Oriented PAM location
	 * @param targetGene Target gene
	 * @param allowNAG Include NAG as a valid PAM sequence
	 * @return The 20nt guide RNA or null if not fully contained in window
	 */
	private static GuideRNA adjacentGuideRNA(Sequence windowChr, int windowStart, int windowEnd, Annotation pam, Gene targetGene, boolean allowNAG) {
		validatePAM(windowChr, pam, allowNAG);
		//logger.debug("Getting guide RNA adjacent to " + pam.toUCSC() + ":" + pam.getOrientation().toString());
		if(pam.getOrientation().equals(Strand.POSITIVE)) {
			if(pam.getStart() - windowStart < 20) {
				// Guide RNA cannot be fully contained in window
				logger.debug("Guide RNA neighboring " + pam.toUCSC() +" not fully contained in " + windowChr.getId() + ":" + windowStart + "-" + windowEnd);
				return null;
			}
			GuideRNA rtrn = new GuideRNA(targetGene, windowChr, pam.getStart() - 20, pam.getStart(), Strand.POSITIVE, allowNAG);
			//logger.debug("Guide RNA neighboring " + pam.toUCSC() + " is " + rtrn.toString());
			return rtrn;
		}
		if(pam.getOrientation().equals(Strand.NEGATIVE)) {
			if(pam.getEnd() + 20 > windowEnd) {
				// Guide RNA cannot be fully contained in window
				logger.debug("Guide RNA neighboring " + pam.toUCSC() +" not fully contained in " + windowChr.getId() + ":" + windowStart + "-" + windowEnd);
				return null;
			}
			GuideRNA rtrn = new GuideRNA(targetGene, windowChr, pam.getEnd(), pam.getEnd() + 20, Strand.NEGATIVE, allowNAG);
			//logger.debug("Guide RNA neighboring " + pam.toUCSC() + " is " + rtrn.toString());
			return rtrn;
		}
		throw new IllegalArgumentException("Strand must be known");
	}
	
	private static void validatePAM(Sequence chr, Annotation pam, boolean allowNAG) {
		if(pam.numBlocks() != 1) {
			throw new IllegalArgumentException("Must have one block");
		}
		Sequence seq = chr.getSubsequence(pam);
		String substr = seq.getSequenceBases().substring(1).toUpperCase();
		if(!allowNAG && !substr.equals("GG")) {
			throw new IllegalArgumentException("Sequence must be NGG. Is " + seq.getSequenceBases());
		}
		if(allowNAG && !(substr.equals("GG") || substr.equals("AG"))) {
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
		if(end - start != 20) {
			throw new IllegalArgumentException("end - start must equal 20");
		}
	}
	
	private static void validateStrand(Strand orientation) {
		if(orientation.equals(Strand.UNKNOWN)) {
			throw new IllegalArgumentException("Strand cannot be unknown");
		}
	}
	
	private static void validateSequence20(Sequence sequence) {
		if(sequence.getLength() != 20) {
			throw new IllegalArgumentException("Sequence length must be 20 " + sequence.getSequenceBases());
		}
	}
	
	private static void validateSequence23(Sequence sequence, boolean allowNAG) {
		if(sequence.getLength() != 23) {
			throw new IllegalArgumentException("Sequence length must be 23: " + sequence.getSequenceBases());
		}
		String substr = sequence.getSequenceBases().substring(21, 23);
		if(!allowNAG && !substr.equals("GG")) {
			throw new IllegalArgumentException("Sequence must end in GG: " + sequence.getSequenceBases());
		}
		if(allowNAG && !(substr.equals("GG") || substr.equals("AG"))) {
			throw new IllegalArgumentException("Sequence must end in GG or AG: " + sequence.getSequenceBases());
		}
	}

	public int hashCode() {
		String s = toBED() + getName() + sequence20.getSequenceBases();
		if (target != null) s = s + target.toBED();
		return s.hashCode();
	}
	
	public boolean equals(Object o) {
		if(!o.getClass().equals(getClass())) return false;
		GuideRNA g = (GuideRNA)o;
		return hashCode() == g.hashCode();
	}

	
	public String toBedWithSequence() {
		return toBED() + "\t" + sequence20.getSequenceBases();
	}

	/**
	 * @author engreitz
	 * For files stored in BED12 plus the 23-mer sequence in column 13
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
	 * For files stored in BED6 with 23-mer sequence in the name column
	 */
	public static class FactoryBED6 implements TabbedReader.Factory<GuideRNA> {
		@Override
		public GuideRNA create(String[] rawFields) throws ParseException {
			Annotation a = new BasicAnnotation.Factory().create(rawFields);
			return new GuideRNA(a.trim(0,3), a.getName());
		}
	}
}
