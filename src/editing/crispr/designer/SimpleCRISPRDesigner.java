package editing.crispr.designer;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeSet;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import editing.crispr.GuideRNA;
import editing.crispr.predicate.GuideDoesNotTargetSequences;

import broad.core.parser.CommandLineParser;
import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.Gene;

/**
 * Design singleton guide RNAs against intervals of interest
 * @author prussell
 *
 */
public class SimpleCRISPRDesigner {
	
	private SingleGuideRNAFinder finder;
	private Map<String, Collection<Gene>> intervals;
	private static Logger logger = Logger.getLogger(SimpleCRISPRDesigner.class.getName());
	
	/**
	 * A format for outputting a guide RNA
	 * @author prussell
	 *
	 */
	private interface OutputFormat {
		public String format(String parentName, GuideRNA guide);
	}
	
	/**
	 * A basic table format for guide RNA info
	 * @author prussell
	 *
	 */
	private class SimpleTableFormat implements OutputFormat {

		public SimpleTableFormat() {}
		
		@SuppressWarnings("unused")
		public String header() {
			String rtrn = "parent_name\t";
			rtrn += "chr\t";
			rtrn += "start\t";
			rtrn += "end\t";
			rtrn += "seq\t";
			rtrn += "seq_pam\t";
			rtrn += "seqRC\t";
			rtrn += "seq_pam_RC";
			return rtrn;
		}
		
		@Override
		public String format(String parentName, GuideRNA guide) {
			String seq = guide.getSequenceString();
			String seqp = guide.getSequenceStringWithPAM();
			
			String rtrn = parentName + "\t";
			rtrn += guide.getChr() + "\t";
			rtrn += guide.getStart() + "\t";
			rtrn += guide.getEnd() + "\t";
			rtrn += seq + "\t";
			rtrn += seqp + "\t";
			rtrn += Sequence.reverseSequence(seq) + "\t";
			rtrn += Sequence.reverseSequence(seqp);
			return rtrn;
		}
		
	}
	
	/**
	 * Bed format for a guide RNA
	 * @author prussell
	 *
	 */
	private class BedFormat implements OutputFormat {

		public BedFormat() {}
		
		@Override
		public String format(String parentName, GuideRNA guide) {
			guide.setName(parentName + ":" + guide.getName());
			return guide.toBED();
		}
		
	}
	
	/**
	 * @param guideFinder A guide RNA finder object
	 * @param bedIntervals Bed file of intervals to target
	 * @throws IOException
	 */
	private SimpleCRISPRDesigner(SingleGuideRNAFinder guideFinder, String bedIntervals) throws IOException {
		finder = guideFinder;
		intervals = BEDFileParser.loadDataByChr(bedIntervals);
	}
	
	/**
	 * Find valid guide RNAs and write to a file
	 * @param outFilesAndFormats Output file names with format for each file name
	 * @throws IOException
	 * @throws InterruptedException
	 */
	@SuppressWarnings("unused")
	private void findAndWriteGuides(Map<String, OutputFormat> outFilesAndFormats) throws IOException, InterruptedException {
		findAndWriteGuides(outFilesAndFormats, -1);
	}
	
	/**
	 * Find valid guide RNAs and write to a file
	 * @param outFilesAndFormats Output file names with format for each file name
	 * @param numGuidesPerInterval Number of guides to get per interval, or pass a negative number if getting all
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private void findAndWriteGuides(Map<String, OutputFormat> outFilesAndFormats, int numGuidesPerInterval) throws IOException, InterruptedException {
		logger.info("");
		logger.info("Finding valid guide RNAs and writing output...");
		Map<String, FileWriter> writers = new HashMap<String, FileWriter>();
		for(String file : outFilesAndFormats.keySet()) {
			writers.put(file, new FileWriter(file));
		}
		for(String chr : intervals.keySet()) {
			logger.info(chr);
			for(Gene region : intervals.get(chr)) {
				// Make sure region only has one block
				if(region.numBlocks() > 1) {
					for(FileWriter w : writers.values()) {
						w.close();
					}
					throw new IllegalArgumentException("Region must have only one block: " + region.toBED());
				}
				Collection<GuideRNA> guides = new ArrayList<GuideRNA>();
				if(numGuidesPerInterval >= 0) {
					guides.addAll(finder.getFilteredGuideRNAs(region.getChr(), region.getStart(), region.getEnd(), numGuidesPerInterval));
				} else {
					guides.addAll(finder.getFilteredGuideRNAs(region.getChr(), region.getStart(), region.getEnd()));
				}
				for(GuideRNA guide : guides) {
					for(String file : outFilesAndFormats.keySet()) {
						FileWriter w = writers.get(file);
						w.write(outFilesAndFormats.get(file).format(region.getName(), guide) + "\n");						
					}
				}
			}
		}
		for(FileWriter w : writers.values()) {
			w.close();
		}
	}
	
	/**
	 * Only get a certain number of guide RNAs close to the 5' end of a group of annotations
	 * Annotations are grouped by same name
	 * @param outFilesAndFormats Output file names with format for each file name
	 * @param numPerNameGroup Number of guides to get from the 5' end of each group
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private void findAndWrite5PrimeGuidesForNameGroup(Map<String, OutputFormat> outFilesAndFormats, int numPerNameGroup) throws IOException, InterruptedException {
		logger.info("");
		logger.info("Writing " + numPerNameGroup + " guide RNAs from the 5' end of each group of intervals sharing the same name...");
		Map<String, FileWriter> writers = new HashMap<String, FileWriter>();
		for(String file : outFilesAndFormats.keySet()) {
			writers.put(file, new FileWriter(file));
		}
		Map<String, TreeSet<Gene>> intervalsByName = new HashMap<String, TreeSet<Gene>>();
		// Organize the intervals by grouping ones with same name
		for(String chr : intervals.keySet()) {
			for(Gene interval : intervals.get(chr)) {
				String name = interval.getName();
				if(!intervalsByName.containsKey(name)) {
					intervalsByName.put(name, new TreeSet<Gene>());
				}
				intervalsByName.get(name).add(interval);
			}
		}
		// Save strand for each name group
		Map<String, Strand> strands = new HashMap<String, Strand>();
		for(String name : intervalsByName.keySet()) {
			Strand strand = intervalsByName.get(name).iterator().next().getOrientation();
			if(strand.equals(Strand.UNKNOWN)) {
				throw new IllegalArgumentException("Strand must be known");
			}
			for(Gene interval : intervalsByName.get(name)) {
				if(!interval.getOrientation().equals(strand)) {
					throw new IllegalArgumentException("All intervals with same name must have same strand");
				}
			}
			strands.put(name, strand);
		}
		for(String name : intervalsByName.keySet()) {
			logger.info(name + "\t" + intervalsByName.get(name).iterator().next().getChr() + "\t" + intervalsByName.get(name).size() + " intervals");
			int numFound = 0;
			Iterator<Gene> iter = strands.get(name).equals(Strand.POSITIVE) ? intervalsByName.get(name).iterator() : intervalsByName.get(name).descendingIterator();
			while(iter.hasNext()) {
				if(numFound >= numPerNameGroup) {
					break;
				}
				Gene region = iter.next();
				// Make sure region only has one block
				if(region.numBlocks() > 1) {
					for(FileWriter w : writers.values()) {
						w.close();
					}
					throw new IllegalArgumentException("Region must have only one block: " + region.toBED());
				}
				Collection<GuideRNA> guides = finder.getFilteredGuideRNAs(region.getChr(), region.getStart(), region.getEnd());
				// Put the guide RNAs in order
				TreeSet<GuideRNA> ordered = new TreeSet<GuideRNA>();
				ordered.addAll(guides);
				Iterator<GuideRNA> gIter = strands.get(name).equals(Strand.POSITIVE) ? ordered.iterator() : ordered.descendingIterator();
				// Only get guides until desired number is reached
				while(gIter.hasNext()) {
					GuideRNA guide = gIter.next();
					for(String file : outFilesAndFormats.keySet()) {
						FileWriter w = writers.get(file);
						w.write(outFilesAndFormats.get(file).format(name, guide) + "\n");					
					}
					numFound++;
					if(numFound >= numPerNameGroup) {
						break;
					}
				}
			}
			if(numFound < numPerNameGroup) {
				logger.warn("Only found " + numFound + " valid guide RNAs for " + name + ".");
			}
		}
		for(FileWriter w : writers.values()) {
			w.close();
		}
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public static void main(String[] args) throws IOException, InterruptedException {
		
		CommandLineParser p = new CommandLineParser();
		
		p.addStringArg("-g", "Genome fasta", true);
		p.addStringArg("-r", "Bed file of single intervals to target", true);
		p.addIntArg("-n", "Number of guides to get per interval, omit if getting all", false, -1);
		
		p.addBooleanArg("-ge", "Apply guide efficacy filter", false, false);
		p.addDoubleArg("-gem", "Max guide efficacy if using filter", false, SingleGuideRNAFinder.MAX_GUIDE_EFFICACY_SCORE);
		
		p.addBooleanArg("-got", "Apply guide off target filter (requires off target bits file)", false, false);
		p.addIntArg("-gotm", "Min guide off target score", false, SingleGuideRNAFinder.MIN_OFF_TARGET_SCORE);
		p.addStringArg("-gotb", "Off target bits file", false);
		
		p.addBooleanArg("-os", "Filter guides that potentially target a sequence in another set", false, false);
		p.addStringArg("-osf", "Fasta file of other sequences to check for potential targetting", false);
		p.addIntArg("-osm16", "For other sequences filter, max number of matches over first 16bp locus with a PAM sequence", false, 12);
		p.addIntArg("-osm4", "For other sequences filter, max number of matches over the last 4bp for locus with a PAM sequence", false, 0);
		p.addBooleanArg("-osnag", "For other sequences filter, count NAG as a possible PAM sequence", false, false);
		
		p.addStringArg("-o", "Output prefix", true);
		p.addBooleanArg("-d", "Debug logging", false, false);
		
		p.addBooleanArg("-5p", "Only get a certain number of guide RNAs from the 5' end of each group of intervals sharing the same name", false, false);
		p.addIntArg("-5pn", "If using -5p, number to get for each name group", false, 2);
		
		p.parse(args);
		
		String genomeFasta = p.getStringArg("-g");
		String regionBed = p.getStringArg("-r");
		int numPerInterval = p.getIntArg("-n");
		
		boolean efficacyFilter = p.getBooleanArg("-ge");
		double maxEfficacyScore = p.getDoubleArg("-gem");
		
		boolean offTargetFilter = p.getBooleanArg("-got");
		int minOffTargetScore = p.getIntArg("-gotm");
		String offTargetBitsFile = p.getStringArg("-gotb");
		
		boolean otherSeqsFilter = p.getBooleanArg("-os");
		String otherSeqsFasta = p.getStringArg("-osf");
		int maxMatch16 = p.getIntArg("-osm16");
		int maxMatch4 = p.getIntArg("-osm4");
		boolean includeNAG = p.getBooleanArg("-osnag");
		
		boolean fivePrime = p.getBooleanArg("-5p");
		int fivePrimeNum = p.getIntArg("-5pn");
		
		String output = p.getStringArg("-o");
		
		if(p.getBooleanArg("-d")) {
			logger.setLevel(Level.DEBUG);
			SingleGuideRNAFinder.logger.setLevel(Level.DEBUG);
			GuideDoesNotTargetSequences.logger.setLevel(Level.DEBUG);
		}
		
		SingleGuideRNAFinder finder = new SingleGuideRNAFinder(genomeFasta);
		if(efficacyFilter) {
			finder.addEfficacyFilter(maxEfficacyScore);
		}
		if(offTargetFilter) {
			if(offTargetBitsFile == null) {
				throw new IllegalArgumentException("If adding off target filter, must provide off target bits file with -gotb");
			}
			finder.addOffTargetFilter(offTargetBitsFile, minOffTargetScore);
		}
		if(otherSeqsFilter) {
			if(otherSeqsFasta == null) {
				throw new IllegalArgumentException("If adding off target filter, must provide fasta file of other sequences with -osf");
			}
			finder.addGuideDoesNotTargetOtherSeqs(otherSeqsFasta, maxMatch16, maxMatch4, includeNAG);
		}
		
		SimpleCRISPRDesigner designer = new SimpleCRISPRDesigner(finder, regionBed);
		
		Map<String, OutputFormat> outFiles = new HashMap<String, OutputFormat>();
		outFiles.put(output + ".bed", designer.new BedFormat());
		outFiles.put(output + ".out", designer.new SimpleTableFormat());
		
		if(fivePrime) {
			designer.findAndWrite5PrimeGuidesForNameGroup(outFiles, fivePrimeNum);
		} else {
			designer.findAndWriteGuides(outFiles, numPerInterval);
		}
		
		logger.info("");
		logger.info("All done.");

	}

}
