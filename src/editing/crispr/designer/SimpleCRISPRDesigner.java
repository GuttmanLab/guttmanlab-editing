package editing.crispr.designer;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import editing.crispr.GuideRNA;
import editing.crispr.predicate.GuideDoesNotTargetSequences;
import general.CommandLineParser;

import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Gene;
import nextgen.core.annotation.Annotation.Strand;

/**
 * Design single guide RNAs against intervals of interest
 * @author prussell
 *
 */
public class SimpleCRISPRDesigner {
	
	protected SingleGuideRNAFinder finder;
	protected Map<String, Collection<Gene>> intervals; // The intervals to target
	private static Logger logger = Logger.getLogger(SimpleCRISPRDesigner.class.getName());
	protected int numPerInterval; // Number of guides to get per interval
	protected String output; // Output file prefix
	
	/**
	 * A format for outputting a guide RNA
	 * @author prussell
	 *
	 */
	protected interface OutputFormat {
		
		/**
		 * @param parentName Parent annotation name
		 * @param guide sgRNA
		 * @return The information formatted in a string
		 */
		public String format(String parentName, GuideRNA guide);
		
		/**
		 * @return File header or null if no header
		 */
		public String header();
		
	}
	
	/**
	 * A basic table format for guide RNA info
	 * @author prussell
	 *
	 */
	protected class SimpleTableFormat implements OutputFormat {

		public SimpleTableFormat() {}
		
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
	protected class BedFormat implements OutputFormat {

		public BedFormat() {}
		
		@Override
		public String format(String parentName, GuideRNA guide) {
			guide.setName(parentName + ":" + guide.getName());
			return guide.toBED();
		}

		@Override
		public String header() {
			return null;
		}
		
	}
	
	/**
	 * @param guideFinder A guide RNA finder object
	 * @param bedIntervals Bed file of intervals to target
	 * @throws IOException
	 */
	public SimpleCRISPRDesigner(SingleGuideRNAFinder guideFinder, String bedIntervals) throws IOException {
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
	protected void findAndWriteGuides(Map<String, OutputFormat> outFilesAndFormats, int numGuidesPerInterval) throws IOException, InterruptedException {
		logger.info("");
		logger.info("Finding valid guide RNAs and writing output...");
		Map<String, FileWriter> writers = new HashMap<String, FileWriter>();
		for(String file : outFilesAndFormats.keySet()) {
			writers.put(file, new FileWriter(file));
			OutputFormat format = outFilesAndFormats.get(file);
			String header = format.header();
			if(header != null) {
				writers.get(file).write(header + "\n");
			}
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
	 * Get a command line parser object for this class
	 * @return Command line parser object
	 */
	protected static CommandLineParser getCommandLineParser() {
		CommandLineParser p = new CommandLineParser();
		
		p.addStringArg("-g", "Reference sequence fasta e.g. genome or gene sequences", true);
		p.addStringArg("-r", "Bed file of single intervals to target within reference sequences. Omit if targeting entire sequences.", false, null);
		p.addIntArg("-n", "Number of guides to get per sequence/interval. Omit if getting all.", false, -1);
		
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
				
		return p;
	}
	
	/**
	 * Instantiate with a command line
	 * @param parser Command line parser object
	 * @param args Command line arguments
	 * @throws IOException
	 */
	public SimpleCRISPRDesigner(CommandLineParser parser, String[] args) throws IOException {
		
		parser.parse(args);
		
		String genomeFasta = parser.getStringArg("-g");
		String regionBed = parser.getStringArg("-r");
		int numInterval = parser.getIntArg("-n");
		
		boolean efficacyFilter = parser.getBooleanArg("-ge");
		double maxEfficacyScore = parser.getDoubleArg("-gem");
		
		boolean offTargetFilter = parser.getBooleanArg("-got");
		int minOffTargetScore = parser.getIntArg("-gotm");
		String offTargetBitsFile = parser.getStringArg("-gotb");
		
		boolean otherSeqsFilter = parser.getBooleanArg("-os");
		String otherSeqsFasta = parser.getStringArg("-osf");
		int maxMatch16 = parser.getIntArg("-osm16");
		int maxMatch4 = parser.getIntArg("-osm4");
		boolean includeNAG = parser.getBooleanArg("-osnag");
		
		
		String out = parser.getStringArg("-o");
		
		if(parser.getBooleanArg("-d")) {
			logger.setLevel(Level.DEBUG);
			SingleGuideRNAFinder.logger.setLevel(Level.DEBUG);
			GuideDoesNotTargetSequences.logger.setLevel(Level.DEBUG);
		}
		
		SingleGuideRNAFinder gfinder = new SingleGuideRNAFinder(genomeFasta);
		if(efficacyFilter) {
			gfinder.addEfficacyFilter(maxEfficacyScore);
		}
		if(offTargetFilter) {
			if(offTargetBitsFile == null) {
				throw new IllegalArgumentException("If adding off target filter, must provide off target bits file with -gotb");
			}
			gfinder.addOffTargetFilter(offTargetBitsFile, minOffTargetScore);
		}
		if(otherSeqsFilter) {
			if(otherSeqsFasta == null) {
				throw new IllegalArgumentException("If adding off target filter, must provide fasta file of other sequences with -osf");
			}
			gfinder.addGuideDoesNotTargetOtherSeqs(otherSeqsFasta, maxMatch16, maxMatch4, includeNAG);
		}
		
		finder = gfinder;
		if(regionBed != null) {
			intervals = BEDFileParser.loadDataByChr(regionBed);
		} else {
			intervals = new TreeMap<String, Collection<Gene>>();
			Map<String, Sequence> refSeqs = finder.getReferenceSequences();
			for(String seqName : refSeqs.keySet()) {
				int len = refSeqs.get(seqName).getLength();
				Collection<Gene> geneAsList = new ArrayList<Gene>();
				geneAsList.add(new Gene(seqName, 0, len, seqName, Strand.POSITIVE));
				intervals.put(seqName, geneAsList);
			}
		}
		
		output = out;
		numPerInterval = numInterval;
		
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public static void main(String[] args) throws IOException, InterruptedException {
		
		CommandLineParser parser = getCommandLineParser();
		SimpleCRISPRDesigner designer = new SimpleCRISPRDesigner(parser, args);
		
		Map<String, OutputFormat> outFiles = new HashMap<String, OutputFormat>();
		outFiles.put(designer.output + ".bed", designer.new BedFormat());
		outFiles.put(designer.output + ".out", designer.new SimpleTableFormat());
		
		designer.findAndWriteGuides(outFiles, designer.numPerInterval);
		
		logger.info("");
		logger.info("All done.");

	}

}
