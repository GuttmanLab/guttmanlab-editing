package editing.crispr.designer;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Gene;
import nextgen.core.annotation.Annotation.Strand;

import broad.core.sequence.Sequence;
import editing.crispr.GuideRNA;
import guttmanlab.core.util.CommandLineParser;
import guttmanlab.core.util.StringParser;

/**
 * Design guide RNAs against groups of intervals with the same name, i.e. exons of a gene that are given a common name in the input file
 * @author prussell
 *
 */
public class NameGroupedCRISPRDesigner extends SimpleCRISPRDesigner {

	private static Logger logger = Logger.getLogger(NameGroupedCRISPRDesigner.class.getName());
	private Map<String, TreeSet<Gene>> intervalsByName; // Intervals grouped by their common name
	private Map<String, Strand> strands; // Strand of each interval group
	
	/**
	 * @param guideFinder A guide RNA finder object
	 * @param bedIntervals Bed file of intervals to target
	 * @throws IOException
	 */
	public NameGroupedCRISPRDesigner(SingleGuideRNAFinder guideFinder, String bedIntervals) throws IOException {
		super(guideFinder, bedIntervals);
		storeIntervalsByName();
	}
	
	/**
	 * @param parser Command line parser object
	 * @param args Command line arguments
	 * @throws IOException
	 */
	public NameGroupedCRISPRDesigner(CommandLineParser parser, String[] args) throws IOException {
		super(parser, args);
		storeIntervalsByName();
	}
	
	/**
	 * Store the groupings of intervals with same name
	 */
	private void storeIntervalsByName() {
		intervalsByName = new HashMap<String, TreeSet<Gene>>();
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
		strands = new HashMap<String, Strand>();
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

	}

	/**
	 * Find guide RNAs for each group of intervals and write to files
	 * @param outFilesAndFormats Output files and formats
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private void findAndWriteGuidesByGroup(Map<String, OutputFormat> outFilesAndFormats) throws IOException, InterruptedException {
		findAndWrite5PrimeGuidesForGroup(outFilesAndFormats, Integer.MAX_VALUE);
	}

	/**
	 * Read super groups (groups of groups) from a file
	 * Line format: <group name> <supergroup name>
	 * @param file File to read from
	 * @return Map of supergroup name to collection of group names making up the supergroup
	 * @throws IOException
	 */
	private static Map<String, Collection<String>> readSuperGroupsFromFile(String file) throws IOException {
		BufferedReader r = new BufferedReader(new FileReader(file));
		StringParser s = new StringParser();
		Map<String, Collection<String>> rtrn = new HashMap<String, Collection<String>>();
		while(r.ready()) {
			s.parse(r.readLine());
			String name = s.asString(0);
			String group = s.asString(1);
			if(!rtrn.containsKey(group)) {
				rtrn.put(group, new HashSet<String>());
			}
			rtrn.get(group).add(name);
		}
		r.close();
		return rtrn;
	}
	
	/**
	 * Find guide RNAs that simulatneously target all groups in super groups and write output
	 * @param superGroupsFile File describing super groups. Line format: <group name> <supergroup name>
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private void findAndWriteGuidesTargetingSupergroups(String superGroupsFile) throws IOException, InterruptedException {
		Map<String, Collection<String>> superGroups = readSuperGroupsFromFile(superGroupsFile);
		FileWriter tableWriter = new FileWriter(output + ".out");
		FileWriter bedWriter = new FileWriter(output + ".bed");
		String tableHeader = "super_group\t";
		tableHeader += "guide_sequence\t";
		tableHeader += "guide_sequence_with_PAM\t";
		tableHeader += "guide_sequence_RC\t";
		tableHeader += "guide_sequence_with_PAM_RC\t";
		tableWriter.write(tableHeader + "\n");
		for(String superGroup : superGroups.keySet()) {
			Map<String, Collection<GuideRNA>> sharedGuides = findGuidesTargetingAllGroups(superGroups.get(superGroup));
			for(String withPAM : sharedGuides.keySet()) {
				// Write table
				String noPAM = withPAM.substring(0, 20);
				String noPamRC = Sequence.reverseSequence(noPAM);
				String withPamRC = Sequence.reverseSequence(withPAM);
				String line = superGroup  + "\t";
				line += noPAM + "\t";
				line += withPAM + "\t";
				line += noPamRC + "\t";
				line += withPamRC + "\t";
				tableWriter.write(line + "\n");
				// Write bed
				for(GuideRNA guide : sharedGuides.get(withPAM)) {
					bedWriter.write(guide.toBED() + "\n");
				}
			}
		}
		tableWriter.close();
		bedWriter.close();
	}
	
	/**
	 * Find guide RNAs that simulatneously target all groups in a super group
	 * @param names The names of the groups in the super group
	 * @return A map of guide RNA sequence (including PAM) to the collection of guide RNA objects targeting each group having this sequence
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private Map<String, Collection<GuideRNA>> findGuidesTargetingAllGroups(Collection<String> names) throws IOException, InterruptedException {
		Map<String, Collection<GuideRNA>> rtrn = new HashMap<String, Collection<GuideRNA>>();
		// Check that all the group names are valid
		for(String name : names) {
			if(!intervalsByName.containsKey(name)) {
				throw new IllegalArgumentException("Group name " + name + " not recognized.");
			}
		}
		Map<String, Collection<GuideRNA>> guidesByGroup = new HashMap<String, Collection<GuideRNA>>();
		for(String name : names) {
			guidesByGroup.put(name, findGuidesForGroup(name));
		}
		// Check if any sets are empty
		for(String name : guidesByGroup.keySet()) {
			if(guidesByGroup.get(name).isEmpty()) {
				logger.warn("Group " + name + " has no valid guide RNAs so returning null");
				return null;
			}
		}
		// Identify guides that are in all sets
		Collection<GuideRNA> firstSet = guidesByGroup.values().iterator().next();
		for(GuideRNA guide : firstSet) {
			String seq23 = guide.getSequenceStringWithPAM();
			boolean targetsAll = true;
			Collection<GuideRNA> guidesWithSameSeq = new ArrayList<GuideRNA>();
			for(Collection<GuideRNA> guides : guidesByGroup.values()) {
				boolean found = false;
				for(GuideRNA otherGuide : guides) {
					if(otherGuide.getSequenceStringWithPAM().equals(seq23)) {
						found = true;
						guidesWithSameSeq.add(otherGuide);
						break;
					}
				}
				if(!found) {
					targetsAll = false;
					break;
				}
			}
			if(targetsAll) {
				rtrn.put(seq23, guidesWithSameSeq);
			}
		}
		if(rtrn.isEmpty()) {
			String nameString = "";
			for(String name : names) {
				nameString += name + " ";
			}
			logger.warn("No guides found that target all of " + nameString);
		}
		return rtrn;
	}
	
	/**
	 * Find all guide RNAs targeting each group
	 * @return Map of group name to collection of guides targeting the group
	 * @throws IOException
	 * @throws InterruptedException
	 */
	@SuppressWarnings("unused")
	private Map<String, Collection<GuideRNA>> findGuidesByGroup() throws IOException, InterruptedException {
		return findGuidesByGroup(Integer.MAX_VALUE);
	}
	
	/**
	 * Find all guide RNAs targeting a particular group
	 * @param name Group name
	 * @return Collection of guide RNAs targeting this group
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private Collection<GuideRNA> findGuidesForGroup(String name) throws IOException, InterruptedException {
		return findGuidesForGroup(name, Integer.MAX_VALUE);
	}
	
	/**
	 * Find a particular number of guide RNAs targeting a particular group
	 * Guides are returned in no particular order, no guarantee on which ones will be returned
	 * If there are insufficient guide RNAs, returns fewer
	 * @param name Group name
	 * @param numToGet Number of guide RNAs to get
	 * @return The specified number of guides or fewer if insufficiently many
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private Collection<GuideRNA> findGuidesForGroup(String name, int numToGet) throws IOException, InterruptedException {
		Collection<GuideRNA> rtrn = new ArrayList<GuideRNA>();
		logger.info(name + "\t" + intervalsByName.get(name).iterator().next().getChr() + "\t" + intervalsByName.get(name).size() + " intervals");
		int numFound = 0;
		Iterator<Gene> iter = strands.get(name).equals(Strand.POSITIVE) ? intervalsByName.get(name).iterator() : intervalsByName.get(name).descendingIterator();
		while(iter.hasNext()) {
			if(numFound >= numToGet) {
				break;
			}
			Gene region = iter.next();
			// Make sure region only has one block
			if(region.numBlocks() > 1) {
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
				numFound++;
				rtrn.add(guide);
				if(numFound >= numToGet) {
					break;
				}
			}
		}
		return rtrn;
	}
	
	/**
	 * Find a particular number of guide RNAs targeting each group
	 * Guides are returned in no particular order, no guarantee on which ones will be returned
	 * If there are insufficient guide RNAs for some group, returns fewer
	 * @param numPerGroup Number of guide RNAs to get per group
	 * @return Map of group name to collection of the specified number of guides or fewer if insufficiently many
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private Map<String, Collection<GuideRNA>> findGuidesByGroup(int numPerGroup) throws IOException, InterruptedException {
		Map<String, Collection<GuideRNA>> rtrn = new HashMap<String, Collection<GuideRNA>>();
		for(String name : intervalsByName.keySet()) {
			Collection<GuideRNA> guides = findGuidesForGroup(name, numPerGroup);
			if(guides.size() < numPerGroup) {
				logger.warn("Found " + guides.size() + " valid guide RNAs for " + name + ".");
			}
			rtrn.put(name, guides);
		}
		return rtrn;
	}
	
	/**
	 * Only get a certain number of guide RNAs close to the 5' end of a group of annotations
	 * Annotations are grouped by same name
	 * @param outFilesAndFormats Output file names with format for each file name
	 * @param numPerGroup Number of guides to get from the 5' end of each group
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private void findAndWrite5PrimeGuidesForGroup(Map<String, OutputFormat> outFilesAndFormats, int numPerGroup) throws IOException, InterruptedException {
		logger.info("");
		logger.info("Writing " + numPerGroup + " guide RNAs from the 5' end of each group of intervals sharing the same name...");
		Map<String, FileWriter> writers = new HashMap<String, FileWriter>();
		for(String file : outFilesAndFormats.keySet()) {
			writers.put(file, new FileWriter(file));
		}
		Map<String, Collection<GuideRNA>> guides = findGuidesByGroup(numPerGroup);
		for(String name : guides.keySet()) {
			for(GuideRNA guide : guides.get(name)) {
				for(String file : outFilesAndFormats.keySet()) {
					FileWriter w = writers.get(file);
					w.write(outFilesAndFormats.get(file).format(name, guide) + "\n");					
				}
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
		
		CommandLineParser parser = getCommandLineParser();
		
		// Add options for getting a certain number of 5'-most guides per group
		parser.addBooleanArg("-5p", "Only get a certain number of guide RNAs from the 5' end of each group of intervals sharing the same name", false, false);
		parser.addIntArg("-5pn", "If using -5p, number to get for each name group", false, 2);
		
		// Add options for getting guides that target several groups at once
		parser.addStringArg("-s", "File of supergroups to find guides that target all groups. Line format: group   supergroup", false, null);
		
		NameGroupedCRISPRDesigner designer = new NameGroupedCRISPRDesigner(parser, args);
		
		boolean fivePrime = parser.getBooleanArg("-5p");
		int fivePrimeNum = parser.getIntArg("-5pn");
		String supergroupsFile = parser.getStringArg("-s");
		
		Map<String, OutputFormat> outFiles = new HashMap<String, OutputFormat>();
		outFiles.put(designer.output + ".bed", designer.new BedFormat());
		outFiles.put(designer.output + ".out", designer.new SimpleTableFormat());
		
		if(fivePrime) {
			designer.findAndWrite5PrimeGuidesForGroup(outFiles, fivePrimeNum);
		} else if(supergroupsFile != null) {
			designer.findAndWriteGuidesTargetingSupergroups(supergroupsFile);
		} else {
			designer.findAndWriteGuidesByGroup(outFiles);
		}
		
		logger.info("");
		logger.info("All done.");

	}

}
