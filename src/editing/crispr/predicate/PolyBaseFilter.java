package editing.crispr.predicate;

import editing.crispr.GuideRNA;


public class PolyBaseFilter {

	public String name = "PolyBase";
	private String basesToFilter;
	private int cutoff;
	private int repeatLength;
	
	/**
	 * 
	 */
	public PolyBaseFilter() {}

	public PolyBaseFilter(String basesToFilter, int cutoff, int repeatLength) {
		this.basesToFilter = basesToFilter;
		this.cutoff = cutoff;
		this.repeatLength = repeatLength;
	}
	
	public String name() {
		return "poly_base_filter";
	}

	public boolean rejectSequence(String s) {
		char[] seq = s.toUpperCase().toCharArray();
		for (char c : basesToFilter.toCharArray()) {
			
			int charMatchesInWindow = 0;
			for (int i = 0; i < seq.length; i++) {
				if (seq[i] == c) {
					charMatchesInWindow++;
				}
				
				if (i >= repeatLength) {
					if (seq[i-repeatLength] == c) {
						charMatchesInWindow--;
					}
				}
				
				if (i >= repeatLength - 1) {
					if (charMatchesInWindow >= cutoff) {
						return true;
					}
				}
			}
		}
		return false;
	}
	
	
	
	@Override
	public String toString() {
		return "poly_base_filter";
	}
	

	public boolean evaluate(GuideRNA g) {
		String seq = g.getSequenceString();
		return !rejectSequence(seq);
	}
	
	public String getPredicateName() {
		return name;
	}

}
