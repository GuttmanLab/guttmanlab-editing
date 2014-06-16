package editing.crispr.predicate;


import org.apache.commons.collections15.Predicate;

import editing.crispr.NickingGuideRNAPair;

public interface GuideRNAPairPredicate extends Predicate<NickingGuideRNAPair> {
	
	/**
	 * Get the name of this predicate
	 * @return Predicate name
	 */
	public String getPredicateName();
	
	/**
	 * Get a short explanation (no spaces) of why the predicate evaluates to false, e.g. for inclusion in a bed name
	 * @return Short string explanation of false value
	 */
	public String getShortFailureMessage(NickingGuideRNAPair g);
	
}
