package ca.mcgill.mcb.pcingola.interval;

import java.util.List;

import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;

/**
 * Interval for a gene, as well as some other information: exons, utrs, cds, etc.
 * 
 * @author pcingola
 *
 */
public class Downstream extends Marker {

	private static final long serialVersionUID = 1636197649250882952L;

	public Downstream(Transcript parent, int start, int end, int strand, String id) {
		super(parent, start, end, strand, id);
		type = EffectType.DOWNSTREAM;
	}

	/**
	 * Upstream sites are no included in transcript (by definition).
	 */
	@Override
	protected boolean isShowWarningIfParentDoesNotInclude() {
		return false;
	}

	@Override
	public List<ChangeEffect> seqChangeEffect(SeqChange snp, ChangeEffect changeEffect) {
		if (!intersects(snp)) return ChangeEffect.emptyResults(); // Sanity check
		int distance = (parent.isStrandPlus() ? snp.getStart() - start : end - snp.getStart()) + 1;
		changeEffect.set(this, EffectType.DOWNSTREAM, distance + " bases");
		return changeEffect.newList();
	}

}
