package ca.mcgill.mcb.pcingola.interval;

import java.util.List;

import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;

/**
 * This is a custom interval (i.e. intervals provided by the user)
 * 
 * @author pcingola
 */
public class Custom extends Marker {

	private static final long serialVersionUID = -6843535415295857726L;

	public Custom(Marker parent, int start, int end, int strand, String id) {
		super(parent, start, end, strand, id);
		type = EffectType.CUSTOM;
	}

	@Override
	public List<ChangeEffect> seqChangeEffect(SeqChange snp, ChangeEffect changeEffec) {
		if (!intersects(snp)) return ChangeEffect.emptyResults(); // Sanity check
		changeEffec.set(this, EffectType.CUSTOM, id);
		return changeEffec.newList();
	}

}
