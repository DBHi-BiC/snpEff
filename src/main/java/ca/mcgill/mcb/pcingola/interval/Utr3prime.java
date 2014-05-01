package ca.mcgill.mcb.pcingola.interval;

import ca.mcgill.mcb.pcingola.interval.SeqChange.ChangeType;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffects;
import chop.cbmi.leipzig.interval.Utr3primeChange;

/**
 * Interval for a UTR (5 prime UTR and 3 prime UTR
 * 
 * @author pcingola
 *
 */
public class Utr3prime extends Utr {

	private static final long serialVersionUID = 5688641008301281991L;

	public Utr3prime() {
		super();
		type = EffectType.UTR_3_PRIME;
	}

	public Utr3prime(Exon parent, int start, int end, int strand, String id) {
		super(parent, start, end, strand, id);
		type = EffectType.UTR_3_PRIME;
	}

	@Override
	public boolean isUtr3prime() {
		return true;
	}

	@Override
	public boolean isUtr5prime() {
		return false;
	}

	@Override
	public boolean seqChangeEffect(SeqChange seqChange, ChangeEffects changeEffects) {
		if (!intersects(seqChange)) return false;

		if (seqChange.includes(this) && (seqChange.getChangeType() == ChangeType.DEL)) {
			changeEffects.add(this, EffectType.UTR_3_DELETED, ""); // A UTR was removed entirely
			return true;
		}

		Transcript tr = (Transcript) findParent(Transcript.class);
		int dist = utrDistance(seqChange, tr);
		changeEffects.add(this, type, dist >= 0 ? dist + " bases from CDS" : "");
		if (dist >= 0) changeEffects.setDistance(dist);

        Utr3primeChange utrChange = new Utr3primeChange(seqChange, this);
        changeEffects.add(utrChange.calculate());

		return true;
	}

	/**
	 * Calculate distance from beginning of 3'UTRs
	 * 
	 * @param snp
	 * @param utr
	 * @return
	 */
	@Override
	int utrDistance(SeqChange seqChange, Transcript tr) {
		int cdsEnd = tr.getCdsEnd();
		if (cdsEnd < 0) return -1;

		if (isStrandPlus()) return seqChange.getStart() - cdsEnd;
		return cdsEnd - seqChange.getEnd();
	}

}
