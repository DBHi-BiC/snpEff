package ca.mcgill.mcb.pcingola.snpEffect;

import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Variant;
import ca.mcgill.mcb.pcingola.interval.Transcript;

/**
 * HGSV notation
 *
 * References: http://www.hgvs.org/
 *
 * @author pcingola
 */
public class Hgvs {

	protected ChangeEffect changeEffect;
	Variant seqChange;
	Marker marker;
	Transcript tr;

	public Hgvs(ChangeEffect changeEffect) {
		this.changeEffect = changeEffect;
		seqChange = changeEffect.getSeqChange();
		marker = changeEffect.getMarker();
		tr = changeEffect.getTranscript();
	}

}
