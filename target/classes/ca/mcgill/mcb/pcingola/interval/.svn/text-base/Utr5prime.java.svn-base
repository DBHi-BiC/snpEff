package ca.mcgill.mcb.pcingola.interval;

import java.util.Collections;
import java.util.List;

import ca.mcgill.mcb.pcingola.codons.CodonTable;
import ca.mcgill.mcb.pcingola.codons.CodonTables;
import ca.mcgill.mcb.pcingola.interval.SeqChange.ChangeType;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;

/**
 * Interval for a UTR (5 prime UTR and 3 prime UTR
 * 
 * @author pcingola
 *
 */
public class Utr5prime extends Utr {

	private static final long serialVersionUID = 3710420226746056364L;

	public Utr5prime() {
		super();
		type = EffectType.UTR_5_PRIME;
	}

	public Utr5prime(Exon parent, int start, int end, int strand, String id) {
		super(parent, start, end, strand, id);
		type = EffectType.UTR_5_PRIME;
	}

	@Override
	public boolean isUtr3prime() {
		return false;
	}

	@Override
	public boolean isUtr5prime() {
		return true;
	}

	@Override
	public List<ChangeEffect> seqChangeEffect(SeqChange seqChange, ChangeEffect changeEffect) {
		if (seqChange.includes(this) && (seqChange.getChangeType() == ChangeType.DEL)) {
			changeEffect.set(this, EffectType.UTR_5_DELETED, ""); // A UTR was removed entirely
			return changeEffect.newList();
		}

		Transcript tint = (Transcript) findParent(Transcript.class);

		String utrDistStr = utrDistance(seqChange, tint);
		String gained = startGained(seqChange, tint);

		if (gained.length() > 0) changeEffect.set(this, EffectType.START_GAINED, gained + ", " + EffectType.UTR_5_PRIME + ": " + utrDistStr);
		else changeEffect.set(this, type, utrDistStr);

		Exon exon = (Exon) findParent(Exon.class);
		//changeEffect.setExon(exon);
		if (exon != null) exon.check(seqChange, changeEffect); // Check that base matches the expected one

		return changeEffect.newList();
	}

	/**
	 * Is a new start codon produced?
	 * @param chars
	 * @param pos
	 * @return New start codon (or empty string if there is no new start codon)
	 */
	String startGained(char[] chars, int pos) {
		CodonTable ctable = CodonTables.getInstance().getTable(getGenome(), getChromosomeName());

		// Analyze all frames
		for (int i = Math.max(0, pos - 2); (i <= pos) && ((i + 2) < chars.length); i++) {
			String codon = "" + chars[i] + chars[i + 1] + chars[i + 2];

			if (ctable.isStart(codon)) return codon.toUpperCase(); // This frame has a start codon?
		}
		return "";
	}

	/**
	 * Did we gain a start codon in this 5'UTR interval?
	 * @param seqChange
	 * @return A new start codon (if gained)
	 */
	String startGained(SeqChange seqChange, Transcript tint) {
		if (!seqChange.isSnp()) return ""; // FIXME: Only SNPs supported! 

		// Get UTRs and sort them
		List<Utr5prime> utrs = tint.get5primeUtrs();
		if (strand >= 0) Collections.sort(utrs, new IntervalComparatorByStart()); // Sort by start position 
		else Collections.sort(utrs, new IntervalComparatorByEnd(true)); // Sort by end position (reversed) 

		// Create UTR sequence
		StringBuffer sb = new StringBuffer();
		for (Utr5prime utr : utrs) {
			Exon eint = (Exon) utr.getParent();
			String utrSeq = eint.getSequence().substring(0, utr.getEnd() - utr.getStart() + 1);
			sb.append(utrSeq);
		}

		// Calculate SNP position relative to UTRs
		int pos = seqChange.distanceFrom(utrs, strand < 0);

		// Change base at SNP position
		char[] chars = sb.toString().toCharArray();
		char snpBase = seqChange.netChange(this).charAt(0);
		chars[pos] = snpBase;

		// Do we gain a new start codon?
		return startGained(chars, pos);
	}

	/**
	 * Calculate distance from the end of 5'UTRs
	 * 
	 * @param seqChange
	 * @param utr
	 * @return
	 */
	@Override
	String utrDistance(SeqChange seqChange, Transcript tint) {
		List<Utr5prime> utrs = tint.get5primeUtrs();
		boolean fromEnd = !(strand < 0); // We want distance from begining of transcript (TSS = End of 5'UTR)
		int dist = seqChange.distanceFrom(utrs, fromEnd) + 1;
		return dist + " bases from TSS";
	}

}
