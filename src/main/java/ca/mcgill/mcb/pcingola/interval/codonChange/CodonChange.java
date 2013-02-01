package ca.mcgill.mcb.pcingola.interval.codonChange;

import java.util.ArrayList;
import java.util.List;

import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;

/**
 * Analyze codon changes based on a SeqChange and a Transcript
 * 
 * @author pcingola
 */
public class CodonChange {

	public static int SHOW_CODONS_AROUND_CHANGE = 0; // Show this many changes around the codon
	public static boolean showCodonChange = true; // This is disabled in some specific test cases
	public static final int CODON_SIZE = 3; // I'll be extremely surprised if you ever need to change this parameter...

	boolean returnNow = false; // Can we return immediately after calculating the first 'codonChangeSingle()'?
	boolean requireNetCdsChange = false;
	SeqChange seqChange;
	Transcript transcript;
	Exon exon = null;
	ChangeEffect changeEffect;
	int codonNum = -1;
	int codonIndex = -1;
    int txPos = -1; //transcript-relative nt position
	String codonsOld = ""; // Old codons (before change)
	String codonsNew = ""; // New codons (after change)
	String aaOld = ""; // Old amino acids (before change)
	String aaNew = ""; // New amino acids (after change)
	String netCdsChange = "";

	public CodonChange(SeqChange seqChange, Transcript transcript, ChangeEffect changeEffect) {
		this.seqChange = seqChange;
		this.transcript = transcript;
		this.changeEffect = changeEffect;
	}

	/**
	 * Calculate all possible codon changes
	 *  
	 * @param seqChange
	 * @param changeEffect
	 * @return
	 */
	public List<ChangeEffect> calculate() {
		ArrayList<ChangeEffect> changes = new ArrayList<ChangeEffect>();

		// Split each seqChange into it's multiple options
		for (int i = 0; i < seqChange.getChangeOptionCount(); i++) {
			ChangeEffect changeEffectNew = changeEffect.clone(); // Create a copy of this result

			// Create a new SeqChange for this option, calculate codonChange for this seqChangeOption and add result to the list
			SeqChange seqChangeNew = seqChange.getSeqChangeOption(i);
			if (seqChangeNew != null) {
				// Create a specific codon change and calculate changes
				CodonChange codonChange = factory(seqChangeNew, transcript, changeEffectNew);
				changes.addAll(codonChange.codonChange()); // Calculate codon change and add them to the list
			}
		}

		return changes;
	}

	/**
	 * Calculate base number in a cds where 'pos' is
	 * @param pos
	 * @return
	 */
	int cdsBaseNumber(int pos) {
		int cdsbn = transcript.cdsBaseNumber(pos, true);

		// Does not intersect the transcript?
		if (cdsbn < 0) {
			// 'pos' before transcript start
			if (pos <= transcript.getCdsStart()) {
				if (transcript.isStrandPlus()) return 0;
				return transcript.cds().length();
			}

			// 'pos' is after CDS end
			if (transcript.isStrandPlus()) return transcript.cds().length();
			return 0;
		}

		return cdsbn;
	}

	/**
	 * Calculate a list of codon changes
	 * @return
	 */
	List<ChangeEffect> codonChange() {
		ArrayList<ChangeEffect> changes = new ArrayList<ChangeEffect>();
		if (!transcript.intersects(seqChange)) return changes;

		List<Exon> exons = transcript.sortedStrand();

		// Get coding start (after 5 prime UTR)
		int cdsStart = transcript.getCdsStart();

		// We may have to calculate 'netCdsChange', which is the effect on the CDS
		netCdsChange = netCdsChange();

		//---
		// Concatenate all exons
		//---
		int firstCdsBaseInExon = 0; // Where the exon maps to the CDS (i.e. which CDS base number does the first base in this exon maps to).
		for (Exon exon : exons) {
			this.exon = exon;

			if (exon.intersects(seqChange)) {
				int cdsBaseInExon; // cdsBaseInExon: base number relative to the beginning of the coding part of this exon (i.e. excluding 5'UTRs)

				if (transcript.isStrandPlus()) {
					int firstSeqChangeBaseInExon = Math.max(seqChange.getStart(), Math.max(exon.getStart(), cdsStart));
					cdsBaseInExon = firstSeqChangeBaseInExon - Math.max(exon.getStart(), cdsStart);
				} else {
					int lastSeqChangeBaseInExon = Math.min(seqChange.getEnd(), Math.min(exon.getEnd(), cdsStart));
					cdsBaseInExon = Math.min(exon.getEnd(), cdsStart) - lastSeqChangeBaseInExon;
				}

				if (cdsBaseInExon < 0) cdsBaseInExon = 0;

				// Get codon number and index within codon (where seqChage is pointing)
				codonNum = (firstCdsBaseInExon + cdsBaseInExon) / CODON_SIZE;
				codonIndex = (firstCdsBaseInExon + cdsBaseInExon) % CODON_SIZE;

                //txPos is cdsBaseinTranscript
                txPos = (firstCdsBaseInExon + cdsBaseInExon);

				// Use appropriate method to calculate codon change
				boolean hasChanged = false; // Was there any change?
				ChangeEffect changeEffectNew = changeEffect.clone();
				changeEffectNew.setMarker(exon); // It is affecting this exon, so we set the marker
				hasChanged = codonChangeSingle(changeEffectNew, exon);

				// Any change? => Add change to list
				if (hasChanged) {
					codonsAround(seqChange, changeEffectNew, codonNum); // Show codons around change (if required)
					changes.add(changeEffectNew);
				}

				// Can we return immediately?
				if (returnNow) return changes;
			}

			if (transcript.isStrandPlus()) firstCdsBaseInExon += Math.max(0, exon.getEnd() - Math.max(exon.getStart(), cdsStart) + 1);
			else firstCdsBaseInExon += Math.max(0, Math.min(cdsStart, exon.getEnd()) - exon.getStart() + 1);
		}

		return changes;
	}

	/**
	 * Calculate the effect of a single change type: SNP, MNP, INS, DEL
	 * @return
	 */
	boolean codonChangeSingle(ChangeEffect changeEffect, Exon exon) {
		throw new RuntimeException("Unimplemented method for this thype of seqChange: " + seqChange.getType());
	}

	/**
	 * Calculate codons changes
	 * @param seqChange
	 * @param changeEffect
	 * @param codonNum
	 */
	void codonsAround(SeqChange seqChange, ChangeEffect changeEffect, int codonNum) {
		if (SHOW_CODONS_AROUND_CHANGE <= 0) return; // Nothing to do?

		String cdsSeq = transcript.cds();
		int changeSizeInCodons = seqChange.size() / 3;

		// Calculate codon positions
		int codonMinBasePos = Math.max(0, codonNum - SHOW_CODONS_AROUND_CHANGE) * CodonChange.CODON_SIZE;
		int codonStartBasePos = codonNum * CodonChange.CODON_SIZE;
		int codonEndBasePos = Math.min(cdsSeq.length(), codonStartBasePos + (1 + changeSizeInCodons) * CodonChange.CODON_SIZE);
		int codonMaxBasePos = Math.min(cdsSeq.length(), codonEndBasePos + SHOW_CODONS_AROUND_CHANGE * CodonChange.CODON_SIZE);

		// Calculate codons around the seqChange
		String codonsLeft = cdsSeq.substring(codonMinBasePos, codonStartBasePos);
		String codonsRight = cdsSeq.substring(codonEndBasePos, codonMaxBasePos);

		changeEffect.setCodonsAround(codonsLeft, codonsRight);
	}

	/**
	 * Calculate new codons
	 * @return
	 */
	String codonsNew() {
		throw new RuntimeException("Unimplemented method for this thype of CodonChange: " + this.getClass().getSimpleName());
	}

	/**
	 * Calculate old codons
	 * @return
	 */
	public String codonsOld() {
		return codonsOld(1);
	}

	/**
	 * Calculate old codons
	 * @return
	 */
	protected String codonsOld(int numCodons) {
		String cds = transcript.cds();
		String codon = "";

		int start = codonNum * CodonChange.CODON_SIZE;
		int end = start + numCodons * CodonChange.CODON_SIZE;

		int len = cds.length();
		if (start >= cds.length()) start = len;
		if (end >= cds.length()) end = len;

		// Capitalize
		codon = cds.substring(start, end);

		// Codon not multiple of three? Add missing bases as 'N'
		if (codon.length() % 3 == 1) codon += "NN";
		else if (codon.length() % 3 == 2) codon += "N";

		return codon;
	}

	/**
	 * Create a specific codon change for a seqChange
	 * @param seqChange
	 * @param transcript
	 * @param changeEffect
	 * @return
	 */
	CodonChange factory(SeqChange seqChange, Transcript transcript, ChangeEffect changeEffect) {
		if (seqChange.isSnp()) return new CodonChangeSnp(seqChange, transcript, changeEffect);
		if (seqChange.isIns()) return new CodonChangeIns(seqChange, transcript, changeEffect);
		if (seqChange.isDel()) return new CodonChangeDel(seqChange, transcript, changeEffect);
		if (seqChange.isMnp()) return new CodonChangeMnp(seqChange, transcript, changeEffect);
		if (seqChange.isInterval()) return new CodonChangeInterval(seqChange, transcript, changeEffect);
		throw new RuntimeException("Unimplemented factory for SeqChange: " + seqChange);
	}

	/**
	 * We may have to calculate 'netCdsChange', which is the effect on the CDS
	 * Note: A deletion or a MNP might affect several exons
	 * @return
	 */
	protected String netCdsChange() {
		if (!requireNetCdsChange) return "";

		if (seqChange.size() > 1) {
			StringBuilder sb = new StringBuilder();
			for (Exon exon : transcript.sortedStrand())
				sb.append(seqChange.netChange(exon));
			return sb.toString();
		}

		return seqChange.netChange(transcript.getStrand());
	}

}
