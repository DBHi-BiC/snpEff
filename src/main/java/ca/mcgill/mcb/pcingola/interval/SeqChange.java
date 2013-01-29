package ca.mcgill.mcb.pcingola.interval;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;
import ca.mcgill.mcb.pcingola.util.GprSeq;

/**
 * A change in a reference sequence
 * 
 * @author pcingola
 */
public class SeqChange extends Marker {

	public enum ChangeType {
		SNP// Single nucleotide polymorphism (i.e. 1 base is changed)
		, MNP // Multiple nucleotide polymorphism (i.e. several bases are changed)
		, INS // Insertion (i.e. some bases added)
		, DEL // Deletion (some bases removed)
		, MIXED // A mixture of insertion, deletions, SNPs and or MNPs
		, Interval
		// Just analyze interval hits. Not a sequence change (e.g. BED input format)
	}

	private static final long serialVersionUID = -2928105165111400441L;;

	Boolean heterozygous; // Is this an heterozygous change?
	ChangeType changeType; // Change type
	String reference; // Reference (i.e. original bases in the genome)
	String change; // Changed bases
	String[] changeOptions; // Available change options (when multiple changes or heterozygous)  
	double quality; // Quality prediction (negative means 'not available')
	double score = Double.NaN; // Score (e.g. BED files)
	int coverage; // Number of reads covering the position (negative means 'not available')
	boolean imprecise = false; // Imprecise variant (coordinates are not exact (E.g. see section "Encoding Structural Variants in VCF" from VCF spec. 4.1)

	public SeqChange(Marker parent, int start, int end, int strand, String id) {
		super(parent, start, end, strand, id);
		reference = change = "";
		changeType = ChangeType.Interval;
		changeOptions = new String[1];
		changeOptions[0] = "";
	}

	/**
	 * This constructor is used when we only have interval data (e.g. when reading a BED file)
	 * @param parent
	 * @param start
	 * @param end
	 * @param id
	 */
	public SeqChange(Marker parent, int start, int end, String id) {
		super(parent, start, end, 1, id);
		reference = change = "";
		changeType = ChangeType.Interval;
		changeOptions = new String[1];
		changeOptions[0] = "";
	}

	/**
	 * Create a sequence change
	 * @param parent : Chromosome interval
	 * @param position : SNP position
	 * @param referenceStr : Reference base (original base)
	 * @param changeStr : SNP base
	 * @param strand : Strand
	 * @param id : Snp ID
	 * @param quality : Snp Quality (negative if not available)
	 * @param coverage : Snp coverage (negative if not available)
	 * @param filterPass : Filter pass or fail
	 * @param filterExplain : Filter details 
	 */
	public SeqChange(Marker parent, int position, String referenceStr, String changeStr, int strand, String id, double quality, int coverage) {
		super(parent, position, position, strand, id);
		init(parent, position, referenceStr, changeStr, strand, id, quality, coverage);
	}

	/**
	 * Return the change (always in positive strand) 
	 * @return
	 */
	public String change() {
		return strand >= 0 ? change : GprSeq.reverseWc(change);
	}

	/**
	 * Distance from the beginning/end of a list of intervals, until this SNP
	 * @param ints
	 * @return
	 */
	public int distanceFrom(List<? extends Marker> ints, boolean fromEnd) {
		// Create a new list of sorted intervals
		ArrayList<Marker> intsSort = new ArrayList<Marker>();
		intsSort.addAll(ints);
		if (fromEnd) Collections.sort(intsSort, new IntervalComparatorByEnd(true));
		else Collections.sort(intsSort, new IntervalComparatorByStart());

		// Calculate distance
		int len = 0, latest = 0;
		for (Marker i : intsSort) {
			if (fromEnd) {
				if (intersects(i)) return len + (i.getEnd() - start);
				latest = i.getStart();
			} else {
				if (intersects(i)) return len + (start - i.getStart());
				latest = i.getEnd();
			}
			len += i.getEnd() - i.getStart();
		}

		return len + (start - latest);
	}

	public String getChange() {
		return change;
	}

	public String getChangeOption(int i) {
		return changeOptions[i];
	}

	public int getChangeOptionCount() {
		return changeOptions.length;
	}

	public ChangeType getChangeType() {
		return changeType;
	}

	public int getCoverage() {
		return coverage;
	}

	public double getQuality() {
		return quality;
	}

	public String getReference() {
		return reference;
	}

	public double getScore() {
		return score;
	}

	/**
	 * Create a new SeqChange for this option
	 * @param i
	 * @return
	 */
	public SeqChange getSeqChangeOption(int i) {
		// Just an interval (i.e. no changes)? => return the 'this' object
		if (changeType == ChangeType.Interval) return this;

		// Not a real change? return null
		// This might happen in heterozygous variants
		// E.g.: "A -> T/A", or "* -> +A/*")
		if (reference.equalsIgnoreCase(changeOptions[i])) return null;

		// Create new change
		return new SeqChange((Marker) parent, start, reference, changeOptions[i], strand, id, quality, coverage);
	}

	@Override
	public int getStrand() {
		return strand;
	}

	@Override
	public int hashCode() {
		int hashCode = getChromosomeName().hashCode();
		hashCode = hashCode * 31 + start;
		hashCode = hashCode * 31 + end;
		hashCode = hashCode * 31 + strand;
		hashCode = hashCode * 31 + id.hashCode();
		hashCode = hashCode * 31 + reference.hashCode();
		hashCode = hashCode * 31 + change.hashCode();
		return hashCode;
	}

	void init(Marker parent, int position, String referenceStr, String changeStr, int strand, String id, double quality, int coverage) {
		reference = referenceStr.toUpperCase();
		change = changeStr.toUpperCase();
		this.quality = quality;
		this.coverage = coverage;
		heterozygous = false;

		// Change type
		changeType = ChangeType.MNP;
		if (change.length() == 1) { // Is it a SNP?
			changeType = ChangeType.SNP;
			if (change.equals("A") || change.equals("C") || change.equals("G") || change.equals("T")) {
				heterozygous = false;
				changeOptions = new String[1];
				changeOptions[0] = change;
			} else {
				heterozygous = true;

				// Reference http://sourceforge.net/apps/mediawiki/samtools/index.php?title=SAM_FAQ#I_do_not_understand_the_columns_in_the_pileup_output.
				// IUB codes: M=A/C, R=A/G, W=A/T, S=C/G, Y=C/T, K=G/T and N=A/C/G/T
				if (change.length() == 1) {
					if (change.equals("N")) { // aNy base
						changeOptions = new String[4];
						changeOptions[0] = "A";
						changeOptions[1] = "C";
						changeOptions[2] = "G";
						changeOptions[3] = "T";
					} else if (change.equals("B")) { // B: not A
						changeOptions = new String[3];
						changeOptions[0] = "C";
						changeOptions[1] = "G";
						changeOptions[2] = "T";
					} else if (change.equals("D")) { // D: not C
						changeOptions = new String[3];
						changeOptions[0] = "A";
						changeOptions[1] = "G";
						changeOptions[2] = "T";
					} else if (change.equals("H")) { // H: not G
						changeOptions = new String[3];
						changeOptions[0] = "A";
						changeOptions[1] = "C";
						changeOptions[2] = "T";
					} else if (change.equals("V")) { // V: not T
						changeOptions = new String[3];
						changeOptions[0] = "A";
						changeOptions[1] = "C";
						changeOptions[2] = "G";
					} else if (change.equals("M")) {
						changeOptions = new String[2];
						changeOptions[0] = "A";
						changeOptions[1] = "C";
					} else if (change.equals("R")) {
						changeOptions = new String[2];
						changeOptions[0] = "A";
						changeOptions[1] = "G";
					} else if (change.equals("W")) { // Weak
						changeOptions = new String[2];
						changeOptions[0] = "A";
						changeOptions[1] = "T";
					} else if (change.equals("S")) { // Strong
						changeOptions = new String[2];
						changeOptions[0] = "C";
						changeOptions[1] = "G";
					} else if (change.equals("Y")) {
						changeOptions = new String[2];
						changeOptions[0] = "C";
						changeOptions[1] = "T";
					} else if (change.equals("K")) {
						changeOptions = new String[2];
						changeOptions[0] = "G";
						changeOptions[1] = "T";
					} else {
						throw new RuntimeException("WARNING: Unkown IUB code for SNP '" + change + "'");
					}
				}
			}
		} else {
			// Homozygous or heterozygous? Hetero if two parts (before and after the slash) are different
			if (change.indexOf(',') >= 0) changeOptions = change.split(",");
			else changeOptions = change.split("/");

			// Heterozygous if any of the changes is different
			if (changeOptions.length > 1) {
				heterozygous = false;
				for (int i = 0; i < (changeOptions.length - 1); i++) {
					heterozygous |= !changeOptions[i].equals(changeOptions[i + 1]);
				}
			}

			// Homozygous with more than one change? => Use only one change (since they are all the same, no need to analyze many times)
			if (!heterozygous && (changeOptions.length > 1)) {
				String changeOptionsNew[] = new String[1];
				changeOptionsNew[0] = changeOptions[0];
				changeOptions = changeOptionsNew;
			}

			// Use the asterisk at the end (e.g. when changeStr is '*/' instead of '/*')
			if (changeOptions[0].equals("*") && (changeOptions.length > 1)) {
				changeOptions[0] = changeOptions[1];
				changeOptions[1] = "*";
				change = changeOptions[0] + "/" + changeOptions[1];
			}

			// Insertions
			if (change.startsWith("+")) {
				changeType = ChangeType.INS;
			} else if (change.startsWith("-")) { // Deletions
				changeType = ChangeType.DEL;
				if (changeOptions[0].length() > 1) end = position + changeOptions[0].length() - 2; // Update 'end' position
			}

			// Insertions and deletions always have '*' as reference
			if ((changeType == ChangeType.INS) || (changeType == ChangeType.DEL)) reference = "*";
		}

		type = EffectType.NONE;

		// Start and end position
		// 	- Start is always the leftmost base
		//	- End is always the rightmost affected base in the reference genome
		start = position;
		if (isIns() || isSnp()) {
			// These changes only affect one position in the reference genome
			end = start;
		} else {
			for (int i = 0; i < changeOptions.length; i++) {
				String ch = changeOptions[i];
				if (ch.startsWith("+") || ch.startsWith("-")) ch = ch.substring(1);
				end = Math.max(end, start + ch.length() - 1);
			}
		}
	}

	/**
	 * Is this a change or are the changes actually the same as the reference
	 * @return
	 */
	public boolean isChange() {
		for (String chg : changeOptions)
			if (!reference.equals(chg)) return true; // Any change option is different? => true
		return false;
	}

	public boolean isDel() {
		return (changeType == ChangeType.DEL);
	}

	public boolean isHeterozygous() {
		return (heterozygous != null) && heterozygous;
	}

	public boolean isHomozygous() {
		return (heterozygous != null) && !heterozygous;
	}

	public boolean isImprecise() {
		return imprecise;
	}

	public boolean isInDel() {
		return (changeType == ChangeType.INS) || (changeType == ChangeType.DEL);
	}

	public boolean isIns() {
		return (changeType == ChangeType.INS);
	}

	public boolean isInterval() {
		return (changeType == ChangeType.Interval);
	}

	public boolean isMnp() {
		return changeType == ChangeType.MNP;
	}

	@Override
	protected boolean isShowWarningIfParentDoesNotInclude() {
		return false;
	}

	public boolean isSnp() {
		return changeType == ChangeType.SNP;
	}

	/**
	 * Return the change (always compared to 'referenceStrand') without any '+' or '-' leading characters
	 * @return
	 */
	public String netChange(int referenceStrand) {
		String netChange = change;
		if (change.startsWith("+") || change.startsWith("-")) netChange = change.substring(1); // Remove leading char

		// Need reverse-WC?
		boolean needRwc = false;
		if ((strand >= 0) && (referenceStrand < 0)) needRwc = true;
		else if ((strand < 0) && (referenceStrand >= 0)) needRwc = true;

		return needRwc ? GprSeq.reverseWc(netChange) : netChange;
	}

	/**
	 * Only the part of the change that overlaps with a marker
	 * Return the change (always in positive strand) without any '+' or '-' leading characters
	 * @return
	 */
	public String netChange(Marker marker) {
		String netChange = change;
		if (change.startsWith("+") || change.startsWith("-")) netChange = change.substring(1); // Remove leading char

		int removeBefore = marker.getStart() - start;
		if (removeBefore > 0) {
			if (removeBefore >= netChange.length()) return ""; // Nothing left
		} else removeBefore = 0;

		int removeAfter = end - marker.getEnd();
		if (removeAfter > 0) {
			if ((removeBefore + removeAfter) >= netChange.length()) return ""; // Nothing left
		} else removeAfter = 0;

		// Use reverse-WC?
		if (strand <= 0) netChange = GprSeq.reverseWc(netChange);

		// Remove leading and trailing parts
		netChange = netChange.substring(removeBefore, netChange.length() - removeAfter);

		return netChange;
	}

	/**
	 * Return the reference (always in positive strand)
	 * @return
	 */
	public String reference() {
		return strand >= 0 ? reference : GprSeq.reverseWc(reference);
	}

	public void setChangeType(ChangeType changeType) {
		this.changeType = changeType;
	}

	public void setHeterozygous(Boolean heterozygous) {
		this.heterozygous = heterozygous;
	}

	public void setImprecise(boolean imprecise) {
		this.imprecise = imprecise;
	}

	public void setScore(double score) {
		this.score = score;
	}

	@Override
	public String toString() {
		if (changeType == ChangeType.Interval) return "chr" + getChromosomeName() + ":" + start + "-" + end;

		return "chr" + getChromosomeName() //
				+ ":" + start //
				+ "_" + getReference() //
				+ "/" + getChange() //
				+ ((id != null) && (id.length() > 0) ? " '" + id + "'" : "");
	}

	/**
	 * Show it required by ENSEMBL's SnpEffectPredictor
	 * @return
	 */
	public String toStringEnsembl() {
		return getChromosomeName() + "\t" + start + "\t" + end + "\t" + reference + "/" + change + "\t+";
	}

}
