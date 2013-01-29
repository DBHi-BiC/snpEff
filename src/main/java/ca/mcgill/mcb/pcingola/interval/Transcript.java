package ca.mcgill.mcb.pcingola.interval;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import ca.mcgill.mcb.pcingola.interval.codonChange.CodonChange;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.stats.ObservedOverExpectedCpG;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * Codon position
 * @author pcingola
 */
class CodonPosition {

	public int codonNum = -1;
	public int codonIndex = -1;
}

/**
 * Interval for a transcript, as well as some other information: exons, utrs, cds, etc.
 * 
 * @author pcingola
 */
public class Transcript extends IntervalAndSubIntervals<Exon> {

	private static final long serialVersionUID = -2665025617916107311L;

	boolean proteinCoding = false; // Is this a protein-coding transcript?
	int cdsStart = -1;
	int cdsEnd = -1;
	String bioType = ""; // Transcript biotype
	String cds = null; // Coding sequence
	ArrayList<SpliceSiteBranch> spliceBranchSites; // Branch splice sites
	ArrayList<Utr> utrs; // UTRs
	ArrayList<Cds> cdss; // CDS information
	ArrayList<Intron> introns; // Intron markers
	Upstream upstream; // Upstream interval
	Downstream downstream; // Downstream interval
	Exon firstCodingExon; // First coding exon. I.e. where transcription start site (TSS) is. 
	int cds2pos[];

	protected Transcript() {
		super();
		spliceBranchSites = new ArrayList<SpliceSiteBranch>();
		utrs = new ArrayList<Utr>();
		cdss = new ArrayList<Cds>();
		type = EffectType.TRANSCRIPT;
	}

	public Transcript(Gene gene, int start, int end, int strand, String id) {
		super(gene, start, end, strand, id);
		spliceBranchSites = new ArrayList<SpliceSiteBranch>();
		utrs = new ArrayList<Utr>();
		cdss = new ArrayList<Cds>();
		type = EffectType.TRANSCRIPT;
	}

	/**
	 * Add a CDS
	 * @param cdsInt
	 */
	public void add(Cds cdsInt) {
		cdss.add(cdsInt);
		cds = null;
	}

	/**
	 * Add a SpliceSiteBranchU12
	 * @param branchU12
	 */
	public void add(SpliceSiteBranchU12 branchU12) {
		spliceBranchSites.add(branchU12);
	}

	/**
	 * Add a UTR
	 * @param utr
	 */
	public void add(Utr utr) {
		utrs.add(utr);
		cds = null;
	}

	/**
	 * Add missing UTRs. See utrFromCds() method.
	 * @param missingUtrs
	 */
	boolean addMissingUtrs(Markers missingUtrs, boolean verbose) {
		missingUtrs.sort(false, strand < 0);

		// Get min/max CDS positions
		int minCds = Integer.MAX_VALUE;
		int maxCds = 0;
		for (Cds c : cdss) {
			minCds = Math.min(minCds, c.getStart());
			maxCds = Math.max(maxCds, c.getEnd());
		}

		if (verbose) System.out.println("Transcript '" + id + "' has missing UTRs. Strand: " + strand + " (minCds: " + minCds + " , maxCds: " + maxCds + "):");

		// Add intervals
		boolean retVal = false;
		for (Marker mu : missingUtrs) {
			Exon exon = intersectingExon(mu);
			if (exon == null) throw new RuntimeException("Cannot find exon for UTR: " + mu);
			Utr toAdd = null;

			if (isStrandPlus()) {
				if (mu.getEnd() <= minCds) toAdd = new Utr5prime(exon, mu.getStart(), mu.getEnd(), strand, mu.getId());
				else if (mu.getStart() >= maxCds) toAdd = new Utr3prime(exon, mu.getStart(), mu.getEnd(), strand, mu.getId());
			} else {
				if (mu.getStart() >= maxCds) toAdd = new Utr5prime(exon, mu.getStart(), mu.getEnd(), strand, mu.getId());
				else if (mu.getEnd() <= minCds) toAdd = new Utr3prime(exon, mu.getStart(), mu.getEnd(), strand, mu.getId());
			}

			// OK?
			if (toAdd != null) {
				add(toAdd);
				if (verbose) System.out.println("\tAdding " + toAdd);
				retVal = true;
			}
		}

		return retVal;
	}

	/**
	 * Adjust transcript coordinates
	 * @return
	 */
	public boolean adjust() {
		boolean changed = false;
		int strandSumTr = 0;
		int newStart = start, newEnd = end;
		if (newStart == 0 && newEnd == 0) {
			newStart = Integer.MAX_VALUE;
			newEnd = Integer.MIN_VALUE;
		}

		int countStrandPlus = 0, countStrandMinus = 0;
		for (Exon exon : sortedStrand()) {
			newStart = Math.min(newStart, exon.getStart());
			newEnd = Math.max(newEnd, exon.getEnd());

			// Conun exon strands
			if (exon.getStrand() > 0) countStrandPlus++;
			else if (exon.getStrand() < 0) countStrandMinus++;
		}

		// UTRs
		for (Utr utr : getUtrs()) {
			newStart = Math.min(newStart, utr.getStart());
			newEnd = Math.max(newEnd, utr.getEnd());
		}

		// Sanity check
		strandSumTr = countStrandPlus - countStrandMinus; // Some exons have incorrect strands, we use the strand indicated by most exons
		int newStrand = strandSumTr >= 0 ? 1 : -1;
		if ((countStrandPlus > 0) && (countStrandMinus > 0)) Gpr.debug("Transcript '" + id + "' has " + countStrandPlus + " exons on the plus and " + countStrandMinus + " exons on the minus strand! This should never happen!");

		// Change transcript strand?
		if (strand != newStrand) {
			changed = true;
			setStart(newStrand); // Change strand
		}

		// Change start?
		if (start != newStart) {
			start = newStart;
			changed = true;
		}

		// Change end?
		if (end != newEnd) {
			end = newEnd;
			changed = true;
		}

		return changed;
	}

	/**
	 * Calculate CDS start and CDS end
	 */
	synchronized void calcCdsStartEnd() {
		// Do we need to calculate these values?
		if (cdsStart < 0) {
			// Calculate coding start (after 5 prime UTR)

			if (utrs.isEmpty()) {
				// No UTRs => Use all exons
				cdsStart = (isStrandPlus() ? end : start); // cdsStart is the position of the first base in the CDS (i.e. the first base after all 5'UTR)
				cdsEnd = (isStrandPlus() ? start : end); // cdsEnd is the position of the last base in the CDS (i.e. the first base before all 3'UTR)

				for (Exon ex : this) {
					if (isStrandPlus()) {
						cdsStart = Math.min(cdsStart, ex.getStart());
						cdsEnd = Math.max(cdsEnd, ex.getEnd());
					} else {
						cdsStart = Math.max(cdsStart, ex.getEnd());
						cdsEnd = Math.min(cdsEnd, ex.getStart());
					}
				}
			} else {
				// We have to take into account UTRs
				cdsStart = (isStrandPlus() ? start : end); // cdsStart is the position of the first base in the CDS (i.e. the first base after all 5'UTR)
				cdsEnd = (isStrandPlus() ? end : start); // cdsEnd is the position of the last base in the CDS (i.e. the first base before all 3'UTR)
				int cdsStartNotExon = cdsStart;

				for (Utr utr : utrs) {
					if (utr instanceof Utr5prime) {
						if (isStrandPlus()) cdsStart = Math.max(cdsStart, utr.getEnd() + 1);
						else cdsStart = Math.min(cdsStart, utr.getStart() - 1);
					} else if (utr instanceof Utr3prime) {
						if (isStrandPlus()) cdsEnd = Math.min(cdsEnd, utr.getStart() - 1);
						else cdsEnd = Math.max(cdsEnd, utr.getEnd() + 1);
					}
				}

				// Make sure cdsStart and cdsEnd lie within an exon
				if (isStrandPlus()) {
					cdsStart = firstExonPositionAfter(cdsStart);
					cdsEnd = lastExonPositionBefore(cdsEnd);
				} else {
					cdsStart = lastExonPositionBefore(cdsStart);
					cdsEnd = firstExonPositionAfter(cdsEnd);
				}

				// We were not able to find cdsStart & cdsEnd within exon limits. 
				// Probably there is something wrong with the database and the transcript does 
				// not have a single coding base (e.g. all of it is UTR).
				if (cdsStart < 0 || cdsEnd < 0) cdsStart = cdsEnd = cdsStartNotExon;
			}
		}
	}

	/**
	 * Retrieve coding sequence
	 */
	public synchronized String cds() {
		if (cds != null) return cds;

		// Concatenate all exons
		List<Exon> exons = sortedStrand();
		StringBuilder sequence = new StringBuilder();
		int utr5len = 0, utr3len = 0;

		// 5 prime UTR length
		for (Utr utr : get5primeUtrs())
			utr5len += utr.size();

		// Append all exon sequences
		boolean missingSequence = false;
		for (Exon exon : exons) {
			missingSequence |= !exon.hasSequence(); // If there is no sequence, we are in trouble
			sequence.append(exon.getSequence());
		}

		if (missingSequence) cds = ""; // One or more exons does not have sequence. Nothing to do
		else {
			// OK, all exons have sequences

			// 3 prime UTR length
			for (Utr utr : get3primeUtrs())
				utr3len += utr.size();

			// Cut 5 prime UTR and 3 prime UTR points
			int subEnd = sequence.length() - utr3len;

			if (utr5len > subEnd) cds = "";
			else cds = sequence.substring(utr5len, subEnd);
		}

		return cds;
	}

	/**
	 * Calculate base number in a CDS where 'pos' maps
	 * 
	 * @returns Base number or '-1' if it does not map to a coding base
	 */
	public synchronized int cdsBaseNumber(int pos, boolean usePrevBaseIntron) {
		// Doesn't hit this transcript?
		if (!intersects(pos)) return -1;

		// Is it in UTR instead of CDS? 
		if (isUtr(pos)) return -1;

		// Calculate cdsStart and cdsEnd (if not already done)
		calcCdsStartEnd();

		// All exons..
		int firstCdsBaseInExon = 0; // Where the exon maps to the CDS (i.e. which CDS base number does the first base in this exon maps to).
		for (Exon eint : sortedStrand()) {
			if (eint.intersects(pos)) {
				int cdsBaseInExon; // cdsBaseInExon: base number relative to the beginning of the coding part of this exon (i.e. excluding 5'UTRs)
				if (strand >= 0) cdsBaseInExon = pos - Math.max(eint.getStart(), cdsStart);
				else cdsBaseInExon = Math.min(eint.getEnd(), cdsStart) - pos;

				cdsBaseInExon = Math.max(0, cdsBaseInExon);

				return firstCdsBaseInExon + cdsBaseInExon;
			} else {
				// Before exon begins?
				if ((isStrandPlus() && (pos < eint.getStart())) // Before exon begins (positive strand)?
						|| (isStrandMinus() && (pos > eint.getEnd()))) // Before exon begins (negative strand)?
					return firstCdsBaseInExon - (usePrevBaseIntron ? 1 : 0);
			}

			if (isStrandPlus()) firstCdsBaseInExon += Math.max(0, eint.getEnd() - Math.max(eint.getStart(), cdsStart) + 1);
			else firstCdsBaseInExon += Math.max(0, Math.min(cdsStart, eint.getEnd()) - eint.getStart() + 1);
		}

		return firstCdsBaseInExon - 1;
	}

	/**
	 * Calculate chromosome position as function of CDS number 
	 * 
	 * @returns An array mapping 'pos[cdsBaseNumber] = chromosmalPos' 
	 */
	public synchronized int[] cdsBaseNumber2ChrPos() {
		if (cds2pos != null) return cds2pos;

		calcCdsStartEnd();

		cds2pos = new int[cds().length()];
		for (int i = 0; i < cds2pos.length; i++)
			cds2pos[i] = -1;

		int cdsMin = Math.min(cdsStart, cdsEnd);
		int cdsMax = Math.max(cdsStart, cdsEnd);

		// For each exon, add CDS position to array
		int cdsBaseNum = 0;
		for (Exon exon : sortedStrand()) {
			int min = isStrandPlus() ? exon.getStart() : exon.getEnd();
			int step = isStrandPlus() ? 1 : -1;
			for (int pos = min; exon.intersects(pos); pos += step)
				if ((cdsMin <= pos) && (pos <= cdsMax)) cds2pos[cdsBaseNum++] = pos;
		}

		return cds2pos;
	}

	/**
	 * Analyze SNPs in this transcript.
	 * Add changeEffect to 'changeEffect'
	 */
	public String codonByCdsBaseNumber(int cdsBaseNumber) {
		int codonNum = cdsBaseNumber / CodonChange.CODON_SIZE;
		int min = codonNum * CodonChange.CODON_SIZE;
		int max = codonNum * CodonChange.CODON_SIZE + CodonChange.CODON_SIZE;
		if ((min >= 0) && (max <= cds().length())) return cds().substring(min, max).toUpperCase();
		return null;
	}

	/**
	 * Collapses exons having gaps of zero (i.e. exons that followed by other exons).
	 * Does the same for CDSs.
	   Does the same for UTRs.
	 */
	public boolean collapseZeroGap() {
		boolean ret = false;
		introns = null; // These need to be recalculated

		//---
		// Collapse Exons
		//---
		Map<Marker, Marker> collapse = MarkerUtil.collapseZeroGap(new Markers().addAll(subintervals())); // Create a map of collapsed exons

		// Replace exons
		for (Marker exon : collapse.keySet()) {
			Exon collapsedExon = (Exon) collapse.get(exon);

			// Is this exon to be replaced?
			if (exon != collapsedExon) {
				ret = true;

				// Replace exon
				remove((Exon) exon);
				if (!containsId(collapsedExon.getId())) add(collapsedExon); // Add collapsedExon. Make sure we don't add it twice (since many exons can be collapsed into one).

				// Change parent exon in UTRs
				for (Marker m : getUtrs()) {
					Utr utr = (Utr) m;
					if (utr.getParent() == exon) utr.setParent(collapsedExon);
				}
			}
		}

		//---
		// Collapse CDS
		//---
		collapse = MarkerUtil.collapseZeroGap(new Markers(cdss));
		cdss = new ArrayList<Cds>(); // Re-create CDSs list
		Markers uniqCollapsedCds = new Markers(collapse.values()).unique(); // Create a set of unique CDSs and add them to CDSs list
		for (Marker cds : uniqCollapsedCds)
			cdss.add((Cds) cds);

		//---
		// Collapse UTRs
		//---
		collapse = MarkerUtil.collapseZeroGap(new Markers(utrs));
		Markers uniqCollapsedUtrs = new Markers(collapse.values()).unique(); // Create a set of unique UTRs, and add them to the list
		utrs = new ArrayList<Utr>(); // Re-generate UTRs list
		for (Marker utr : uniqCollapsedUtrs)
			utrs.add((Utr) utr);

		return ret;
	}

	/**
	 * Calculate CpG bias: number of CpG / expected[CpG]
	 * @return
	 */
	public double cpgExonBias() {
		ObservedOverExpectedCpG oe = new ObservedOverExpectedCpG();
		return oe.oe(this);
	}

	/**
	 * Count total CpG in this transcript's exons
	 * @return
	 */
	public int cpgExons() {
		ObservedOverExpectedCpG oe = new ObservedOverExpectedCpG();
		return oe.observed(this);
	}

	/**
	 * Creates a list of UP/DOWN stream regions (for each transcript)
	 * Upstream (downstream) stream is defined as upDownLength before (after) transcript
	 */
	public void createUpDownStream(int upDownLength) {
		Chromosome chr = getChromosome();
		int min = chr.getStart(), max = chr.getEnd();

		// Create up/down stream intervals and add them to the list
		if (isStrandPlus()) {
			upstream = new Upstream(this, Math.max(start - upDownLength, min), Math.max(start - 1, min), 1, id);
			downstream = new Downstream(this, Math.min(end + 1, max), Math.min(end + upDownLength, max), 1, id);
		} else {
			upstream = new Upstream(this, Math.min(end + 1, max), Math.min(end + upDownLength, max), 1, id);
			downstream = new Downstream(this, Math.max(start - upDownLength, min), Math.max(start - 1, min), 1, id);
		}
	}

	/**
	 * Deletes redundant exons (i.e. exons that are totally included in other exons).
	 * Does the same for CDSs.
	   Does the same for UTRs.
	 */
	public boolean deleteRedundant() {
		boolean ret = false;
		introns = null; // These need to be recalculated

		//---
		// Delete redundant exons
		//---
		Map<Marker, Marker> includedIn = MarkerUtil.redundant(subintervals());
		for (Marker exon : includedIn.keySet()) {
			ret = true;
			remove((Exon) exon);

			// Change parent exon in UTRs
			for (Marker m : getUtrs()) {
				Utr utr = (Utr) m;
				if (utr.getParent() == exon) utr.setParent(includedIn.get(exon));
			}
		}

		//---
		// Delete redundant CDS
		//---
		includedIn = MarkerUtil.redundant(cdss);
		for (Marker cds : includedIn.keySet())
			cdss.remove(cds);

		//---
		// Delete redundant UTRs
		//---
		includedIn = MarkerUtil.redundant(utrs);
		for (Marker utr : includedIn.keySet())
			utrs.remove(utr);

		return ret;
	}

	/**
	 * Return the Exon that hits position 'pos'
	 * @param pos
	 * @return An exon intersecting 'pos' (null if not found)
	 */
	public Exon findExon(int pos) {
		// Is it in UTR instead of CDS? 
		for (Exon exon : this)
			if (exon.intersects(pos)) return exon;
		return null;
	}

	/**
	 * Find all splice sites.
	 * 
	 * @param createIfMissing : If true, create canonical splice sites if they are missing.
	 * 
	 * @return
	 */
	public List<SpliceSite> findSpliceSites(boolean createIfMissing) {
		List<SpliceSite> list = new LinkedList<SpliceSite>();

		// For each gene, transcript and exon
		ArrayList<Exon> exons = (ArrayList<Exon>) sortedStrand();

		if (exons.size() > 0) {
			for (int i = 0; i < exons.size(); i++) {

				Exon exon = exons.get(i);
				Exon prev = (i >= 1 ? exons.get(i - 1) : null);
				Exon next = (i < exons.size() - 1 ? exons.get(i + 1) : null);

				//---
				// Distance to previous exon
				//---
				if (prev != null) {
					int dist = 0;
					if (strand >= 0) dist = exon.getStart() - prev.getEnd() - 1;
					else dist = prev.getStart() - exon.getEnd() - 1;

					// Acceptor splice site: before exon start, but not before first exon
					SpliceSite ss = exon.getSpliceSiteAcceptor();
					if ((ss == null) && createIfMissing) ss = exon.createSpliceSiteAcceptor(Math.min(SpliceSite.CORE_SPLICE_SITE_SIZE, dist));
					if (ss != null) list.add(ss);
				}

				//---
				// Distance to next exon
				//---
				if (next != null) {
					int dist = 0;
					if (strand >= 0) dist = next.getStart() - exon.getEnd() - 1;
					else dist = exon.getStart() - next.getEnd() - 1;

					// Donor splice site: after exon end, but not after last exon
					SpliceSite ss = exon.getSpliceSiteDonor();
					if ((ss == null) && createIfMissing) ss = exon.createSpliceSiteDonor(Math.min(SpliceSite.CORE_SPLICE_SITE_SIZE, dist));
					if (ss != null) list.add(ss);
				}

				// Sanity check
				int rank = i + 1;
				if (exon.getRank() != rank) {
					String msg = "Rank numbers do not march: " + rank + " != " + exon.getRank() + "\n\tTranscript: " + this;
					throw new RuntimeException(msg);
				}
			}
		}

		return list;
	}

	/**
	 * Return the UTR that hits position 'pos'
	 * @param pos
	 * @return An UTR intersecting 'pos' (null if not found)
	 */
	public Utr findUtr(int pos) {
		// Is it in UTR instead of CDS? 
		for (Utr utr : utrs)
			if (utr.intersects(pos)) return utr;
		return null;
	}

	/**
	 * Find the first position after 'pos' within an exon
	 * @param pos
	 * @return
	 */
	int firstExonPositionAfter(int pos) {
		for (Exon ex : sorted()) {
			if (pos <= ex.getStart()) return ex.getStart();
			if (pos <= ex.getEnd()) return pos;
		}

		System.err.println("WARNING: Cannot find first exonic position after " + pos + " for transcript '" + id + "'");
		return -1;
	}

	/**
	 * Create a list of 3 prime UTRs
	 */
	public List<Utr3prime> get3primeUtrs() {
		ArrayList<Utr3prime> list = new ArrayList<Utr3prime>();
		for (Utr utr : utrs)
			if (utr instanceof Utr3prime) list.add((Utr3prime) utr);
		return list;
	}

	/**
	 * Create a list of 5 prime UTRs
	 */
	public List<Utr5prime> get5primeUtrs() {
		ArrayList<Utr5prime> list = new ArrayList<Utr5prime>();
		for (Utr utr : utrs)
			if (utr instanceof Utr5prime) list.add((Utr5prime) utr);
		return list;
	}

	public String getBioType() {
		return bioType;
	}

	/**
	 * Get all CDSs
	 * @return
	 */
	public List<Cds> getCds() {
		return cdss;
	}

	public int getCdsEnd() {
		calcCdsStartEnd();
		return cdsEnd;
	}

	public int getCdsStart() {
		calcCdsStartEnd();
		return cdsStart;
	}

	public Downstream getDownstream() {
		return downstream;
	}

	/**
	 * Get first coding exon
	 * @return
	 */
	public synchronized Exon getFirstCodingExon() {
		if (firstCodingExon == null) {
			// Get transcription start position
			long cstart = getCdsStart();

			// Pick exon intersecting cdsStart (TSS)
			for (Exon exon : sortedStrand())
				if (exon.intersects(cstart)) firstCodingExon = exon;

			// Sanity check
			if (firstCodingExon == null) throw new RuntimeException("Error: Cannot find first coding exon for transcript:\n" + this);
		}
		return firstCodingExon;
	}

	public ArrayList<SpliceSiteBranch> getSpliceBranchSites() {
		return spliceBranchSites;
	}

	public Upstream getUpstream() {
		return upstream;
	}

	/**
	 * Get all UTRs
	 * @return
	 */
	public List<Utr> getUtrs() {
		return utrs;
	}

	/**
	 * Return the first exon that intersects 'interval' (null if not found)
	 * @param interval
	 * @return
	 */
	public Exon intersectingExon(Marker interval) {
		for (Exon ei : this)
			if (ei.intersects(interval)) return ei;
		return null;
	}

	/**
	 * Get all introns (lazy init)
	 * @return
	 */
	public synchronized ArrayList<Intron> introns() {
		if (introns == null) {
			introns = new ArrayList<Intron>();

			Exon exBefore = null;
			for (Exon ex : sortedStrand()) {
				if (exBefore != null) {
					// Create intron
					Intron intron;
					int rank = introns.size() + 1;

					// Find intron start and end
					int start, end;
					if (isStrandPlus()) {
						start = exBefore.getEnd() + 1;
						end = ex.getStart() - 1;
					} else {
						start = ex.getEnd() + 1;
						end = exBefore.getStart() - 1;
					}

					int size = end - start + 1;
					if (size > 0) {
						// Add intron to list
						intron = new Intron(this, start, end, strand, id + "_intron_" + rank);
						intron.setRank(rank);
						introns.add(intron);
					}
				}
				exBefore = ex;
			}
		}
		return introns;
	}

	/**
	 * Return the intron size for intron number 'intronNum'
	 * 
	 * Note: Intron number 'N' is the intron between exon number N and exon number N+1
	 * Note: Numbering is zero-based (not to be confused with exon 'ranking', which is one-based)
	 * 
	 * @param intronNum
	 * @return
	 */
	public int intronSize(int intronNum) {
		if (intronNum >= (numChilds() - 1)) return -1;
		ArrayList<Exon> exons = (ArrayList<Exon>) sortedStrand();
		Exon exon = exons.get(intronNum);
		Exon next = exons.get(intronNum + 1);
		return (isStrandPlus() ? (next.getStart() - exon.getEnd()) : (exon.getStart() - next.getEnd())) - 1;
	}

	@Override
	protected boolean isAdjustIfParentDoesNotInclude(Marker parent) {
		return true;
	}

	/**
	 * Is this seqChange in the CDS part of this transcript?
	 * @param seqChange
	 * @return
	 */
	boolean isCds(SeqChange seqChange) {
		calcCdsStartEnd();

		int cs = cdsStart;
		int ce = cdsEnd;

		if (isStrandMinus()) {
			cs = cdsEnd;
			ce = cdsStart;
		}

		return (seqChange.getEnd() >= cs) && (seqChange.getStart() <= ce);
	}

	public boolean isProteinCoding() {
		return proteinCoding;
	}

	/**
	 * Does this 'pos' hit a UTR?
	 * @param pos
	 * @return
	 */
	public boolean isUtr(int pos) {
		return findUtr(pos) != null;
	}

	/**
	 * Find the last position before 'pos' within an exon
	 * @param pos
	 * @return
	 */
	int lastExonPositionBefore(int pos) {
		int last = -1;
		for (Exon ex : sorted()) {
			if (pos < ex.getStart()) {
				// Nothing found?
				if (last < 0) {
					System.err.println("WARNING: Cannot find last exonic position before " + pos + " for transcript '" + id + "'");
					return -1;
				}
				return last;
			} else if (pos <= ex.getEnd()) return pos;
			last = ex.getEnd();
		}

		System.err.println("WARNING: Cannot find last exonic position before " + pos + " for transcript '" + id + "'");
		return -1;
	}

	/**
	 * Retrieve coding sequence AND the UTRs (mRNA = 5'UTR + CDS + 3'UTR)
	 * I.e. Concatenate all exon sequences
	 */
	public String mRna() {
		List<Exon> exons = sortedStrand();

		// Concatenate all exons
		StringBuilder sequence = new StringBuilder();
		for (Exon eint : exons)
			sequence.append(eint.getSequence());

		return sequence.toString();
	}

	/**
	 * Protein sequence (amino acid sequence produced by this transcripts) 
	 * @return
	 */
	public String protein() {
		if (!Config.get().isTreatAllAsProteinCoding() && !isProteinCoding()) return "";
		return codonTable().aa(cds());
	}

	/**
	 * Assign ranks to exons
	 */
	public boolean rankExons() {
		boolean changed = false;
		int rank = 1;
		for (Exon exon : sortedStrand()) {
			if (rank != exon.getRank()) {
				exon.setRank(rank);
				changed = true;
			}
			rank++;
		}
		return changed;
	}

	/**
	 * Get some details about the effect on this transcript
	 * @param seqChange
	 * @return
	 */
	@Override
	public List<ChangeEffect> seqChangeEffect(SeqChange seqChange, ChangeEffect changeEffect) {
		if (!intersects(seqChange)) return ChangeEffect.emptyResults(); // Sanity check

		// Create a list of changes
		ArrayList<ChangeEffect> changeEffectList = new ArrayList<ChangeEffect>();

		//---
		// Hits a UTR region?
		//---
		boolean included = false;
		for (Utr utr : utrs)
			if (utr.intersects(seqChange)) {
				// Calculate the effect
				List<ChangeEffect> chEffList = utr.seqChangeEffect(seqChange, changeEffect.clone());
				if (!chEffList.isEmpty()) changeEffectList.addAll(chEffList);
				included |= utr.includes(seqChange); // Is this seqChange fully included in the UTR?
			}
		if (included) return changeEffectList; // SeqChange fully included in the UTR? => We are done.

		//---
		// Hits a SpliceSiteBranch region?
		//---
		included = false;
		for (SpliceSiteBranch ssbranch : spliceBranchSites)
			if (ssbranch.intersects(seqChange)) {
				// Calculate the effect
				List<ChangeEffect> chEffList = ssbranch.seqChangeEffect(seqChange, changeEffect.clone());
				if (!chEffList.isEmpty()) changeEffectList.addAll(chEffList);
				included |= ssbranch.includes(seqChange); // Is this seqChange fully included branch site?
			}
		if (included) return changeEffectList; // SeqChange fully included in the Branch site? => We are done.

		// Does it hit an intron?
		for (Intron intron : introns())
			if (intron.intersects(seqChange)) {
				ChangeEffect cheff = changeEffect.clone();
				cheff.set(intron, EffectType.INTRON, "");
				changeEffectList.add(cheff);
			}

		//---
		// Analyze non-coding transcripts (or 'interval' seqChanges)
		//---
		if ((!Config.get().isTreatAllAsProteinCoding() && !isProteinCoding()) || seqChange.isInterval()) {
			// Do we have exon information for this transcript?
			if (!subintervals().isEmpty()) {
				// Add all exons
				for (Exon exon : this) {
					if (exon.intersects(seqChange)) {
						ChangeEffect cheff = changeEffect.clone();
						cheff.set(exon, EffectType.EXON, "");
						changeEffectList.add(cheff);
					}
				}
			} else {
				// No exons annotated? Just mark it as hitting a transcript
				ChangeEffect cheff = changeEffect.clone();
				cheff.set(this, EffectType.TRANSCRIPT, "");
				changeEffectList.add(cheff);
			}

			return changeEffectList;
		}

		//---
		// This is a protein coding transcript.
		// We analyze codon replacement effect
		//---
		if (isCds(seqChange)) {
			// Get codon change effect 
			CodonChange codonChange = new CodonChange(seqChange, this, changeEffect);
			changeEffectList.addAll(codonChange.calculate());
		}

		return changeEffectList;
	}

	/**
	 * Parse a line from a serialized file
	 * @param line
	 * @return
	 */
	@Override
	public void serializeParse(MarkerSerializer markerSerializer) {
		super.serializeParse(markerSerializer);
		bioType = markerSerializer.getNextField();
		proteinCoding = markerSerializer.getNextFieldBoolean();
		upstream = (Upstream) markerSerializer.getNextFieldMarker();
		downstream = (Downstream) markerSerializer.getNextFieldMarker();

		for (Marker m : markerSerializer.getNextFieldMarkers())
			utrs.add((Utr) m);

		for (Marker m : markerSerializer.getNextFieldMarkers())
			cdss.add((Cds) m);

		for (Marker m : markerSerializer.getNextFieldMarkers())
			spliceBranchSites.add((SpliceSiteBranchU12) m);
	}

	/**
	 * Create a string to serialize to a file
	 * @return
	 */
	@SuppressWarnings({ "unchecked", "rawtypes" })
	@Override
	public String serializeSave(MarkerSerializer markerSerializer) {
		return super.serializeSave(markerSerializer) //
				+ "\t" + bioType //
				+ "\t" + proteinCoding //
				+ "\t" + markerSerializer.save(upstream) //
				+ "\t" + markerSerializer.save(downstream) //
				+ "\t" + markerSerializer.save((Iterable) utrs)//
				+ "\t" + markerSerializer.save((Iterable) cdss)//
				+ "\t" + markerSerializer.save((Iterable) spliceBranchSites)//
		;
	}

	public void setBioType(String bioType) {
		this.bioType = bioType;
	}

	public void setProteinCoding(boolean proteinCoding) {
		this.proteinCoding = proteinCoding;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		sb.append(getChromosomeName() + ":" + start + "-" + end);
		sb.append(", strand: " + (strand >= 0 ? "+" : "-"));
		if ((id != null) && (id.length() > 0)) sb.append(", id:" + id);
		if ((bioType != null) && (bioType.length() > 0)) sb.append(", bioType:" + bioType);
		if (isProteinCoding()) sb.append(", Protein");

		if (numChilds() > 0) {
			sb.append("\n");
			for (Utr utr : get5primeUtrs())
				sb.append("\t\t5'UTR   :\t" + utr + "\n");

			sb.append("\t\tExons:\n");
			for (Exon eint : sorted())
				sb.append("\t\t" + eint + "\n");

			for (Utr utr : get3primeUtrs())
				sb.append("\t\t3'UTR   :\t" + utr + "\n");

			// We may show CDS
			if (isProteinCoding()) {
				sb.append("\t\tCDS     :\t" + cds() + "\n");
				sb.append("\t\tProtein :\t" + protein() + "\n");
			}
		}

		return sb.toString();
	}

	/**
	 * Show a transcript as an ASCII Art
	 * @return
	 */
	public String toStringAsciiArt() {
		char art[] = new char[size()];
		for (int i = start, j = 0; i <= end; i++, j++) {
			Utr utr = findUtr(i);
			if (utr != null) art[j] = utr.isUtr5prime() ? '5' : '3';
			else {
				Exon exon = findExon(i);
				if (exon != null) art[j] = exon.isStrandPlus() ? '>' : '<';
				else art[j] = '-';
			}
		}
		return new String(art);
	}

	/**
	 * Calculate UTR regions from CDSs
	 */
	public boolean utrFromCds(boolean verbose) {
		if (cdss.size() <= 0) return false; // Cannot do this if we don't have CDS information

		// All exons minus all UTRs and CDS should give us the missing UTRs
		Markers exons = new Markers();
		Markers minus = new Markers();

		// Add all exons
		for (Exon e : this)
			exons.add(e);

		// Add all UTRs and CDSs to the 'minus' set
		for (Utr uint : getUtrs())
			minus.add(uint);

		for (Cds cint : cdss)
			minus.add(cint);

		Markers missingUtrs = exons.minus(minus); // Perform interval minus
		if (missingUtrs.size() > 0) return addMissingUtrs(missingUtrs, verbose); // Anything left? => There was a missing UTR
		return false;
	}
}
