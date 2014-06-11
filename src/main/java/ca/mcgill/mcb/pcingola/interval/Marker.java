package ca.mcgill.mcb.pcingola.interval;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import ca.mcgill.mcb.pcingola.codons.CodonTable;
import ca.mcgill.mcb.pcingola.codons.CodonTables;
import ca.mcgill.mcb.pcingola.serializer.MarkerSerializer;
import ca.mcgill.mcb.pcingola.serializer.TxtSerializable;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffects;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * An interval intended as a mark
 * 
 * @author pcingola
 */
public class Marker extends Interval implements TxtSerializable {

	private static final long serialVersionUID = 7878886900660027549L;
	protected EffectType type = EffectType.NONE;

	protected Marker() {
		super();
		type = EffectType.NONE;
	}

	public Marker(Marker parent, int start, int end, int strand, String id) {
		super(parent, start, end, strand, id);

		// Adjust parent if child is not included?
		if ((parent != null) && !parent.includes(this)) {
			String err = "";
			if (isShowWarningIfParentDoesNotInclude()) err = "WARNING: " + this.getClass().getSimpleName() + " is not included in parent " + parent.getClass().getSimpleName() + ". " //
					+ "\t" + this.getClass().getSimpleName() + " '" + getId() + "'  [ " + getStart() + " , " + getEnd() + " ]" //
					+ "\t" + parent.getClass().getSimpleName() + " '" + parent.getId() + "' [ " + parent.getStart() + " , " + parent.getEnd() + " ]";

			// Adjust parent?
			if (isAdjustIfParentDoesNotInclude(parent)) {
				parent.adjust(this);
				if (isShowWarningIfParentDoesNotInclude()) err += "\t=> Adjusting " + parent.getClass().getSimpleName() + " '" + parent.getId() + "' to [ " + parent.getStart() + " , " + parent.getEnd() + " ]";
			}

			// Show an error?
			if (isShowWarningIfParentDoesNotInclude()) System.err.println(err);
		}
	}

	/**
	 * Adjust [start,end] to include child
	 * @param child
	 */
	protected void adjust(Marker child) {
		start = Math.min(start, child.getStart());
		end = Math.max(end, child.getEnd());
	}

	/**
	 * Apply a SeqChange to a marker.
	 * 
	 * Create a new marker 
	 * 		newMarker = marker.apply( seqChange )
	 * 
	 * such that  
	 *		seqChange = Diff( newMarker , marker )		// Differences in seuquence 
	 * 
	 * @param seqChange
	 * @return A new marker after applying seqChange 
	 */
	public Marker apply(Variant seqChange) {
		// SeqChange after this marker: No effect
		if (end < seqChange.getStart()) return this;

		// We can only handle one change at a time
		if (seqChange.isVariantMultiple()) throw new RuntimeException("Cannot apply multiple changes!\n\tseqChange.isChangeMultiple() = " + seqChange.isVariantMultiple() + "\n\tSeqChange : " + seqChange);
		// Negative strand changes are a pain. We will eventually get rid of them...(they do not make sense any more)
		if (seqChange.isStrandMinus()) throw new RuntimeException("Only seqChenges in postive strand are accepted!\n\tSeqChange : " + seqChange);

		int lenChange = seqChange.lengthChange();
		if (lenChange == 0) return this;

		// SeqChange after marker end: Nothing to do
		if (end < seqChange.getStart()) return this;

		// We are not ready for mixed changes
		if (seqChange.isIns()) return applyIns(seqChange, lenChange);
		if (seqChange.isDel()) return applyDel(seqChange, lenChange);

		// TODO : Mixed changes are not supported
		throw new RuntimeException("Seqchange type not supported: " + seqChange.getChangeType() + "\n\t" + seqChange);
	}

	/**
	 * Apply a SeqChange to a marker. SeqChange is a deletion
	 * @param seqChange
	 * @return
	 */
	protected Marker applyDel(Variant seqChange, int lenChange) {
		Marker m = clone();

		// SeqChange Before start: Adjust coordinates
		if (seqChange.getEnd() < start) {
			m.start += lenChange;
			m.end += lenChange;
		} else if (seqChange.includes(this)) {
			// TODO: Should I create a zero length marker?
			//       This brings a few problems: 
			//			- We don't have any way to represent a zero length marker (if start=end then length is 1 base)
			//			- How does a zero length marker behaves in an INSERTION?
			//			- How does a zero length marker intersect?
			return null; // SeqChange completely includes this marker => The whole marker deleted
		} else if (includes(seqChange)) m.end += lenChange; // This marker completely includes seqChange. But seqChange does not include marker (i.e. they are not equal). Only 'end' coordinate needs to be updated
		else {
			// SeqChange is partially included in this marker.
			// This is treated as three different type of deletions:
			//		1- One before the marker
			//		2- One inside the marker
			//		3- One after the marker 
			// Note that type 1 and 3 cannot exists at the same time, otherwise the deletion would fully include the marker (previous case)

			// Deletion after the marker
			if (end < seqChange.getEnd()) {
				// Actually this does not affect the coordinates, so we don't care about this part
			}

			// Deletion matching the marker
			int istart = Math.max(start, m.getStart());
			int iend = Math.min(end, m.getEnd());
			if (iend < istart) throw new RuntimeException("This should never happen!"); // Sanity check

			end -= (iend - istart); // Update end coordinate

			// Deletion before the marker
			if (seqChange.getStart() < start) {
				// Update cooredinates shifting the marker to the left
				int delta = start - seqChange.getStart();
				m.start -= delta;
				m.end -= delta;
			}
		}

		return m;
	}

	/**
	 * Apply a SeqChange to a marker. SeqChange is an insertion
	 * @param seqChange
	 * @return
	 */
	public Marker applyIns(Variant seqChange, int lenChange) {
		Marker m = clone();

		if (seqChange.getStart() <= start) {
			// Insertion point before marker start? => Adjust both coordinates
			m.start += lenChange;
			m.end += lenChange;
		} else if (seqChange.getStart() <= end) {
			// Insertion point after start, but before end? => Adjust end coordinate
			m.end += lenChange;
		} else {
			// Insertion point after end, no effect on marker coordinates
		}

		return m;
	}

	@Override
	public Marker clone() {
		Marker m = (Marker) super.clone();
		if (m.getStart() != start) throw new RuntimeException("Cloned start does not match!");
		if (m.getEnd() != end) throw new RuntimeException("Cloned end does not match!");
		if (m.getStrand() != strand) throw new RuntimeException("Cloned strand does not match!");
		if (m.getParent() != parent) throw new RuntimeException("Cloned parent does not match!");
		return m;
	}

	/**
	 * Get a suitable codon table
	 * @return
	 */
	public CodonTable codonTable() {
		return CodonTables.getInstance().getTable(getGenome(), getChromosomeName());
	}

	/**
	 * Compare by start and end
	 */
	@Override
	public int compareTo(Interval i2) {
		// Compare chromosome names
		Marker m2 = (Marker) i2;

		Chromosome chr1 = getChromosome();
		Chromosome chr2 = m2.getChromosome();

		if ((chr1 != null) && (chr2 != null)) {
			// Non-null: Compare chromosomes
			int compChromo = chr1.compareChromoName(chr2);
			if (compChromo != 0) return compChromo;
		} else if ((chr1 == null) && (chr2 != null)) return 1; // One chromosome is null
		else if ((chr1 != null) && (chr2 == null)) return -1;

		// Compare by start position
		if (start > i2.start) return 1;
		if (start < i2.start) return -1;

		// Compare by end position
		if (end > i2.end) return 1;
		if (end < i2.end) return -1;

		// Compare by ID
		if ((id == null) && (i2.getId() == null)) return 0;
		if ((id != null) && (i2.getId() == null)) return -1;
		if ((id == null) && (i2.getId() != null)) return 1;
		return id.compareTo(i2.getId());
	}

	/**
	 * How far apart are these intervals?
	 * @return  Distance or -1 if they are not comparable (i.e. different chromosomes)
	 */
	public int distance(Marker interval) {
		if (!interval.getChromosomeName().equals(getChromosomeName())) return -1;

		if (intersects(interval)) return 0;

		if (start > interval.getEnd()) return start - interval.getEnd();
		if (interval.getStart() > end) return interval.getStart() - end;

		throw new RuntimeException("This should never happen!");
	}

	/**
	 * Distance from the beginning/end of a list of intervals, until this SNP
	 * It count the number of bases in 'markers'
	 * @param markers
	 * @return
	 */
	public int distanceBases(List<? extends Marker> markers, boolean fromEnd) {

		// Create a new list of sorted intervals
		ArrayList<Marker> markersSort = new ArrayList<Marker>();
		markersSort.addAll(markers);
		if (fromEnd) Collections.sort(markersSort, new IntervalComparatorByEnd(true));
		else Collections.sort(markersSort, new IntervalComparatorByStart());

		// Calculate distance
		int len = 0, latest = -1;
		for (Marker m : markersSort) {

			// Initialize
			if (latest < 0) {
				if (fromEnd) latest = m.getEnd() + 1;
				else latest = m.getStart() - 1;
			}

			if (fromEnd) {
				if (intersects(m)) return len + (m.getEnd() - start);
				else if (start > m.getEnd()) return len - 1 + (latest - start);

				latest = m.getStart();
			} else {
				if (intersects(m)) return len + (start - m.getStart());
				else if (start < m.getStart()) return len - 1 + (start - latest);

				latest = m.getEnd();
			}

			len += m.size();
		}

		if (fromEnd) return len - 1 + (latest - start);
		return len - 1 + (start - latest);
	}

	/**
	 * Go up (parent) until we find an instance of 'clazz'
	 */
	@SuppressWarnings("rawtypes")
	public Interval findParent(Class clazz) {
		if (this.getClass().equals(clazz)) return this;
		if ((parent != null) && (parent instanceof Marker)) return ((Marker) parent).findParent(clazz);
		return null;
	}

	public Chromosome getChromosome() {
		return (Chromosome) findParent(Chromosome.class);
	}

	/**
	 * Find chromosome name
	 * @return
	 */
	public String getChromosomeName() {
		Chromosome chromo = (Chromosome) findParent(Chromosome.class);
		if (chromo != null) return chromo.getId();
		return "";
	}

	/**
	 * Find chromosome and return it's number
	 * 
	 * @return Chromosome number if found, -1 otherwise
	 */
	public double getChromosomeNum() {
		Chromosome chromo = (Chromosome) findParent(Chromosome.class);
		if (chromo != null) return chromo.chromosomeNum;
		return -1;
	}

	/**
	 * Find genome 
	 * @return
	 */
	public Genome getGenome() {
		return (Genome) findParent(Genome.class);
	}

	/**
	 * Find genome name
	 * @return
	 */
	public String getGenomeName() {
		Genome genome = (Genome) findParent(Genome.class);
		if (genome != null) return genome.getId();
		return "";
	}

	@Override
	public Marker getParent() {
		return (Marker) parent;
	}

	public EffectType getType() {
		return type;
	}

	@Override
	public int hashCode() {
		int hashCode = getChromosomeName().hashCode();
		hashCode = hashCode * 31 + start;
		hashCode = hashCode * 31 + end;
		hashCode = hashCode * 31 + strand;
		if (id != null) hashCode = hashCode * 31 + id.hashCode();
		return hashCode;
	}

	public String idChain() {
		return idChain(";", true);
	}

	public String idChain(String separator, boolean useGeneId) {
		return idChain(separator, useGeneId, null);
	}

	/**
	 * A list of all IDs and parent IDs until chromosome
	 * @param m
	 * @return
	 */
	public String idChain(String separator, boolean useGeneId, ChangeEffect changeEffect) {
		StringBuilder sb = new StringBuilder();

		for (Marker m = this; (m != null) && !(m instanceof Chromosome) && !(m instanceof Genome); m = m.getParent()) {
			if (sb.length() > 0) sb.append(separator);

			switch (m.getType()) {
			case EXON:
				Transcript tr = (Transcript) m.getParent();
				Exon ex = (Exon) m;
				sb.append("exon_" + ex.getRank() + "_" + tr.numChilds() + "_" + ex.getSpliceType());
				break;

			case INTRON:
				Intron intron = (Intron) m;
				sb.append("intron_" + intron.getRank() + "_" + intron.getSpliceType());
				break;

			case GENE:
				Gene g = (Gene) m;
				sb.append(useGeneId ? m.getId() : g.getGeneName());
				sb.append(separator + g.getBioType());
				break;

			case TRANSCRIPT:
				sb.append(m.getId());
				sb.append(separator + ((Transcript) m).getBioType());
				break;

			case DOWNSTREAM:
				if ((changeEffect != null) && (changeEffect.getSeqChange() != null)) {
					Downstream downstream = (Downstream) m;
					sb.append(downstream.distanceToTr(changeEffect.getSeqChange()));
				}
				break;

			case UPSTREAM:
				if ((changeEffect != null) && (changeEffect.getSeqChange() != null)) {
					Upstream upstream = (Upstream) m;
					sb.append(upstream.distanceToTr(changeEffect.getSeqChange()));
				}
				break;

			case CHROMOSOME:
			case INTERGENIC:
				sb.append(m.getId());
				break;

			default:
				break;
			}
		}

		// Empty? Add ID
		if (sb.length() <= 0) sb.append(getId());

		// Prepend type
		sb.insert(0, this.getClass().getSimpleName() + separator);

		return sb.toString();
	}

	/**
	 * Is 'interval' completely included in 'this'?
	 * @return  return true if 'this' includes 'interval' 
	 */
	public boolean includes(Marker interval) {
		if (!interval.getChromosomeName().equals(getChromosomeName())) return false;
		return (start <= interval.start) && (interval.end <= end);
	}

	/**
	 * Intersect of two markers
	 * @param m
	 * @return A new marker which is the intersect of the two
	 */
	public Marker intersect(Marker m) {
		if (!getChromosomeName().equals(m.getChromosomeName())) return null;

		int istart = Math.max(start, m.getStart());
		int iend = Math.min(end, m.getEnd());
		if (iend < istart) return null;
		return new Marker(getParent(), istart, iend, strand, "");
	}

	/**
	 * Do the intervals intersect?
	 * @return  return true if this intersects 'interval' 
	 */
	public boolean intersects(Marker interval) {
		if (!interval.getChromosomeName().equals(getChromosomeName())) return false;
		return (interval.getEnd() >= start) && (interval.getStart() <= end);
	}

	/**
	 * How much do intervals intersect?
	 * @return  number of bases these intervals intersect
	 */
	public int intersectSize(Marker interval) {
		if (!interval.getChromosomeName().equals(getChromosomeName())) return 0;

		int start = Math.max(this.start, interval.getStart());
		int end = Math.min(this.end, interval.getEnd());

		if (end < start) return 0;
		return (end - start) + 1;
	}

	/**
	 * Adjust parent if it does not include child? 
	 * @return
	 */
	protected boolean isAdjustIfParentDoesNotInclude(Marker parent) {
		return false;
	}

	/**
	 * Show an error if parent does not include child? 
	 * @return
	 */
	protected boolean isShowWarningIfParentDoesNotInclude() {
		return false;
	}

	/**
	 * Return the difference between two markers
	 * 
	 * @param interval
	 * @return A set of 'markers'. Note that the result can have zero, one or two markers
	 */
	public Markers minus(Marker interval) {
		Markers ints = new Markers();
		if (intersects(interval)) {
			if ((interval.getStart() <= getStart()) && (getEnd() <= interval.getEnd())) {
				// 'this' is included in 'interval' => Nothing left
			} else if ((interval.getStart() <= getStart()) && (interval.getEnd() < getEnd())) {
				// 'interval' overlaps left part of 'this' => Include right part of 'this'
				ints.add(new Marker(getParent(), interval.getEnd() + 1, getEnd(), getStrand(), getId()));
			} else if ((getStart() < interval.getStart()) && (getEnd() <= interval.getEnd())) {
				// 'interval' overlaps right part of 'this' => Include left part of 'this'
				ints.add(new Marker(getParent(), getStart(), interval.getStart() - 1, getStrand(), getId()));
			} else if ((getStart() < interval.getStart()) && (interval.getEnd() < getEnd())) {
				// 'interval' overlaps middle of 'this' => Include left and right part of 'this'
				ints.add(new Marker(getParent(), getStart(), interval.getStart() - 1, getStrand(), getId()));
				ints.add(new Marker(getParent(), interval.getEnd() + 1, getEnd(), getStrand(), getId()));
			} else throw new RuntimeException("Interval intersection not analysed. This should nbever happen!");
		} else ints.add(this); // No intersection => Just add 'this' interval

		return ints;
	}

	/**
	 * Query all genomic regions that intersect 'marker' (this makes sense in Gene, Transcript, Exon, etc.)
	 */
	public Markers query(Marker marker) {
		return null;
	}

	/**
	 * Parse a line (form a file)
	 * Format: "chromosome \t start \t end \t id \n" 
	 */
	public void readTxt(String line, int lineNum, Genome genome, int positionBase) {
		line = line.trim(); // Remove spaces

		// Ignore empty lines and comment lines
		if ((line.length() > 0) && (!line.startsWith("#"))) {
			// Parse line
			String fields[] = line.split("\\s+");

			// Is line OK?
			if (fields.length >= 3) {
				Chromosome chromo = genome.getChromosome(fields[0].trim());
				if (chromo == null) System.err.println("WARNING: Chromosome '" + fields[0] + "' not found in genome '" + genome.getGenomeName() + "', version '" + genome.getVersion() + "'!\n\tLine: " + lineNum + "\t'" + line + "'");
				parent = chromo;
				start = Gpr.parseIntSafe(fields[1]) - positionBase;
				end = Gpr.parseIntSafe(fields[2]) - positionBase;
				id = "";

				if (fields.length >= 4) {
					// Join all ids using a space character (and remove all spaces)
					for (int t = 3; t < fields.length; t++)
						id += fields[t].trim() + " ";

					id = id.trim();
				}
			} else throw new RuntimeException("Error line " + lineNum + " (number of fields is " + fields.length + "):\t" + line);
		}
	}

	/**
	 * Calculate the effect of this seqChange
	 * @param seqChange : Sequence change
	 * @param changeEffect
	 * @return
	 */
	public boolean seqChangeEffect(Variant seqChange, ChangeEffects changeEffects) {
		if (!intersects(seqChange)) return false;
		changeEffects.add(this, type, "");
		return true;
	}

	/**
	 * Calculate the effect of this seqChange
	 * @param seqChange : Sequence change
	 * @param changeEffects
	 * @param seqChangeRef : Before analyzing results, we have to change markers using seqChangerRef to create a new reference 'on the fly'
	 * @return
	 */
	public boolean seqChangeEffect(Variant seqChange, ChangeEffects changeEffects, Variant seqChangerRef) {
		if (!intersects(seqChange)) return false;// Sanity check

		if (seqChangerRef != null) {
			Marker m = apply(seqChangerRef);

			// Has the marker been deleted?
			// Then there is no effect over this marker (it does not exist any more)
			if (m == null) return false;

			return m.seqChangeEffect(seqChange, changeEffects);
		}

		return seqChangeEffect(seqChange, changeEffects);
	}

	/**
	 * Parse a line from a serialized file
	 * @param line
	 * @return
	 */
	@Override
	public void serializeParse(MarkerSerializer markerSerializer) {
		type = EffectType.valueOf(markerSerializer.getNextField());
		markerSerializer.getNextFieldInt();
		parent = new MarkerParentId(markerSerializer.getNextFieldInt()); // Create a 'fake' parent. It will be replaced after all objects are in memory.
		start = markerSerializer.getNextFieldInt();
		end = markerSerializer.getNextFieldInt();
		id = markerSerializer.getNextField();
		strand = (byte) markerSerializer.getNextFieldInt();
	}

	/**
	 * Create a string to serialize to a file
	 * @return
	 */
	@Override
	public String serializeSave(MarkerSerializer markerSerializer) {
		return type //
				+ "\t" + markerSerializer.getIdByMarker(this) //
				+ "\t" + (parent != null ? markerSerializer.getIdByMarker((Marker) parent) : -1) //
				+ "\t" + start //
				+ "\t" + end //
				+ "\t" + id //
				+ "\t" + strand //
		;
	}

	/**
	 * To string as a simple "chr:start-end" format
	 * @return
	 */
	public String toStr() {
		return getClass().getSimpleName() + "_" + getChromosomeName() + ":" + (start + 1) + "-" + (end + 1);
	}

	@Override
	public String toString() {
		return getChromosomeName() + "\t" + start + "-" + end //
				+ " " //
				+ type + ((id != null) && (id.length() > 0) ? " '" + id + "'" : "");
	}

	/**
	 * Union of two markers
	 * @param m
	 * @return A new marker which is the union of the two 
	 */
	public Marker union(Marker m) {
		if (!getChromosomeName().equals(m.getChromosomeName())) return null;

		int ustart = Math.min(start, m.getStart());
		int uend = Math.max(end, m.getEnd());
		return new Marker(getParent(), ustart, uend, strand, "");
	}

}
