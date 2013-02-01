package ca.mcgill.mcb.pcingola.interval;

import java.util.List;

import ca.mcgill.mcb.pcingola.codons.CodonTable;
import ca.mcgill.mcb.pcingola.codons.CodonTables;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * An interval intended as a mark
 * 
 * @author pcingola
 */
public class Marker extends Interval {

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
					+ "\t" + this.getClass().getSimpleName() + " '" + this.getId() + "'  [ " + this.getStart() + " , " + this.getEnd() + " ]" //
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
		hashCode = hashCode * 31 + id.hashCode();
		return hashCode;
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
	 * how much do intervals intersect?
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
	 * @param seqChange
	 * @param changeEffect
	 * @return
	 */
	public List<ChangeEffect> seqChangeEffect(SeqChange seqChange, ChangeEffect changeEffect) {
		if (!intersects(seqChange)) return ChangeEffect.emptyResults(); // Sanity check
		changeEffect.set(this, type, "");
		return changeEffect.newList();
	}

	/**
	 * Parse a line from a serialized file
	 * @param line
	 * @return
	 */
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
