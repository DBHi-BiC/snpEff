package ca.mcgill.mcb.pcingola.snpEffect.factory;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import ca.mcgill.mcb.pcingola.fileIterator.FastaFileIterator;
import ca.mcgill.mcb.pcingola.interval.Cds;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.interval.IntervalComparatorByEnd;
import ca.mcgill.mcb.pcingola.interval.IntervalComparatorByStart;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.GprSeq;

/**
 * This class creates a SnpEffectPredictor from a file (or a set of files) and a configuration
 * 
 * @author pcingola
 */
public abstract class SnpEffPredictorFactory {

	// Show a mark every
	public static final int MARK = 1000;

	// Debug mode?
	public static boolean debug = false;

	boolean verbose = false;
	boolean readSequences = true; // Do not read sequences from GFF file (this is only used for debugging)
	String fileName;
	String fastaFile; // Only used for debugging or testing
	String line;
	int lineNum;
	Config config;
	Genome genome;
	SnpEffectPredictor snpEffectPredictor;
	int inOffset; // This amount is subtracted to all position coordinates
	int totalSeqsAdded = 0, totalSeqsIgnored = 0; // Number of sequences added and ignored
	HashMap<String, Integer> exonsByChromo;
	HashMap<String, Marker> markersById;
	HashMap<String, Gene> genesById;
	HashMap<String, Transcript> transcriptsById;

	public SnpEffPredictorFactory(Config config, int inOffset) {
		this.config = config;
		this.inOffset = inOffset;

		genome = config.getGenome();
		snpEffectPredictor = new SnpEffectPredictor(config.getGenome());
		exonsByChromo = new HashMap<String, Integer>();
		markersById = new HashMap<String, Marker>();
		genesById = new HashMap<String, Gene>();
		transcriptsById = new HashMap<String, Transcript>();
	}

	protected void add(Cds cds) {
		if (debug) Gpr.debug("Adding CDS:\t" + cds);
		Transcript tr = (Transcript) cds.getParent();
		tr.add(cds);

		addMarker(cds, false);
	}

	protected void add(Chromosome chromo) {
		genome.add(chromo);
	}

	protected void add(Exon exon) {
		if (debug) Gpr.debug("Adding Exon:\t" + exon);
		Transcript tr = (Transcript) exon.getParent();

		// Make sure the same exon was not added before
		Exon oldex = tr.get(exon.getId());
		if (oldex != null) {
			if (oldex.includes(exon)) return; // Redundant, just ignore it.

			// Create a new exon with same info and different 'id'
			exon = new Exon(tr, exon.getStart(), exon.getEnd(), exon.getStrand(), exon.getId() + "_" + tr.subintervals().size(), exon.getRank());
		}

		// Add exon
		tr.add(exon);
		addMarker(exon, false);
	}

	/**
	 * Add a Gene 
	 * @param gene
	 */
	protected void add(Gene gene) {
		if (debug) Gpr.debug("Adding Gene:\t" + gene.getId());
		snpEffectPredictor.add(gene);

		if (genesById.containsKey(gene.getId())) throw new RuntimeException("Gene  '" + gene.getId() + "' already exists");
		genesById.put(gene.getId(), gene);
	}

	/**
	 * Add a generic Marker
	 * @param marker
	 */
	protected void add(Marker marker) {
		addMarker(marker, false);
	}

	/**
	 * Add a transcript
	 * @param tr
	 */
	protected void add(Transcript tr) {
		if (debug) Gpr.debug("Adding Transcript:\t" + tr.getId());
		Gene gene = (Gene) tr.getParent();
		gene.add(tr);

		if (transcriptsById.containsKey(tr.getId())) throw new RuntimeException("Transcript  '" + tr.getId() + "' already exists");
		transcriptsById.put(tr.getId(), tr);
	}

	/**
	 * Add sequences to exon intervals
	 */
	protected void addExonSequences(String chr, String chrSeq) {
		int seqsAdded = 0, seqsIgnored = 0;
		System.out.print("\t\tAdding genomic sequences to exons: ");

		// Find and add sequences for all exons in this chromosome
		for (Gene gene : genome.getGenes()) {
			if (gene.getChromosomeName().equalsIgnoreCase(chr)) { // Same chromosome? => go on
				for (Transcript tr : gene) {
					for (Exon exon : tr) {
						int ssStart = exon.getStart();
						int ssEnd = exon.getEnd() + 1; // String.substring does not include the last character in the interval (so we have to add 1)

						if ((ssStart < 0) || (ssEnd > chrSeq.length())) {
							System.err.println("Ignoring exon outside chromosome range (chromo length: " + chrSeq.length() + "). Exon: " + exon);
							seqsIgnored++;
						} else {
							try {
								String seq = chrSeq.substring(ssStart, ssEnd).toUpperCase();
								// Reverse strand? => reverse complement of the sequence
								if (exon.isStrandMinus()) seq = GprSeq.reverseWc(seq);
								exon.setSequence(seq);
								seqsAdded++;
							} catch (Throwable t) {
								t.printStackTrace();
								throw new RuntimeException("Error trying to add sequence to exon:\n\tChromosome sequence length: " + chrSeq.length() + "\n\tExon: " + exon);
							}
						}
					}
				}
			}
		}

		System.out.println("\tDone (" + seqsAdded + " sequences added, " + seqsIgnored + " ignored).");
		totalSeqsAdded += seqsAdded;
		totalSeqsIgnored += seqsIgnored;
	}

	/**
	 * Add a marker to the collection
	 * @param marker
	 */
	protected void addMarker(Marker marker, boolean unique) {
		String key = marker.getId();
		if (unique && markersById.containsKey(key)) throw new RuntimeException("Marker '" + key + "' already exists");
		markersById.put(key, marker);
	}

	/**
	 * Adjust chromosome length using gene information
	 * This is used when the sequence is not available (which makes sense on test-cases and debugging only)
	 */
	protected void adjustChromosomes() {
		HashMap<String, Integer> lenByChr = new HashMap<String, Integer>();
		for (Gene gene : config.getGenome().getGenes()) {
			String chrName = gene.getChromosomeName();
			Integer len = lenByChr.get(chrName);
			int max = Math.max(gene.getEnd(), (len != null ? len : 0));
			lenByChr.put(chrName, max);
		}

		// Set length
		for (String chrName : lenByChr.keySet()) {
			config.getGenome().getChromosome(chrName).setEnd(lenByChr.get(chrName));
		}
	}

	/**
	 * Perform some actions before reading sequences
	 */
	protected void beforeExonSequences() {
		// Sometimes we have to guess exon info from CDS info (not the best case scenario, but there are a lot of crappy genome annotations around)
		exonsFromCds();

		// Some annotation formats split exons in two parts (e.g. stop-codon not part of exon in GTF).
		deleteRedundant();

		// Some annotations introduce zero size introns
		collapseZeroGap();
	}

	/**
	 * Get (or create) a chromosome and set it's length
	 * @param chromoName
	 * @param len
	 */
	void chromoLen(String chromoName, int len) {
		Chromosome chromo = getOrCreateChromosome(chromoName);
		chromo.setLength(len);
	}

	/**
	 * Collapse exons having zero size introns between them
	 */
	protected void collapseZeroGap() {
		System.out.print("\n\tCollapsing zero length introns (if needed): ");

		int count = 0;
		for (Gene gene : genome.getGenes())
			for (Transcript tr : gene)
				if (tr.collapseZeroGap()) {
					count++;
					System.out.print('.');
				}

		System.out.println("\n\tTotal collapsed transcripts: " + count);
	}

	/**
	 * Count number of exons by chromosome
	 */
	@SuppressWarnings("unused")
	void countExonsByChromo() {
		exonsByChromo = new HashMap<String, Integer>();

		for (Gene gint : genome.getGenes()) {
			Chromosome chromo = gint.getChromosome();
			for (Transcript tint : gint) {
				for (Exon eint : tint) {
					// Get current count
					String chromoName = chromo.getId();
					Integer count = exonsByChromo.get(chromoName);

					// Increment
					if (count == null) count = 1;
					else count++;

					// Store
					exonsByChromo.put(chromoName, count);
				}
			}
		}
	}

	public abstract SnpEffectPredictor create();

	/**
	 * Consolidate transcripts: 
	 * If two exons are one right next to the other, join them
	 * E.g. exon1:1234-2345, exon2:2346-2400 => exon:1234-2400
	 * This happens mostly in GTF files, where the stop-codon is specified separated from the exon info.
	 */
	protected void deleteRedundant() {
		System.out.print("\n\tDeleting redundant exons (if needed): ");
		int count = 0;
		for (Gene gene : genome.getGenes())
			for (Transcript tr : gene)
				if (tr.deleteRedundant()) {
					count++;
					System.out.print('.');
				}

		System.out.println("\n\tTotal transcripts with deleted exons: " + count);
	}

	/**
	 * Error: Throw a runtime exception (show some details)
	 * @param msg
	 */
	void error(String msg) {
		throw new RuntimeException("FATAL ERROR: " + msg + ". File '" + fileName + "' line " + lineNum + "\n\t'" + line + "'\n");
	}

	/**
	 * Create exons from CDS info
	 */
	protected void exonsFromCds() {
		System.out.print("\n\tCreate exons from CDS (if needed): ");

		int count = 0;
		for (Gene gene : genome.getGenes()) {
			for (Transcript tr : gene) {
				// CDS length
				int lenCds = 0;
				for (Cds cds : tr.getCds())
					lenCds += cds.size();

				// Exon length
				int lenExons = 0;
				for (Exon ex : tr)
					lenExons += ex.size();

				// Cds length larger than exons? => something is missing
				if (lenCds > lenExons) {
					exonsFromCds(tr);
					count++;
				}
			}
		}
		System.out.println("\n\tExons created for " + count + " transcripts.");
	}

	/**
	 * Create exons from CDS info
	 * WARNING: We might end up with redundant exons if some exons existed before this process
	 * 
	 * @param tr : Transcript with CDS info, but no exons
	 */
	protected void exonsFromCds(Transcript tr) {
		List<Cds> cdss = tr.getCds();

		// First: Check and adjust strand info
		int trStrand = tr.getStrand() >= 0 ? 1 : -1;
		int cdsStrand = 0;
		for (Cds cds : cdss)
			cdsStrand += cds.getStrand();
		cdsStrand = cdsStrand >= 0 ? 1 : -1;
		if (cdsStrand != trStrand) {
			System.out.print(cdsStrand >= 0 ? '+' : '-');
			tr.setStrand(cdsStrand);
		}

		// Sort CDS by strand
		if (tr.getStrand() >= 0) Collections.sort(cdss, new IntervalComparatorByStart()); // Sort by start position 
		else Collections.sort(cdss, new IntervalComparatorByEnd(true)); // Sort by end position (reversed) 

		// Add cds as exons
		// WARNING: We might end up with redundant exons if some exons existed before this process
		int rank = 1;
		for (Cds cds : cdss) {
			// Create exon and add it to transcript
			String id = "Exon_" + cds.getChromosomeName() + "_" + cds.getStart() + "_" + cds.getEnd();
			if (tr.get(id) == null) { // Don't add an exon twice
				Exon exon = new Exon(tr, cds.getStart(), cds.getEnd(), trStrand, id, rank);
				tr.add(exon);
			}

			rank++;
			System.out.print('.');
		}
	}

	protected Gene findGene(String id) {
		Gene gene = genesById.get(id);
		if (gene != null) return gene;
		return genesById.get("Gene_" + id); // Alternative gene ID
	}

	protected Gene findGene(String geneId, String id) {
		Gene gene = findGene(geneId);
		if (gene != null) return gene;
		return genesById.get("Gene_" + id); // Alternative gene ID
	}

	protected Marker findMarker(String id) {
		return markersById.get(id);
	}

	protected Transcript findTranscript(String id) {
		Transcript tr = transcriptsById.get(id);
		if (tr != null) return tr;
		return transcriptsById.get("Transcript_" + id); // Alternative transcript ID
	}

	protected Transcript findTranscript(String trId, String id) {
		Transcript tr = findTranscript(trId);
		if (tr != null) return tr;
		return transcriptsById.get("Transcript_" + id);
	}

	/**
	 * Finish up procedure to ensure consistency
	 */
	void finishUp(boolean verbose) {
		int showEvery = 100;

		// Adjust genes: recalculate start, end, strand, etc.
		int i = 1;
		System.out.print("\n\tAdjusting genes: ");
		for (Gene gene : genome.getGenes())
			if (gene.adjust()) Gpr.showMark(i++, showEvery);

		// Adjust chromosome sizes
		System.out.print("\n\tAdjusting chromosome sizes: ");
		for (Gene gene : genome.getGenes()) {
			Chromosome chr = gene.getChromosome();

			if (gene.getEnd() > chr.getEnd()) {
				chr.setLength(gene.getEnd() + 1);
				Gpr.showMark(i++, showEvery);
			}
		}

		// Remove suspicious transcripts
		removerSuspiciousTranscripts(showEvery);

		// Adjust transcripts: recalculate start, end, strand, etc.
		i = 1;
		System.out.print("\n\tAdjusting transcripts: ");
		for (Gene gene : genome.getGenes())
			for (Transcript tr : gene)
				if (tr.adjust()) Gpr.showMark(i++, showEvery);

		// Adjust exons: Most file formats don't have exon rank information.
		i = 1;
		System.out.print("\n\tRanking exons: ");
		for (Gene gene : genome.getGenes())
			for (Transcript tr : gene)
				if (tr.rankExons()) Gpr.showMark(i++, showEvery);

		// If some UTRs are missing: calculate UTR information from CDS whenever possible
		System.out.print("\n\tCreate UTRs from CDS (if needed): ");
		utrFromCds(verbose, showEvery);

		// Remove empty chromosomes
		removeEmptyChromos();

		// Done
		System.out.println("");
	}

	/**
	 * Get a chromosome. If it doesn't exist, create it
	 * @param chromoName
	 * @return
	 */
	protected Chromosome getOrCreateChromosome(String chromoName) {
		Chromosome chromo = genome.getChromosome(chromoName);

		// Not found? => Create a new one
		if (chromo == null) {
			chromo = new Chromosome(genome, 0, 0, 1, chromoName);
			genome.add(chromo);
		}

		return chromo;
	}

	/**
	 * Does this chromosome have any exons?
	 * @param chromoName
	 * @return
	 */
	boolean hasExons(String chromoName) {
		Integer count = exonsByChromo.get(chromoName);
		return (count != null) && (count > 0);
	}

	/**
	 * Parse a string as a 'position'.
	 * Note: It subtracts 'inOffset' so that all coordinates are zero-based
	 * 
	 * @param posStr
	 * @return
	 */
	protected int parsePosition(String posStr) {
		return Gpr.parseIntSafe(posStr) - inOffset;
	}

	/**
	 * Read exon sequences from a FASTA file
	 * @param fastaFile
	 */
	protected void readExonSequences() {
		List<String> files = config.getFileListGenomeFasta();

		// Force a specific file?
		if (fastaFile != null) {
			files.clear();
			files.add(fastaFile);
		}

		// Try all files in the list until one is available
		for (String file : files) {
			System.out.println("\t\tTrying FASTA file: '" + file + "'");

			if (Gpr.canRead(file)) {
				// Read fasta sequence
				FastaFileIterator ffi = new FastaFileIterator(file);
				for (String seq : ffi) {
					String chromo = ffi.getName();
					System.out.println("\t\tReading sequence '" + chromo + "', length: " + seq.length());
					Chromosome chromoInt = getOrCreateChromosome(chromo);
					chromoInt.setLength(seq.length()); // Set chromosome length
					addExonSequences(chromo, seq); // Add all sequences
				}
				return;
			}
		}

		throw new RuntimeException("Cannot find reference sequence.");
	}

	/**
	 * Remove empty chromosomes
	 */
	void removeEmptyChromos() {
		System.out.println("\n\tRemove empty chromosomes: ");
		ArrayList<Chromosome> chrToDelete = new ArrayList<Chromosome>();
		for (Chromosome chr : config.getGenome())
			if (chr.size() <= 1) chrToDelete.add(chr);

		for (Chromosome chr : chrToDelete) {
			System.out.println("\t\tRemoving empty chromosome: '" + chr.getId() + "'");
			config.getGenome().remove(chr);
		}

		if (chrToDelete.size() > 0) {
			System.out.print("\t\tChromosome left: ");
			for (Chromosome chr : config.getGenome())
				System.out.print(chr.getId() + " ");
			System.out.println("");
		}
	}

	/** 
	 * Remove suspicious transcripts. 
	 * Transcripts whose first exon has a non-zero frame are removed. According 
	 * to ENSEMBL, there is no enough evidence to define these transcripts properly.
	 * 
	 * @param showEvery
	 */
	void removerSuspiciousTranscripts(int showEvery) {
		// Update Exon.frame data (if not available)
		System.out.print("\n\tUpdating frame info from CDSs to Exons.");
		int i = 0;
		for (Gene gene : genome.getGenes())
			for (Transcript tr : gene) {
				for (Exon ex : tr) {
					// No frame info? => try to find matching CDS
					if (ex.getFrame() < 0) {
						for (Cds cds : tr.getCds()) {
							// Found same cds as exon? => Copy frame info
							if ((ex.getStart() == cds.getStart()) && (ex.getEnd() == cds.getEnd())) ex.setFrame(cds.getFrame());
						}
					}
				}
			}

		// Mark Transcripts for removal
		i = 1;
		System.out.print("\n\tFiltering out suspicious transcripts (first exon has non-zero frame): ");
		LinkedList<Transcript> trToDelete = new LinkedList<Transcript>();
		for (Gene gene : genome.getGenes())
			for (Transcript tr : gene) {
				List<Exon> exons = tr.sortedStrand();
				if ((exons != null) && !exons.isEmpty()) {
					Exon exon = exons.get(0); // Get first exon
					if (exon.getFrame() > 0) {
						trToDelete.add(tr); // First exon is not zero? => Mark for deletion
						Gpr.showMark(i++, showEvery);
					}
				}
			}

		// Remove transcripts
		for (Transcript tr : trToDelete) {
			Gene gene = (Gene) tr.getParent();
			gene.remove(tr);
		}

		System.out.print((trToDelete.size() > 0 ? "\n" : "") + "\tTotal: " + trToDelete.size() + " removed.");
	}

	public void setFastaFile(String fastaFile) {
		this.fastaFile = fastaFile;
	}

	public void setFileName(String fileName) {
		this.fileName = fileName;
	}

	/**
	 * Read sequences?
	 * Note: This is only used for debugging and testing
	 */
	public void setReadSequences(boolean readSequences) {
		this.readSequences = readSequences;
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	String unquote(String qstr) {
		return qstr.replaceAll("\"", "");
	}

	/**
	 * Create missing UTRs from CDS information
	 */
	void utrFromCds(boolean verbose, int showEvery) {
		int i = 1;
		for (Gene gint : genome.getGenes())
			for (Transcript tint : gint) {
				boolean show = tint.utrFromCds(verbose);
				if (show && !verbose) Gpr.showMarkStderr(i++, showEvery);
			}
	}

	/**
	 * Warning: Show a warning message (show some details)
	 * @param msg
	 */
	void warning(String msg) {
		if (verbose) System.err.println("WARNING: " + msg + ". File '" + fileName + "' line " + lineNum + "\t'" + line + "'");
	}

}
