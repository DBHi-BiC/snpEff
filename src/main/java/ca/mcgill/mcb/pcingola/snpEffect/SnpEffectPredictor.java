package ca.mcgill.mcb.pcingola.snpEffect;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import net.sf.samtools.util.RuntimeEOFException;
import ca.mcgill.mcb.pcingola.interval.Cds;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.interval.Intergenic;
import ca.mcgill.mcb.pcingola.interval.Intron;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.MarkerSerializer;
import ca.mcgill.mcb.pcingola.interval.Markers;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.SpliceSite;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.interval.Utr;
import ca.mcgill.mcb.pcingola.interval.tree.IntervalForest;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.ErrorType;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * Predicts effects of SNPs 
 * 
 * Note: Actually tries to predict any kind of SeqChange, not only SNPs . It is called SnpEffectPredictor for 'historical reasons'.
 * 
 * @author pcingola
 *
 */
public class SnpEffectPredictor implements Serializable {
	private static final long serialVersionUID = 4519418862303325081L;

	public static final int DEFAULT_UP_DOWN_LENGTH = 5000;

	boolean useChromosomes = true;
	int upDownStreamLength = DEFAULT_UP_DOWN_LENGTH;
	Genome genome;
	Markers markers; // All other markers are stored here (e.g. custom markers, intergenic, etc.)
	IntervalForest intervalForest;

	/**
	 * Load predictor from a binary file
	 */
	public static SnpEffectPredictor load(Config config) {
		String snpEffPredFile = config.getFileSnpEffectPredictor();

		// Sanity check
		if (!Gpr.canRead(snpEffPredFile)) throw new RuntimeException("\tERROR: Cannot read file '" + snpEffPredFile + "'.\n\tYou can try to download the database by running the following command:\n\t\tjava -jar snpEff.jar download " + config.getGenome().getVersion() + "\n");

		// Load markers from file
		MarkerSerializer ms = new MarkerSerializer();
		Markers markers = ms.load(snpEffPredFile);

		// Find genome
		Genome genome = null;
		for (Marker m : markers)
			if (m instanceof Genome) genome = (Genome) m;
		if (genome == null) throw new RuntimeException("Genome not found. This should never happen!");

		// Create predictor
		SnpEffectPredictor snpEffectPredictor = new SnpEffectPredictor(genome);

		// Add genes
		for (Marker m : markers)
			if (m instanceof Gene) {
				Gene gene = (Gene) m;
				snpEffectPredictor.add(gene);
			}

		// Add 'other' markers
		for (Marker m : markers)
			if (!(m instanceof Genome) //
					&& !(m instanceof Chromosome) //
					&& !(m instanceof Gene) //
					&& !(m instanceof Transcript) //
					&& !(m instanceof Exon) //
					&& !(m instanceof Cds) //
					&& !(m instanceof Utr) //
					&& !(m instanceof SpliceSite) //
			) snpEffectPredictor.add(m);

		return snpEffectPredictor;
	}

	public SnpEffectPredictor(Genome genome) {
		this.genome = genome;
		markers = new Markers();
	}

	/**
	 * Add a gene interval
	 * @param gene
	 */
	public void add(Gene gene) {
		genome.getGenes().add(gene);
	}

	/** 
	 * Add a marker
	 * 
	 * Note: Markers have to be added BEFORE building the interval trees. 
	 *       Interval trees are built the first time you call snpEffect(snp) method.
	 * 
	 * @param marker
	 */
	public void add(Marker marker) {
		markers.add(marker);
	}

	/**
	 * Create interval trees (forest)
	 */
	public void buildForest() {
		intervalForest = new IntervalForest();

		// Add all chromosomes to forest
		if (useChromosomes) {
			for (Chromosome chr : genome)
				intervalForest.add(chr);
		}

		// Add all genes to forest
		for (Gene gene : genome.getGenes())
			intervalForest.add(gene);

		//---
		// Add to markers to 'markers'
		//---
		// Add up-down stream intervals
		for (Marker upDownStream : genome.getGenes().createUpDownStream(upDownStreamLength))
			add(upDownStream);

		// Add splice site intervals
		for (Marker spliceSite : genome.getGenes().findSpliceSites(true))
			add(spliceSite);

		// Intergenic markers
		for (Intergenic intergenic : genome.getGenes().createIntergenic())
			add(intergenic);

		intervalForest.add(markers); // Add all 'markers' to forest (includes custom intervals)

		// Build interval forest
		intervalForest.build();
	}

	/**
	 * Obtain a gene interval
	 * @param geneIntervalId
	 * @return
	 */
	public Gene getGene(String geneIntervalId) {
		return genome.getGenes().get(geneIntervalId);
	}

	public Genome getGenome() {
		return genome;
	}

	public IntervalForest getIntervalForest() {
		return intervalForest;
	}

	public Markers getMarkers() {
		return markers;
	}

	public int getUpDownStreamLength() {
		return upDownStreamLength;
	}

	/**
	 * Return a collection of intervals thet intercept marker
	 */
	public Markers intersects(Marker marker) {
		return intervalForest.query(marker);
	}

	/**
	 * Is the chromosome missing in this marker?
	 * @param marker
	 * @return
	 */
	boolean isChromosomeMissing(Marker marker) {
		// Missing chromosome in marker?
		if (marker.getChromosome() == null) return true;

		// Missing chromosome in genome?
		String chrName = marker.getChromosomeName();
		Chromosome chr = genome.getChromosome(chrName);
		if (chr == null) return true;

		// Chromosome length is 1 or less?
		if (chr.size() < 1) return true;

		// Tree not found in interval forest?
		if (!intervalForest.hasTree(chrName)) return true;

		// OK, we have the chromosome
		return false;
	}

	/**
	 * Remove all non-canonical transcripts
	 * @return : Number of transcripts removed
	 */
	public int keepTranscripts(Set<String> trIds) {
		int total = 0;
		for (Gene g : genome.getGenes())
			total += g.keepTranscripts(trIds);
		return total;
	}

	/**
	 * Dump to sdtout
	 */
	public void print() {
		System.out.println(genome);

		// Show genes
		for (Gene gene : genome.getGenes().sorted())
			System.out.println(gene);

		// Show other inervals
		for (Marker marker : markers)
			System.out.println(marker);
	}

	/**
	 * Name of the regions hit by a marker
	 * @param marker
	 * @return A set of region names
	 */
	public Set<String> regions(Marker marker, boolean showGeneDetails, boolean compareTemplate) {
		return regions(marker, showGeneDetails, compareTemplate, null);
	}

	/**
	 * Name of the regions hit by a marker
	 * @param marker
	 * @param showGeneDetails
	 * @param compareTemplate
	 * @param id : Only use genes or transcripts matching this ID
	 * @return
	 */
	public Set<String> regions(Marker marker, boolean showGeneDetails, boolean compareTemplate, String id) {
		if (Config.get().isErrorOnMissingChromo() && isChromosomeMissing(marker)) throw new RuntimeEOFException("Chromosome missing for marker: " + marker);

		boolean hitChromo = false;
		HashSet<String> hits = new HashSet<String>();

		Markers intersects = intersects(marker);
		if (intersects.size() > 0) {
			for (Marker markerInt : intersects) {

				if (markerInt instanceof Chromosome) {
					hitChromo = true; // OK (we have to hit a chromosome, otherwise it's an error
					hits.add(markerInt.getClass().getSimpleName()); // Add marker name to the list
				} else if (markerInt instanceof Gene) {
					// Analyze Genes
					Gene gene = (Gene) markerInt;
					regionsAddHit(hits, gene, marker, showGeneDetails, compareTemplate);

					// For all transcripts...
					for (Transcript tr : gene) {
						if ((id == null) || gene.getId().equals(id) || tr.getId().equals(id)) { // Mathes ID? (...or no ID to match)

							// Does it intersect this transcript?
							if (tr.intersects(marker)) {
								regionsAddHit(hits, tr, marker, showGeneDetails, compareTemplate);

								// Does it intersect a UTR? 
								for (Utr utr : tr.getUtrs())
									if (utr.intersects(marker)) regionsAddHit(hits, utr, marker, showGeneDetails, compareTemplate);

								// Does it intersect an exon?
								for (Exon ex : tr)
									if (ex.intersects(marker)) regionsAddHit(hits, ex, marker, showGeneDetails, compareTemplate);

								// Does it intersect an intron?
								for (Intron intron : tr.introns())
									if (intron.intersects(marker)) regionsAddHit(hits, intron, marker, showGeneDetails, compareTemplate);
							}
						}
					}
				} else {
					// No ID to match?
					if (id == null) regionsAddHit(hits, markerInt, marker, showGeneDetails, compareTemplate);
					else {
						// Is ID from transcript?
						Transcript tr = (Transcript) markerInt.findParent(Transcript.class);
						if ((tr != null) && (tr.getId().equals(id))) {
							regionsAddHit(hits, markerInt, marker, showGeneDetails, compareTemplate); // Transcript ID matches => count
						} else {
							// Is ID from gene?
							Gene gene = (Gene) markerInt.findParent(Gene.class);
							if ((gene != null) && (gene.getId().equals(id))) regionsAddHit(hits, markerInt, marker, showGeneDetails, compareTemplate); // Gene ID matches => count
						}
					}
				}
			}
		}

		if (!hitChromo) throw new RuntimeException("ERROR: Out of chromosome range. " + marker);
		return hits;
	}

	/**
	 * Add into to a hash
	 * @param hits
	 * @param marker
	 * @param hit2add
	 * @param showGeneDetails
	 * @param compareTemplate
	 */
	void regionsAddHit(HashSet<String> hits, Marker hit2add, Marker marker, boolean showGeneDetails, boolean compareTemplate) {
		String hitStr = hit2add.getClass().getSimpleName();

		if (compareTemplate) {
			Gene gene = (Gene) hit2add.findParent(Gene.class);
			if (gene != null) hitStr += (hit2add.isStrandPlus() == marker.isStrandPlus()) ? "_TEMPLATE_STRAND" : "_NON_TEMPLATE_STRAND";
		}

		if (showGeneDetails && (hit2add instanceof Gene)) {
			Gene gene = (Gene) hit2add;
			hitStr += "[" + gene.getBioType() + ", " + gene.getGeneName() + ", " + (gene.isProteinCoding() ? "protein" : "not-protein") + "]";
		}

		hits.add(hitStr); // Add marker name to the list
	}

	/**
	 * Regions hit by a marker
	 * @return
	 */
	public Set<Marker> regionsMarkers(Marker marker) {
		if (Config.get().isErrorOnMissingChromo() && isChromosomeMissing(marker)) throw new RuntimeEOFException("Chromosome missing for marker: " + marker);

		boolean hitChromo = false;
		HashSet<Marker> hits = new HashSet<Marker>();

		Markers intersects = intersects(marker);
		if (intersects.size() > 0) {
			for (Marker markerInt : intersects) {
				hits.add(markerInt);

				if (markerInt instanceof Chromosome) {
					hitChromo = true; // OK (we have to hit a chromosome, otherwise it's an error
				} else if (markerInt instanceof Gene) {
					// Analyze Genes
					Gene gene = (Gene) markerInt;

					// For all transcripts...
					for (Transcript tr : gene) {

						if (tr.intersects(marker)) { // Does it intersect this transcript?
							hits.add(tr);

							for (Utr utr : tr.getUtrs())
								if (utr.intersects(marker)) hits.add(utr);

							for (Exon ex : tr)
								if (ex.intersects(marker)) hits.add(ex);

							for (Intron intron : tr.introns())
								if (intron.intersects(marker)) hits.add(intron);
						}
					}
				}
			}
		}

		if (!hitChromo) throw new RuntimeException("ERROR: Out of chromosome range. " + marker);
		return hits;
	}

	/**
	 * Remove all non-canonical transcripts
	 */
	public void removeNonCanonical() {
		for (Gene g : genome.getGenes())
			g.removeNonCanonical();
	}

	/**
	 * Save predictor to a binary file (specified by the configuration)
	 */
	public void save(Config config) {
		String databaseFile = config.getFileSnpEffectPredictor();
		MarkerSerializer markerSerializer = new MarkerSerializer();
		markerSerializer.save(databaseFile, this);
	}

	/**
	 * Predict the effect of a seqChange
	 * @param seqChange
	 */
	public List<ChangeEffect> seqChangeEffect(SeqChange seqChange) {
		// No change? => Nothing to predict
		if ((!seqChange.isChange()) && (!seqChange.isInterval())) return ChangeEffect.emptyResults();

		ChangeEffect results = new ChangeEffect(seqChange);

		// Chromosome missing?
		if (Config.get().isErrorOnMissingChromo() && isChromosomeMissing(seqChange)) {
			results.addError(ErrorType.ERROR_CHROMOSOME_NOT_FOUND);
			return results.newList();
		}

		// Which intervals does seqChange intersect?
		Markers intersects = intersects(seqChange);

		// Show all results
		boolean hitChromo = false, hitSomething = false;
		ArrayList<ChangeEffect> resultsList = new ArrayList<ChangeEffect>();
		if (intersects.size() > 0) {
			for (Marker marker : intersects) {
				if (marker instanceof Chromosome) hitChromo = true; // Do we hit any chromosome?
				else { // Analyze all markers 
					results = new ChangeEffect(seqChange);
					List<ChangeEffect> resList = marker.seqChangeEffect(seqChange, results);
					if (!resList.isEmpty()) resultsList.addAll(resList);
					hitSomething = true;
				}
			}
		}

		// Any errors or intergenic (i.e. did not hit any gene)
		if (!hitChromo) {
			if (Config.get().isErrorChromoHit()) {
				results.addError(ErrorType.ERROR_OUT_OF_CHROMOSOME_RANGE);
				return results.newList();
			}
		} else if (!hitSomething) {
			if (Config.get().isOnlyRegulation()) {
				results.effectType = EffectType.NONE;
				return results.newList();
			} else {
				results.effectType = EffectType.INTERGENIC;
				return results.newList();
			}
		}

		return resultsList;
	}

	public void setUpDownStreamLength(int upDownStreamLength) {
		this.upDownStreamLength = upDownStreamLength;
	}

	public void setUseChromosomes(boolean useChromosomes) {
		this.useChromosomes = useChromosomes;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(genome.getVersion() + "\n");
		for (Chromosome chr : genome)
			sb.append(chr + "\n");
		sb.append(genome.getGenes());
		return sb.toString();
	}

}
