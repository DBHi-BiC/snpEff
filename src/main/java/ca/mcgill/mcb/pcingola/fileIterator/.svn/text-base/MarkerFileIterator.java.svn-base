package ca.mcgill.mcb.pcingola.fileIterator;

import java.io.BufferedReader;

import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * Opens a Marker file and iterates over all markers
 * 
 * @author pcingola
 */
public abstract class MarkerFileIterator<M extends Marker> extends FileIterator<M> {

	boolean createChromos = false; // Create chromosomes if not found
	protected Genome genome;
	protected boolean ignoreChromosomeErrors = true; // If true, do not throw an exception when a chromosome is not found. Just ignore the line
	protected int inOffset;

	public MarkerFileIterator(BufferedReader reader, int inOffset) {
		super(reader);
		this.inOffset = inOffset;
		this.genome = new Genome("genome");
		this.createChromos = true;
	}

	public MarkerFileIterator(String fileName, Genome genome, int inOffset) {
		super(fileName);
		this.inOffset = inOffset;
		this.genome = genome;
	}

	public MarkerFileIterator(String fileName, int inOffset) {
		super(fileName);
		this.inOffset = inOffset;
		this.genome = new Genome("genome");
		this.createChromos = true;
	}

	/**
	 * Find chromosome 'chromoName'. If it does not exists and 'createChromos' is true, the chromosome is created
	 * @param chromoName
	 * @return
	 */
	public Chromosome getChromosome(String chromoName) {
		if (createChromos) return genome.getOrCreateChromosome(chromoName);
		return genome.getChromosome(chromoName);
	}

	public Genome getGenome() {
		return genome;
	}

	public boolean isIgnoreChromosomeErrors() {
		return ignoreChromosomeErrors;
	}

	/**
	 * Parse a string as a 'position'.
	 * Note: It subtracts 'inOffset' so that all coordinates are zero-based
	 * 
	 * @param posStr
	 * @return
	 */
	public int parsePosition(String posStr) {
		return Gpr.parseIntSafe(posStr) - inOffset;
	}

	/**
	 * Sanity check
	 */
	public void sanityCheckChromo(String chromoName, Chromosome chromo) {
		if (chromo == null) {
			if (ignoreChromosomeErrors) {
				System.err.println("WARNING: Chromosome '" + chromoName + "' not found. File '" + fileName + "', line " + lineNum);
				return;
			}
			throw new RuntimeException("ERROR: Chromosome '" + chromoName + "' not found! File '" + fileName + "', line " + lineNum);
		}
	}

	public void setCreateChromos(boolean createChromos) {
		this.createChromos = createChromos;
	}

	public void setIgnoreChromosomeErrors(boolean ignoreChromosomeErrors) {
		this.ignoreChromosomeErrors = ignoreChromosomeErrors;
	}

	public void setInOffset(int inOffset) {
		this.inOffset = inOffset;
	}

}
