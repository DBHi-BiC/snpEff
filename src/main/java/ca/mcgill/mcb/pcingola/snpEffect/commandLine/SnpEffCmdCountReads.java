package ca.mcgill.mcb.pcingola.snpEffect.commandLine;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import net.sf.samtools.AbstractBAMFileIndex;
import net.sf.samtools.BAMIndexMetaData;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.interval.Intron;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.stats.CountByKey;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Count reads from a BAM file given a list of intervals
 * 
 * @author pablocingolani
 */
public class SnpEffCmdCountReads extends SnpEff {

	public static int SHOW_EVERY = 10000;
	public static boolean debug = true;

	List<String> samFileNames;
	ArrayList<CountByKey<Marker>> countReadsByFile;
	ArrayList<CountByKey<Marker>> countBasesByFile;
	SnpEffectPredictor snpEffectPredictor;
	boolean verbose = false;

	public SnpEffCmdCountReads() {
		samFileNames = new ArrayList<String>();
		countReadsByFile = new ArrayList<CountByKey<Marker>>();
		countBasesByFile = new ArrayList<CountByKey<Marker>>();
	}

	/**
	 * Count all reads in a BAM file 
	 * Note: It uses the BAM index
	 * 
	 * @param samReader
	 * @return
	 */
	int countTotalReads(String samFileName) {
		try {
			if (verbose) Timer.showStdErr("Counting reads on file: " + samFileName);
			SAMFileReader samReader = new SAMFileReader(new File(samFileName));
			AbstractBAMFileIndex index = (AbstractBAMFileIndex) samReader.getIndex();
			int count = 0;
			for (int i = 0; i < index.getNumberOfReferences(); i++) {
				BAMIndexMetaData meta = index.getMetaData(i);
				count += meta.getAlignedRecordCount();
			}
			samReader.close();
			if (verbose) Timer.showStdErr("Total " + count + " reads.");
			return count;
		} catch (Exception e) {
			// Error? (e.g. no index)
			System.err.println("ERROR! BAM file not indexed?");
			return -1;
		}
	}

	/**
	 * A list of all IDs and parent IDs until chromosome
	 * @param m
	 * @return
	 */
	String idChain(Marker marker) {
		StringBuilder sb = new StringBuilder();

		for (Marker m = marker; (m != null) && !(m instanceof Chromosome) && !(m instanceof Genome); m = m.getParent()) {

			switch (m.getType()) {
			case EXON:
				if (sb.length() > 0) sb.append(";");
				sb.append("exon_" + ((Exon) m).getRank());
				break;

			case INTRON:
				if (sb.length() > 0) sb.append(";");
				sb.append("intron_" + ((Intron) m).getRank());
				break;

			case CHROMOSOME:
			case INTERGENIC:
			case GENE:
			case TRANSCRIPT:
				if (sb.length() > 0) sb.append(";");
				sb.append(m.getId());
				break;

			default:
				break;
			}
		}

		// Empty? Add ID
		if (sb.length() <= 0) sb.append(marker.getId());

		// Prepend type
		sb.insert(0, marker.getClass().getSimpleName() + "\t");

		return sb.toString();
	}

	/**
	 * Parse
	 * @param args
	 */
	@Override
	public void parseArgs(String[] args) {
		// Parse command line arguments
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-v")) verbose = true;
			else if ((genomeVer == null) || genomeVer.isEmpty()) genomeVer = args[i];
			else samFileNames.add(args[i]);
		}

		// Sanity check
		if ((genomeVer == null) || genomeVer.isEmpty()) usage("Missing genome version");
		if (samFileNames.size() < 1) usage("Missing SAM/BAM file/s");
	}

	public void print() {
		// Show title
		System.out.print("chr\tstart\tend\ttype\tIDs");
		for (int j = 0; j < countReadsByFile.size(); j++)
			System.out.print("\tReads:" + samFileNames.get(j) + "\tBases:" + samFileNames.get(j));
		System.out.print("\n");

		// Retrieve all possible keys, sort them
		HashSet<Marker> keys = new HashSet<Marker>();
		for (CountByKey<Marker> cbt : countReadsByFile)
			keys.addAll(cbt.keySet());

		ArrayList<Marker> keysSorted = new ArrayList<Marker>(keys.size());
		keysSorted.addAll(keys);
		Collections.sort(keysSorted);

		// Show results
		for (Marker key : keysSorted) {
			// Show 'key' information in first columns
			System.out.print(key.getChromosomeName() //
					+ "\t" + (key.getStart() + 1) //
					+ "\t" + (key.getEnd() + 1) //
					+ "\t" + idChain(key) //
			);

			// Show counter data
			for (int idx = 0; idx < countReadsByFile.size(); idx++)
				System.out.print("\t" + countReadsByFile.get(idx).get(key) + "\t" + countBasesByFile.get(idx).get(key));
			System.out.print("\n");
		}
		System.out.print("\n");
	}

	/**
	 * Run
	 */
	@Override
	public boolean run() {
		// Load Config
		if (verbose) Timer.showStdErr("Reading configuration file '" + configFile + "'");
		config = new Config(genomeVer, configFile); // Read configuration
		if (verbose) Timer.showStdErr("done");

		// Load database
		if (verbose) Timer.showStdErr("Reading database for genome '" + genomeVer + "'");
		config.loadSnpEffectPredictor(); // Read snpEffect predictor
		if (verbose) Timer.showStdErr("done");

		// Build forest
		if (verbose) Timer.showStdErr("Building interval forest");
		snpEffectPredictor = config.getSnpEffectPredictor(); // Read snpEffect predictor
		snpEffectPredictor.buildForest();
		if (verbose) Timer.showStdErr("done");

		runCountIntervals();
		print();

		return true;
	}

	/**
	 * Count reads onto intervals
	 */
	void runCountIntervals() {
		Genome genome = config.getGenome();

		// Iterate over all BAM/SAM files
		for (String samFileName : samFileNames) {
			try {
				if (verbose) Timer.showStdErr("Reading reads file '" + samFileName + "'");
				CountByKey<Marker> countReads = new CountByKey<Marker>();
				CountByKey<Marker> countBases = new CountByKey<Marker>();

				// Open file
				int readNum = 1;
				SAMFileReader sam = new SAMFileReader(new File(samFileName));
				for (SAMRecord samRecord : sam) {
					try {

						if (!samRecord.getReadUnmappedFlag()) { // Mapped?
							Chromosome chr = genome.getChromosome(samRecord.getReferenceName());
							if (chr != null) {
								// Create a marker
								Marker read = new Marker(chr, samRecord.getAlignmentStart(), samRecord.getAlignmentEnd(), 1, "");

								// Find all intersects
								Set<Marker> regions = snpEffectPredictor.regionsMarkers(read);
								for (Marker m : regions) {
									countReads.inc(m); // Count reads
									countBases.inc(m, m.intersectSize(read)); // Count number bases that intersect
								}
							}
						}

						if (verbose) Gpr.showMark(readNum, SHOW_EVERY);
						readNum++;
					} catch (Exception e) {
						e.printStackTrace();
					}

				}
				sam.close();

				// Add count to list
				countReadsByFile.add(countReads);
				countBasesByFile.add(countBases);
			} catch (Exception e) {
				e.printStackTrace();
			}

			if (verbose) {
				System.err.println("");
				Timer.showStdErr("Finished reding file " + samFileName);
			}
		}
		if (verbose) Timer.showStdErr("Done.");
	}

	@Override
	public void usage(String message) {
		if (message != null) System.err.println("Error: " + message + "\n");
		System.err.println("snpEff version " + VERSION);
		System.err.println("Usage: snpEff countReads [-v] genome readsFile_1 readsFile_2 ...  readsFile_N");
		System.err.println("\treadsFile : A file contianing the reads. Either BAM or SAM format.");
		System.exit(-1);
	}
}
