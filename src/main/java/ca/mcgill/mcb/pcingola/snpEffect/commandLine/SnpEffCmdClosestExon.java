package ca.mcgill.mcb.pcingola.snpEffect.commandLine;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;

import ca.mcgill.mcb.pcingola.fileIterator.SeqChangeBedFileIterator;
import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Markers;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.interval.tree.IntervalForest;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;

/**
 * Command line: Find closes exon to each variant
 * 
 * Note: Transcripts are ordered by 'mRNA' length in order 
 *       for exon form longer transcripts to appear first 
 *       (same exon form shorter transcripts is omitted).  
 *       
 * @author pcingola
 */
public class SnpEffCmdClosestExon extends SnpEff {

	public static final String CLOSEST_EXON = "CLOSEST_EXON";
	public static final String INFO_LINE = "##INFO=<ID=" + CLOSEST_EXON + ",Number=4,Type=String,Description=\"Closest exon: Distance (bases), exons Id, transcript Id, gene name\">";

	boolean verbose = false;
	boolean bedFormat = false;
	String configFile = Config.DEFAULT_CONFIG_FILE;
	String genomeVer = "";
	String inFile = "";
	Config config;
	IntervalForest intervalForest;

	public SnpEffCmdClosestExon() {
		command = "closestExon";
	}

	public SnpEffCmdClosestExon(Config config) {
		command = "closestExon";
		this.config = config;
		inFile = config.getFileNameProteins();
	}

	/**
	 * Update header
	 * @param vcf
	 */
	void addHeaderLines(VcfFileIterator vcf) {
		vcf.getVcfHeader().addLine("##SnpEffVersion=\"" + SnpEff.VERSION + "\"");
		vcf.getVcfHeader().addLine("##SnpEffCmd=\"" + commandLineStr(false) + "\"");
		vcf.getVcfHeader().addLine(INFO_LINE);
	}

	/**
	 * Iterate over VCF file, find closest exons and annotate vcf lines
	 */
	void bedIterate() {
		// Open file
		SeqChangeBedFileIterator bfi = new SeqChangeBedFileIterator(inFile, config.getGenome(), 0);
		bfi.setCreateChromos(true); // Any 'new' chromosome in the input file will be created (otherwise an error will be thrown)

		for (SeqChange bed : bfi) {
			try {
				// Find closest exon
				Exon exon = (Exon) findClosestExons(bed);

				String id = bed.getId();

				// Update INFO fields if any exon was found
				if (exon != null) {
					int dist = exon.distance(bed);
					Transcript tr = (Transcript) exon.getParent();
					Gene gene = (Gene) tr.getParent();
					id = (id.isEmpty() ? "" : bed.getId() + ";") + dist + "," + exon.getId() + "," + tr.getId() + "," + gene.getGeneName();
				}

				// Show output
				System.out.println(bed.getChromosomeName() //
						+ "\t" + bed.getStart() // BED format: Zero-based position
						+ "\t" + (bed.getEnd() + 1) // BED format: End base is not included
						+ "\t" + id //
				);

			} catch (Exception e) {
				e.printStackTrace(); // Show exception and move on...
			}
		}
	}

	/**
	 * Create an interval forest containing all exons.
	 * 
	 * Note: Transcripts are ordered by 'mRNA' length in order 
	 *       for exon form longer transcripts to appear first 
	 *       (same exon form shorter transcripts is omitted).  
	 * 
	 * @return
	 */
	IntervalForest createForest() {
		if (verbose) Timer.showStdErr("Creating interval forest...");

		// Build an 'exon' forest. Forget about everything else.
		intervalForest = new IntervalForest();
		HashSet<String> exons = new HashSet<String>();

		// Compare by mRNA length
		Comparator<Transcript> mRnaLenComp = new Comparator<Transcript>() {

			@Override
			public int compare(Transcript t1, Transcript t2) {
				return t2.mRna().length() - t1.mRna().length();
			}
		};

		// For all genes, find transcripts and add exons to forest
		int countAdded = 0, countSkipped = 0;
		for (Gene gene : config.getGenome().getGenes()) {

			// Sort transcripts by length
			ArrayList<Transcript> transcripts = new ArrayList<Transcript>();
			for (Transcript tr : gene)
				transcripts.add(tr);

			// Sort by mRna length
			Collections.sort(transcripts, mRnaLenComp);

			// Add exons form transcripts (longer transcripts first) 
			for (Transcript tr : transcripts) {
				for (Exon exon : tr) {
					// Create a key
					String key = exon.getChromosomeName() + ":" + exon.getStart() + "-" + exon.getEnd();

					// Already in the set? => Don't add exon
					if (!exons.contains(key)) {
						intervalForest.add(exon); // Add exon and key
						exons.add(key);
						countAdded++;
					} else countSkipped++;
				}
			}
		}

		if (verbose) Timer.showStdErr("Done. Added " + countAdded + " exons. Skipped " + countSkipped + " (redundant exons).");
		return intervalForest;
	}

	/**
	 * Find closest exon for this interval
	 * @param inputInterval
	 */
	Marker findClosestExons(Marker inputInterval) {
		int initialExtension = 1000;

		Chromosome chr = inputInterval.getChromosome();
		if (chr.size() > 0) {
			// Extend interval to capture 'close' exons
			for (int extend = initialExtension; extend < chr.size(); extend *= 2) {
				int start = Math.max(inputInterval.getStart() - extend, 0);
				int end = inputInterval.getEnd() + extend;
				Marker extended = new Marker(chr, start, end, 1, "");

				// Find all exons that intersect with the interval
				Markers markers = intervalForest.query(extended);
				int minDist = Integer.MAX_VALUE;
				Marker minDistMarker = null;
				for (Marker m : markers) {
					int dist = m.distance(inputInterval);
					if (dist < minDist) {
						minDistMarker = m;
						minDist = dist;
					}

					// Zero distance? Cannot be lower than this => return
					if (minDist <= 0) return minDistMarker;
				}

				// Found something?
				if (minDistMarker != null) return minDistMarker;
			}
		} else Gpr.debug("CHROMOSOME LENGTH IS " + chr.size());

		// Nothing found
		return null;
	}

	/**
	 * Parse command line arguments
	 */
	@Override
	public void parseArgs(String[] args) {
		this.args = args;
		for (int i = 0; i < args.length; i++) {

			// Argument starts with '-'?
			if (args[i].startsWith("-")) {
				if ((args[i].equals("-c") || args[i].equalsIgnoreCase("-config"))) {
					if ((i + 1) < args.length) configFile = args[++i];
					else usage("Option '-c' without config file argument");
				} else if (args[i].equals("-bed")) {
					bedFormat = true;
				} else if (args[i].equals("-v") || args[i].equalsIgnoreCase("-verbose")) {
					verbose = true;
				} else usage("Unknow option '" + args[i] + "'");
			} else if (genomeVer.isEmpty()) genomeVer = args[i];
			else if (inFile.isEmpty()) inFile = args[i];
			else usage("Unknow parameter '" + args[i] + "'");
		}

		// Check: Do we have all required parameters?
		if (genomeVer.isEmpty()) usage("Missing genomer_version parameter");
		if (inFile.isEmpty()) usage("Missing protein_file parameter");
	}

	/**
	 * Run command
	 */
	@Override
	public boolean run() {
		// Load config
		if (config == null) {
			if (verbose) Timer.showStdErr("Reading configuration...");
			config = new Config(genomeVer, configFile); // Read configuration
			if (verbose) Timer.showStdErr("done");
		}

		if (verbose) Timer.showStdErr("Loading predictor...");
		config.loadSnpEffectPredictor();
		if (verbose) Timer.showStdErr("done");

		createForest();

		if (verbose) Timer.showStdErr("Reading file '" + inFile + "'");
		if (bedFormat) bedIterate();
		else vcfIterate();
		if (verbose) Timer.showStdErr("done");

		return true;
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	/**
	 * Show usage and exit
	 */
	@Override
	public void usage(String message) {
		if (message != null) System.err.println("Error: " + message + "\n");
		System.err.println("snpEff version " + SnpEff.VERSION);
		System.err.println("Usage: snpEff closestExon [options] genome_version file.vcf");
		System.err.println("\nOptions:");
		System.err.println("\t-bed          : Input format is BED. Default: VCF");
		System.err.println("\t-c , -config  : Specify config file");
		System.err.println("\t-noLog        : Do not report usage statistics to server");
		System.err.println("\t-v , -verbose : Verbose mode");
		System.exit(-1);
	}

	/**
	 * Iterate over VCF file, find closest exons and annotate vcf lines
	 */
	void vcfIterate() {
		// Open file
		VcfFileIterator vcf = new VcfFileIterator(inFile, config.getGenome());
		vcf.setCreateChromos(true); // Any 'new' chromosome in the input file will be created (otherwise an error will be thrown)

		boolean header = true;
		for (VcfEntry ve : vcf) {
			try {
				if (header) {
					// Update and show header
					addHeaderLines(vcf);
					String headerStr = vcf.getVcfHeader().toString();
					if (!headerStr.isEmpty()) System.out.println(headerStr);
					header = false;
				}

				// Find closest exon
				Exon exon = (Exon) findClosestExons(ve);

				// Update INFO fields if any exon was found
				if (exon != null) {
					int dist = exon.distance(ve);
					Transcript tr = (Transcript) exon.getParent();
					Gene gene = (Gene) tr.getParent();
					String value = dist + "," + exon.getId() + "," + tr.getId() + "," + gene.getGeneName();
					ve.addInfo(CLOSEST_EXON, value);
				}

				// Show output
				System.out.println(ve);
			} catch (Exception e) {
				e.printStackTrace(); // Show exception and move on...
			}
		}
	}

}
