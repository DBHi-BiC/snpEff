package ca.mcgill.mcb.pcingola.snpEffect.commandLine;

import java.util.HashMap;

import ca.mcgill.mcb.pcingola.codons.CodonTable;
import ca.mcgill.mcb.pcingola.codons.CodonTables;
import ca.mcgill.mcb.pcingola.fileIterator.FastaFileIterator;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Command line: Calculate coding sequences from a file and compare them to the ones calculated from our data structures
 * 
 * Note: This is done in order to see potential incompatibility 
 *       errors between genome sequence and annotation.
 * 
 * @author pcingola
 */
public class SnpEffCmdCds extends SnpEff {

	public static boolean debug = false;
	public static boolean onlyOneError = false; // This is used in some test-cases
	public static double maxErrorPercentage = 0.01; // Maximum allowed error is 1% (otherwise test fails)

	//	boolean verbose = false;
	int totalErrors = 0;
	int totalOk = 0;
	int totalWarnings = 0;
	int totalNotFound = 0;
	String configFile = Config.DEFAULT_CONFIG_FILE;
	String genomeVer = "";
	String cdsFile = "";
	Config config;
	HashMap<String, String> cdsByTrId;

	public SnpEffCmdCds() {
	}

	public SnpEffCmdCds(Config config) {
		this.config = config;
		cdsFile = config.getFileNameCds();
	}

	public SnpEffCmdCds(String genomeVer, String configFile, String cdsFile) {
		this.configFile = configFile;
		this.genomeVer = genomeVer;
		this.cdsFile = cdsFile;
	}

	/**
	 * Compare all CDS
	 */
	double cdsCompare() {
		int i = 1;

		if (!verbose) {
			// Show labels
			System.err.println("\tLabels:");
			System.err.println("\t\t'+' : OK");
			System.err.println("\t\t'.' : Missing");
			System.err.println("\t\t'*' : Error");
			System.err.print("\t");
		}

		// Compare all genes
		for (Gene gint : config.getGenome().getGenes())
			for (Transcript tint : gint) {
				String cds = tint.cds().toUpperCase();
				String mRna = tint.mRna().toUpperCase();
				String cdsReference = cdsByTrId.get(tint.getId());

				if (cdsReference != null) cdsReference = cdsReference.toUpperCase();

				if (cdsReference == null) {
					if (verbose) System.err.println("\nWARNING:Cannot find reference CDS for transcript '" + tint.getId() + "'");
					else System.out.print('.');
					totalNotFound++;
				} else if (cds.isEmpty()) {
					if (verbose) System.err.println("\nWARNING:Empty CDS for transcript '" + tint.getId() + "'");
					else System.out.print('.');
					totalNotFound++;
				} else if (cds.equals(cdsReference)) {
					totalOk++;
					if (!verbose) System.out.print('+');

					// Sanity check: Start and stop codons
					if ((cds != null) && (cds.length() >= 3)) {
						CodonTable ctable = CodonTables.getInstance().getTable(config.getGenome(), tint.getChromosomeName());

						// Check start codon
						String startCodon = cds.substring(0, 3);
						if (!ctable.isStart(startCodon)) {
							if (verbose) System.err.println("\nWARNING: CDS for transcript '" + tint.getId() + "' does not start with a start codon:\t" + startCodon + "\t" + cds);
							totalWarnings++;
						}

						// Check stop codon
						String stopCodon = cds.substring(cds.length() - 3, cds.length());
						if (!ctable.isStop(stopCodon)) {
							if (verbose) System.err.println("\nWARNING: CDS for transcript '" + tint.getId() + "' does not end with a stop codon:\t" + stopCodon + "\t" + cds);
							totalWarnings++;
						}
					}
				} else if (mRna.equals(cdsReference)) { // May be the file has mRNA instead of CDS?
					totalOk++;
					if (!verbose) System.out.print('+');
				} else if ((mRna.length() < cdsReference.length()) // CDS longer than mRNA? May be it is actually an mRNA + poly-A tail (instead of a CDS)
						&& cdsReference.substring(mRna.length()).replace('A', ' ').trim().isEmpty() // May be it is an mRNA and it has a ploy-A tail added
						&& cdsReference.substring(0, mRna.length()).equals(mRna) // Compare cutting poly-A tail
				) {
					// OK, it was a mRNA +  polyA
					totalOk++;
					if (!verbose) System.out.print('+');
				} else if ((mRna.length() > cdsReference.length()) // PolyA in the reference? 
						&& mRna.substring(cdsReference.length()).replace('A', ' ').trim().isEmpty() // 
						&& mRna.substring(0, cdsReference.length()).equals(mRna) // 
				) {
					// OK, it was a mRNA +  polyA
					totalOk++;
					if (!verbose) System.out.print('+');
				} else {
					if (verbose || onlyOneError) {
						// Create a string indicating differences
						String diffMrna = SnpEffCmdProtein.diffStr(mRna, cdsReference);
						int diffMrnaCount = SnpEffCmdProtein.diffCount(mRna, cdsReference);

						String diffCds = SnpEffCmdProtein.diffStr(cds, cdsReference);
						int diffCdsCount = SnpEffCmdProtein.diffCount(cds, cdsReference);

						System.err.println("\nERROR:CDS do not match for transcript " + tint.getId() + "\tStrand:" + tint.getStrand() + "\tExons: " + tint.numChilds());

						if (diffMrnaCount < diffCdsCount) {
							System.err.println(String.format("\tsnpEff mRNA (%6d) : '%s'", mRna.length(), mRna.toLowerCase()));
							System.err.println(String.format("\tdiff        (%6d) : '%s'", diffMrnaCount, diffMrna));
						} else {
							System.err.println(String.format("\tsnpEff CDS  (%6d) : '%s'", cds.length(), cds.toLowerCase()));
							System.err.println(String.format("\tdiff        (%6d) : '%s'", diffCdsCount, diffCds));
						}

						System.err.println(String.format("\tReference   (%6d) : '%s'", cdsReference.length(), cdsReference.toLowerCase()));
						System.err.println("Transcript details:\n" + tint);
					} else System.out.print('*');

					totalErrors++;

					if (onlyOneError) {
						System.err.println("Transcript details:\n" + tint);
						throw new RuntimeException("DIE");
					}
				}

				// Show a mark
				if (!verbose && (i % 100 == 0)) System.out.print("\n\t");
				i++;
			}

		double perc = ((double) totalErrors) / ((double) (totalErrors + totalOk));
		System.out.println("\nCDS check:\t" + config.getGenome().getVersion() + "\tOK: " + totalOk + "\tWarnings: " + totalWarnings + "\tNot found: " + totalNotFound + "\tErrors: " + totalErrors + "\tError percentage: " + (100 * perc) + "%");
		return perc;
	}

	/**
	 * Parse command line arguments
	 */
	@Override
	public void parseArgs(String[] args) {
		for (int i = 0; i < args.length; i++) {

			// Argument starts with '-'?
			if (args[i].startsWith("-")) {
				if ((args[i].equals("-c") || args[i].equalsIgnoreCase("-config"))) {
					if ((i + 1) < args.length) configFile = args[++i];
					else usage("Option '-c' without config file argument");
				} else if (args[i].equals("-v") || args[i].equalsIgnoreCase("-verbose")) {
					verbose = true;
				} else usage("Unknow option '" + args[i] + "'");
			} else if (genomeVer.isEmpty()) genomeVer = args[i];
			else if (cdsFile.isEmpty()) cdsFile = args[i];
			else usage("Unknow parameter '" + args[i] + "'");
		}

		// Check: Do we have all required parameters?
		if (genomeVer.isEmpty()) usage("Missing genomer_version parameter");
		if (cdsFile.isEmpty()) usage("Missing cds_file parameter");
	}

	/**
	 * Read a file that has all CDS
	 */
	void readCdsFile() {
		cdsByTrId = new HashMap<String, String>();

		if (cdsFile.endsWith("txt") || cdsFile.endsWith("txt.gz")) readCdsFileTxt();
		else readCdsFileFasta();
	}

	/**
	 * Read CDSs from a file
	 * Format: Tab-separated format, containing "sequence \t transcriptId"
	 */
	void readCdsFileFasta() {
		// Load file
		FastaFileIterator ffi = new FastaFileIterator(cdsFile);
		for (String seq : ffi) {
			String trId = ffi.getName();

			// Repeated transcript Id? => Check that CDS is the same 
			if ((cdsByTrId.get(trId) != null) && (!cdsByTrId.get(trId).equals(seq))) System.err.println("ERROR: Different CDS for the same transcript ID. This should never happen!!!\n\tLine number: " + ffi.getLineNum() + "\n\tTranscript ID:\t" + trId + "\n\tCDS:\t\t" + cdsByTrId.get(trId) + "\n\tCDS (new):\t" + seq);

			cdsByTrId.put(trId, seq); // Add it to the hash
			if (debug) Gpr.debug("Adding cdsByTrId{'" + trId + "'} :\t" + seq);
		}
	}

	/**
	 * Read CDSs from a file
	 * Format: Tab-separated format, containing "sequence \t transcriptId"
	 */
	void readCdsFileTxt() {
		// Load file
		String cdsData = Gpr.readFile(cdsFile);
		String cdsLines[] = cdsData.split("\n");

		// Parse each line
		int lineNum = 1;
		for (String cdsLine : cdsLines) {
			// Split tab separated fields
			String field[] = cdsLine.split("\\s+");

			// Parse fields
			if (field.length >= 2) {
				// OK Parse fields
				String seq = field[1].trim();
				String trId = field[0].trim();

				// Repeated transcript Id? => Check that CDS is the same 
				if ((cdsByTrId.get(trId) != null) && (!cdsByTrId.get(trId).equals(seq))) System.err.println("ERROR: Different CDS for the same transcript ID. This should never happen!!!\n\tLine number: " + lineNum + "\n\tTranscript ID:\t" + trId + "\n\tCDS:\t\t" + cdsByTrId.get(trId) + "\n\tCDS (new):\t" + seq);

				cdsByTrId.put(trId, seq); // Add it to the hash
			}

			lineNum++;
		}
	}

	/**
	 * Run command
	 */
	@Override
	public boolean run() {
		if (verbose) Timer.showStdErr("Checking database using CDS sequences");

		// Load config
		if (config == null) {
			if (verbose) Timer.showStdErr("Reading configuration...");
			config = new Config(genomeVer, configFile); // Read configuration
			if (verbose) Timer.showStdErr("done");
		}

		// Read CDS form file
		if (verbose) Timer.showStdErr("Reading CDSs from file '" + cdsFile + "'...");
		readCdsFile(); // Load CDS
		if (verbose) Timer.showStdErr("done (" + cdsByTrId.size() + " CDSs).");

		// Load predictor
		if (config.getSnpEffectPredictor() == null) {
			if (verbose) Timer.showStdErr("Reading database...");
			config.loadSnpEffectPredictor(); // Read snpEffect predictor
			if (verbose) Timer.showStdErr("done");
		}

		// Compare CDS
		if (verbose) Timer.showStdErr("Comparing CDS...");
		cdsCompare();
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
		System.err.println("Usage: snpEff cds [options] genome_version cds_file");
		System.err.println("\nOptions:");
		System.err.println("\t-c , -config            : Specify config file");
		System.err.println("\t-noLog                  : Do not report usage statistics to server");
		System.err.println("\t-v , -verbose           : Verbose mode");
		System.exit(-1);
	}
}
