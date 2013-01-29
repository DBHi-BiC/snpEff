package ca.mcgill.mcb.pcingola.snpEffect.commandLine;

import ca.mcgill.mcb.pcingola.interval.Cds;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.interval.Utr;
import ca.mcgill.mcb.pcingola.interval.tree.IntervalTree;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Command line program: Build database
 * 
 * @author pcingola
 */
public class SnpEffCmdDump extends SnpEff {

	public enum DumpFormat {
		Simple, Bed
	}

	DumpFormat dumpFormat = DumpFormat.Simple;
	String chrStr = "";

	public SnpEffCmdDump() {
		super();
	}

	/**
	 * Parse command line arguments
	 * @param args
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
				} else if (args[i].equals("-v") || args[i].equalsIgnoreCase("-verbose")) {
					verbose = true;
					quiet = false;
				} else if (args[i].equals("-q") || args[i].equalsIgnoreCase("-quiet")) {
					quiet = true;
					verbose = false;
				} else if (args[i].equalsIgnoreCase("-chr")) chrStr = args[++i];
				else if ((args[i].equals("-if") || args[i].equalsIgnoreCase("-inOffset"))) {
					if ((i + 1) < args.length) inOffset = Gpr.parseIntSafe(args[++i]);
				} else if (args[i].equals("-1")) inOffset = outOffset = 1;
				else if (args[i].equals("-0")) inOffset = outOffset = 0;
				else if (args[i].equals("-bed")) {
					dumpFormat = DumpFormat.Bed;
					inOffset = outOffset = 0;
				} else if (args[i].equals("-h") || args[i].equalsIgnoreCase("-help")) {
					usage(null);
					System.exit(0);
				} else usage("Unknow option '" + args[i] + "'");
			} else if (genomeVer.length() <= 0) genomeVer = args[i];
			else usage("Unknow parameter '" + args[i] + "'");
		}

		// Check: Do we have all required parameters?
		if (genomeVer.isEmpty()) usage("Missing genomer_version parameter");
	}

	/**
	 * Show all intervals in BED format
	 * References: http://genome.ucsc.edu/FAQ/FAQformat.html#format1
	 */
	void printBed() {
		for (IntervalTree tree : config.getSnpEffectPredictor().getIntervalForest()) {
			for (Marker i : tree) {
				printBed(i);

				// Show gene specifics
				if (i instanceof Gene) {
					Gene g = (Gene) i;

					// Show transcripts: UTR and Exons
					for (Transcript t : g) {
						printBed(t);

						for (Cds c : t.getCds())
							printBed(c);

						for (Utr u : t.getUtrs())
							printBed(u);

						for (Exon e : t)
							printBed(e);
					}
				}
			}
		}
	}

	/**
	 * Show a marker in BED format
	 * @param marker
	 */
	void printBed(Marker marker) {
		String chr = chrStr + marker.getChromosome().getId();
		int start = marker.getStart() + outOffset; // The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0. 
		int end = marker.getEnd() + outOffset + 1; // The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature.
		String name = marker.getClass().getSimpleName() + "_" + marker.getId();
		System.out.println(chr + "\t" + start + "\t" + end + "\t" + name);
	}

	/**
	 * Run according to command line options
	 */
	@Override
	public boolean run() {
		//---
		// Dump database
		//---
		config = new Config(genomeVer, configFile);

		// Read database
		if (verbose) Timer.showStdErr("Reading database for genome '" + genomeVer + "' (this might take a while)");
		config.loadSnpEffectPredictor(); // Read snpEffect predictor
		if (verbose) Timer.showStdErr("done");

		// Build forest
		if (verbose) Timer.showStdErr("Building interval forest");
		config.getSnpEffectPredictor().buildForest();
		if (verbose) Timer.showStdErr("Done.");

		// Dump database
		if (dumpFormat == DumpFormat.Simple) config.getSnpEffectPredictor().print();
		else if (dumpFormat == DumpFormat.Bed) printBed();
		else throw new RuntimeException("Unimplemented format '" + dumpFormat + "'");

		return true;
	}

	/**
	 * Show 'usage;' message and exit with an error code '-1'
	 * @param message
	 */
	@Override
	public void usage(String message) {
		if (message != null) System.err.println("Error: " + message + "\n");
		System.err.println("snpEff version " + VERSION);
		System.err.println("Usage: snpEff dump [options] genome_version");
		System.err.println("\t-bed                    : Dump in BED format (implies -0)");
		System.err.println("\t-chr <string>           : Prepend 'string' to chromosome name (e.g. 'chr1' instead of '1')");
		System.err.println("\nGeneric options:");
		System.err.println("\t-0                      : File positions are zero-based");
		System.err.println("\t-1                      : File positions are one-based");
		System.err.println("\t-c , -config            : Specify config file");
		System.err.println("\t-h , -help              : Show this help and exit");
		System.err.println("\t-noLog                  : Do not report usage statistics to server");
		System.err.println("\t-q , -quiet             : Quiet mode (do not show any messages or errors)");
		System.err.println("\t-v , -verbose           : Verbose mode");
		System.exit(-1);
	}
}
