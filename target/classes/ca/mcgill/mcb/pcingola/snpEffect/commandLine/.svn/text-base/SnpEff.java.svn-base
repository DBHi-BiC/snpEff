package ca.mcgill.mcb.pcingola.snpEffect.commandLine;

import java.util.ArrayList;
import java.util.HashMap;

import ca.mcgill.mcb.pcingola.Config2DownloadTable;
import ca.mcgill.mcb.pcingola.Pcingola;
import ca.mcgill.mcb.pcingola.logStatsServer.LogStats;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.spliceSites.SpliceAnalysis;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * Command line program
 * 
 * @author pcingola
 */
public class SnpEff implements CommandLine {

	/**
	 *  Available gene database formats
	 */
	public enum GeneDatabaseFormat {
		BIOMART, GFF3, GFF2, GTF22, REFSEQ, KNOWN_GENES, GENBANK, EMBL
	}

	/**
	 *  Available input formats
	 */
	public enum InputFormat {
		TXT, PILEUP, VCF, BED
	}

	/**
	 *  Available output formats
	 */
	public enum OutputFormat {
		TXT, VCF, BED, BEDANN, GATK
	}

	public static final int COMMAND_LINE_WIDTH = 40;

	public static final String SOFTWARE_NAME = "SnpEff";
	public static final String BUILD = "2012-12-04";
	public static final String VERSION_MAJOR = "3.1";
	public static final String REVISION = "h";
	public static final String VERSION_SHORT = VERSION_MAJOR + REVISION;
	public static final String VERSION = SOFTWARE_NAME + " " + VERSION_SHORT + " (build " + BUILD + "), by " + Pcingola.BY;

	public static final String DEFAULT_SUMMARY_FILE = "snpEff_summary.html";
	public static final String DEFAULT_SUMMARY_GENES_FILE = "snpEff_genes.txt";

	protected String command = "";
	protected String[] args; // Arguments used to invoke this command
	protected String[] shiftArgs;
	protected boolean verbose = false; // Be verbose
	protected boolean quiet = false; // Be quiet
	protected boolean log = true; // Log to server (statistics)
	protected boolean multiThreaded = false; // Use multiple threads
	protected int numWorkers = Gpr.NUM_CORES; // Max number of threads (if multi-threaded version is available)
	protected int inOffset = 1; // By default positions are 1-based
	protected int outOffset = 1;
	protected String configFile; // Config file
	protected String genomeVer; // Genome version
	protected Config config; // Configuration

	/**
	 * Main
	 */
	public static void main(String[] args) {
		// Parse
		SnpEff snpEff = new SnpEff();
		snpEff.parseArgs(args);

		// Run
		boolean ok = snpEff.run();
		int retCode = ok ? 0 : -1;

		// Exit
		System.exit(retCode);
	}

	public SnpEff() {
		genomeVer = ""; // Genome version
		configFile = Config.DEFAULT_CONFIG_FILE; // Config file
	}

	/**
	 * 	Command line argument list (try to fit it into COMMAND_LINE_WIDTH)
	 * 
	 * @param splitLines
	 * @return
	 */
	String commandLineStr(boolean splitLines) {
		StringBuilder argsList = new StringBuilder();
		argsList.append("SnpEff " + command + " ");
		int size = argsList.length();

		for (String arg : args) {
			argsList.append(arg);
			size += arg.length();
			if (splitLines && (size > COMMAND_LINE_WIDTH)) {
				argsList.append(" \n");
				size = 0;
			} else {
				argsList.append(" ");
				size++;
			}
		}

		return argsList.toString();
	}

	/**
	 * Show an error (if not 'quiet' mode)
	 * @param message
	 */
	public void error(Throwable e, String message) {
		if (verbose && (e != null)) e.printStackTrace();
		if (!quiet) System.err.println(message);
	}

	/**
	 * Show an error message and exit
	 * @param message
	 */
	public void fatalError(String message) {
		System.err.println(message);
		System.exit(-1);
	}

	/**
	 * Parse command line arguments
	 * @param args
	 */
	@Override
	public void parseArgs(String[] args) {
		this.args = args;
		if (args.length <= 0) usage("Missing command");

		if (args[0].equalsIgnoreCase("build") //
				|| args[0].equalsIgnoreCase("dump") //
				|| args[0].equalsIgnoreCase("cds") //
				|| args[0].equalsIgnoreCase("eff") //
				|| args[0].equalsIgnoreCase("download") //
				|| args[0].equalsIgnoreCase("protein") //
				|| args[0].equalsIgnoreCase("closestExon") //
				|| args[0].equalsIgnoreCase("test") //
				|| args[0].equalsIgnoreCase("cfg2table") //
				|| args[0].equalsIgnoreCase("spliceAnalysis") //
				|| args[0].equalsIgnoreCase("countReads") //
		) {
			command = args[0].toLowerCase();

			// Copy all args except initial 'command'
			ArrayList<String> argsList = new ArrayList<String>();
			for (int i = 1; i < args.length; i++) {
				if (args[i].equalsIgnoreCase("-noLog")) log = false; // This option is always available (to allow privacy in all commands)
				else argsList.add(args[i]);

				if (args[i].equals("-v")) verbose = true; // Make this option available here as well
				if (args[i].equals("-q")) quiet = true; // Make this option available here as well
			}
			shiftArgs = argsList.toArray(new String[0]);
		} else {
			command = "eff"; // Default command is 'eff'
			shiftArgs = args;
		}
	}

	/**
	 * Additional values to be reported
	 * @return
	 */
	public HashMap<String, String> reportValues() {
		HashMap<String, String> reportValues = new HashMap<String, String>();
		return reportValues;
	}

	/**
	 * Run according to command line options
	 */
	@Override
	public boolean run() {
		boolean ok = false;
		SnpEff snpEff = null;

		if (command.equals("build")) {
			//---
			// Build database
			//---
			snpEff = new SnpEffCmdBuild();
			snpEff.parseArgs(shiftArgs);
		} else if (command.equals("dump")) {
			//---
			// Dump database
			//---
			snpEff = new SnpEffCmdDump();
			snpEff.parseArgs(shiftArgs);
		} else if (command.equals("download")) {
			//---
			// Download database
			//---
			snpEff = new SnpEffCmdDownload();
			snpEff.parseArgs(shiftArgs);
		} else if (command.equals("cds")) {
			//---
			// CDS test
			//---
			snpEff = new SnpEffCmdCds();
			snpEff.parseArgs(shiftArgs);
		} else if (command.equals("eff")) {
			//---
			// Align to reference genome
			//---
			snpEff = new SnpEffCmdEff();
			snpEff.parseArgs(shiftArgs);
		} else if (command.equals("protein")) {
			//---
			// Protein test
			//---
			snpEff = new SnpEffCmdProtein();
			snpEff.parseArgs(shiftArgs);
		} else if (command.equals("closestexon")) {
			//---
			// Find closest exon
			//---
			snpEff = new SnpEffCmdClosestExon();
			snpEff.parseArgs(shiftArgs);
		} else if (command.equals("cfg2table")) {
			// Create download table and galaxy list from config file
			snpEff = new Config2DownloadTable();
			snpEff.parseArgs(shiftArgs);
		} else if (command.equals("test")) {
			//---
			// Test command (only for testing weird stuff)
			//---
			snpEff = new SnpEffCmdTest();
			snpEff.parseArgs(shiftArgs);
		} else if (command.equals("spliceanalysis")) {
			//---
			// Splice site analysis
			//---
			snpEff = new SpliceAnalysis();
			snpEff.parseArgs(shiftArgs);
		} else if (command.equals("countreads")) {
			//---
			// Count reads site analysis
			//---
			snpEff = new SnpEffCmdCountReads();
			snpEff.parseArgs(shiftArgs);
		} else throw new RuntimeException("Unknown command '" + command + "'");

		//---
		// Run
		//---
		String err = "";
		try {
			ok = snpEff.run();
		} catch (Throwable t) {
			err = t.getMessage();
			t.printStackTrace();
		}

		// Report to server (usage statistics) 
		if (log) {
			LogStats logStats = LogStats.report(SOFTWARE_NAME, VERSION_SHORT, VERSION, ok, verbose, args, err, snpEff.reportValues());
			if (!quiet && logStats.isNewVersion()) {
				System.err.println("New version available: " //
						+ "\n\tNew version  : " + logStats.getLatestVersion() // 
						+ "\n\tRelease date : " + logStats.getLatestReleaseDate() //
						+ "\n\tDownload URL : " + logStats.getLatestUrl() //
						+ "\n\nTo update run:\n\tjava snpEff.jar download -v snpeff\n" //
				);
			}
		}

		return ok;
	}

	/**
	 * Show 'usage' message and exit with an error code '-1'
	 * @param message
	 */
	@Override
	public void usage(String message) {
		if (message != null) System.err.println("Error: " + message + "\n");
		System.err.println("snpEff version " + VERSION);
		System.err.println("Usage: snpEff [command] [options] [files]");
		System.err.println("\nAvailable vommands: ");
		System.err.println("   [eff]           : Calculate effect of variants. Default (no command or 'eff').");
		System.err.println("   download        : Download a SnpEff database.");
		System.err.println("   build           : Build a SnpEff database.");
		System.err.println("   dump            : Dump to STDOUT a SnpEff database (mostly used for debugging).");
		System.err.println("   cds             : Compare CDS sequences calculated form a SnpEff database to the one in a FASTA file. Used for checking databases correctness.");
		System.err.println("   protein         : Compare protein sequences calculated form a SnpEff database to the one in a FASTA file. Used for checking databases correctness.");
		System.err.println("   closestExon     : Calculate closes exon/s given a set of genomic positions or intervals.");
		System.err.println("   spliceAnalysis  : Perform an analysis of splice sites. Experimental feature.");
		System.err.println("   countReads      : Count how many reads (from a BAM file) overlap with each genomic interval. Experimental feature.");
		System.err.println("\nRun 'java -jar snpEff.jar command' for help on each specifig command");
		System.exit(-1);
	}
}
