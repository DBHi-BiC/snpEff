package ca.mcgill.mcb.pcingola.snpEffect.commandLine;

import ca.mcgill.mcb.pcingola.fileIterator.LineFileIterator;
import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;

/**
 * Command line: Test
 * 
 * Note: Used for testing weird stuff
 * 
 * @author pcingola
 */
public class SnpEffCmdTest extends SnpEff {

	String vcfFile;
	int testNum;

	public SnpEffCmdTest() {
	}

	/**
	 * Parse command line arguments
	 */
	@Override
	public void parseArgs(String[] args) {
		if (args.length != 2) usage("Missing file.vcf");
		testNum = Gpr.parseIntSafe(args[0]);
		vcfFile = args[1];
	}

	/**
	 * Run command
	 */
	@Override
	public boolean run() {
		Timer.showStdErr("Test " + testNum + " : Start");
		String testType = "";
		VcfFileIterator vcfFileIterator;

		switch (testNum) {
		case 1:
			testType = "VcfFileIterator";
			vcfFileIterator = new VcfFileIterator(vcfFile);
			for (VcfEntry vcfEntry : vcfFileIterator)
				System.out.println(vcfEntry);
			break;

		case 2:
			testType = "LineFileIterator";
			LineFileIterator lfi = new LineFileIterator(vcfFile);
			for (String line : lfi)
				System.out.println(line);

			break;

		case 3:
			testType = "VcfFileIterator, no output";
			vcfFileIterator = new VcfFileIterator(vcfFile);
			long s = 0;
			for (VcfEntry vcfEntry : vcfFileIterator)
				s += vcfEntry.getStart();
			Gpr.debug("Meaningless calculation: " + s);
			break;

		case 4:
			testType = "Print, no read";
			vcfFileIterator = new VcfFileIterator(vcfFile);
			VcfEntry vcfEntry = vcfFileIterator.next(); // Read one entry
			// Print it millons of times
			for (int i = 0; i < 1000000; i++)
				System.out.println(vcfEntry);
			vcfFileIterator.close();
			break;

		default:
			usage("Unknown test " + testNum);
		}

		Timer.showStdErr("Test " + testType + " : End");
		return true;
	}

	/**
	 * Show usage and exit
	 */
	@Override
	public void usage(String message) {
		if (message != null) System.err.println("Error: " + message + "\n");
		System.err.println("snpEff version " + SnpEff.VERSION);
		System.err.println("Usage: snpEff test testNumber file.vcf");
		System.exit(-1);
	}
}
