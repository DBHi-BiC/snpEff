package ca.mcgill.mcb.pcingola.coverage;

import java.io.Serializable;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import ca.mcgill.mcb.pcingola.fileIterator.SamFileIterator;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.sam.SamEntry;
import ca.mcgill.mcb.pcingola.sam.SamHeaderRecord;
import ca.mcgill.mcb.pcingola.sam.SamHeaderRecordSq;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Base by base coverage (one chromsome)
 * 
 * @author pcingola
 */
public class Coverage implements Serializable {

	private static int SHOW_EVERY = 10 * 1000;
	private static final long serialVersionUID = 1150158182247576450L;

	HashMap<String, CoverageChr> coverageByName;

	/**
	 * Calculate coverage from a SAM file
	 * @param samFile
	 * @param verbose
	 * @return
	 */
	public static Coverage calculateFromSam(String samFile, boolean verbose) {
		int i = 1;
		boolean header = true;
		String chrPrev = "";
		SamFileIterator sfi = new SamFileIterator(samFile);
		Coverage coverage = new Coverage();
		Pattern patternCigar = Pattern.compile("(\\d+)([A-Z])");

		if( verbose ) Timer.showStdErr("Processing file '" + samFile + "'");

		for( SamEntry se : sfi ) {
			if( header ) {
				header = false;
				// Create coverage records
				for( SamHeaderRecord rec : sfi.getHeaders().getRecords("SQ") ) {
					SamHeaderRecordSq sq = (SamHeaderRecordSq) rec;
					coverage.createChr(sq.getSequenceName(), sq.getLength());
				}
			}

			// Find coverage object
			String chrName = se.getRname();

			// Get start position
			int start = se.getPos() - 1; // One-based coordinates
			String cigar = se.getCigar();
			Matcher matcher = patternCigar.matcher(cigar);
			while(matcher.find()) {
				int len = Gpr.parseIntSafe(matcher.group(1));
				String op = matcher.group(2);
				if( op.equals("M") ) coverage.inc(chrName, start, start + len - 1);
				if( op.equals("M") || op.equals("D") || op.equals("N") || op.equals("EQ") || op.equals("X") || op.equals("P") ) start += len;
			}

			// Show mark
			if( verbose ) {
				if( !chrName.equals(chrPrev) ) {
					System.err.println("");
					Timer.showStdErr(chrName + "\t");
				}
				chrPrev = chrName;
				Gpr.showMark(i++, SHOW_EVERY);
			}
		}

		return coverage;
	}

	public Coverage() {
		coverageByName = new HashMap<String, CoverageChr>();
	}

	/**
	 * Calculate average coverage per base
	 * 
	 * @param m : A marker interval
	 * @return
	 */
	public double avgCoverage(Marker m) {
		String chr = m.getChromosomeName();
		CoverageChr cchr = get(chr);
		if( cchr == null ) throw new RuntimeException("Chromosome '" + chr + "' not found!");
		return cchr.avgCoverage(m.getStart(), m.getEnd());
	}

	/**
	 * Calculate total coverage per base
	 * 
	 * @param m : A marker interval
	 * @return
	 */
	public long coverage(Marker m) {
		String chr = m.getChromosomeName();
		CoverageChr cchr = get(chr);
		if( cchr == null ) throw new RuntimeException("Chromosome '" + chr + "' not found!");
		return cchr.coverage(m.getStart(), m.getEnd());
	}

	/**
	 * Create new chromosome coverage
	 * @param chr
	 * @param len
	 */
	public void createChr(String chr, int len) {
		coverageByName.put(Chromosome.simpleName(chr), new CoverageChr(len));
	}

	public CoverageChr get(String chr) {
		return coverageByName.get(Chromosome.simpleName(chr));
	}

	/**
	 * Increment a region
	 * @param start
	 * @param end
	 */
	public void inc(String chr, int start, int end) {
		// Increment a part of the chromosome
		chr = Chromosome.simpleName(chr);
		get(chr).inc(start, end);
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		return sb.toString();
	}
}
