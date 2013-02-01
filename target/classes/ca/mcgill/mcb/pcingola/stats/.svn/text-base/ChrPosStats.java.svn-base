package ca.mcgill.mcb.pcingola.stats;

import ca.mcgill.mcb.pcingola.stats.plot.GooglePlotInt;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * How many changes per position do we have in a chromosome.
 * Summary by dividing the chromosome into MAX_BINS bins
 * 
 * @author pcingola
 */
public class ChrPosStats {

	public static int MAX_BINS = 300; // Max number of points to show in a plot (plots to sample)
	String chrName; // Chromosome name
	int chrLength; // Chromosome length;
	int factor = 1;
	int count[];
	int total;

	public ChrPosStats(String chrName, int chrLength) {
		this.chrName = chrName;
		this.chrLength = chrLength;
		total = 0;

		// Multiplier factor 
		factor = 1;
		for (factor = 1; (chrLength / factor) > MAX_BINS; factor *= 10);

		int len = chrLength / factor + 1; // May be the chromosome is smaller then 'MAX_POINTS' (e.g. when you have small contigs)

		// Initialize count
		count = new int[len];
		for (int i = 0; i < count.length; i++)
			count[i] = 0;
	}

	String factorStr() {
		if (factor > 1000000000) return factor / 1000000000 + "Gb";
		if (factor > 1000000) return factor / 1000000 + "Mb";
		if (factor > 1000) return factor / 1000 + "Kb";
		return factor + "b";
	}

	public int getTotal() {
		return total;
	}

	int[] posArray() {
		int pos[] = new int[count.length];
		for (int i = 0; i < pos.length; i++)
			pos[i] = i * factor;
		return pos;
	}

	/**
	 * Use 'num' as a sample
	 * @param num
	 */
	public void sample(int position) {
		// Ignore counts for zero or one-length chromosomes
		if (chrLength <= 1) {
			Gpr.debug("Warning: Chromosome '" + chrName + "' has length " + chrLength);
			return;
		}

		int i = position / factor;
		if ((i >= 0) && (i < count.length)) {
			count[i]++;
			total++;
		} else Gpr.debug("Error counting samples on chromosome '" + chrName + "'. Position '" + position + "' => count[" + i + "]  (count.length: " + count.length + ", factor: " + factor + ", chrLength: " + chrLength + ").");
	}

	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append("Chromosome: " + chrName + "\n");
		sb.append("\tPosition :\t");

		// Position
		int pos[] = posArray();
		for (int i = 0; i < count.length; i++)
			sb.append(pos[i] + "\t");
		sb.append("\n");

		// Counts
		sb.append("\tCount    :\t");
		for (int i = 0; i < count.length; i++)
			sb.append(count[i] + "\t");
		sb.append("\n");

		return sb.toString();
	}

	/**
	 * Create a histogram plot using Google charts
	 * @return
	 */
	public String toStringHistoPlot(String title, String xAxisLabel, String yAxisLabel) {
		int pos[] = posArray(); // Create data arrays

		GooglePlotInt gghisto = new GooglePlotInt(pos, count, title, xAxisLabel, yAxisLabel + "/" + factorStr()); // Create histogram
		return gghisto.toURLString();
	}
}
