package ca.mcgill.mcb.pcingola.snpEffect.testCases;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import junit.framework.Assert;
import junit.framework.TestCase;
import ca.mcgill.mcb.pcingola.binseq.DnaSequence;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.MarkerUtil;
import ca.mcgill.mcb.pcingola.interval.Markers;
import ca.mcgill.mcb.pcingola.interval.Variant;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.interval.tree.IntervalForest;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.GprSeq;

public class TestCasesIntervals extends TestCase {

	int maxLen = 100;
	boolean verbose = true;
	boolean compareCdsTestsEnable = false;
	Random rand;
	Genome genome;

	// Markers intervals;

	// Config config;

	public TestCasesIntervals() {
		super();
		initRand();
		String genomeName = "testCaseHg";
		genome = new Genome(genomeName);

		// Create chromosomes
		for (int i = 1; i < 22; i++)
			genome.add(new Chromosome(genome, 0, 0, 1, "" + i));
		genome.add(new Chromosome(genome, 0, 0, 1, "X"));
		genome.add(new Chromosome(genome, 0, 0, 1, "Y"));
	}

	/** 
	 * Compare each result. If one matches, we consider it OK
	 * @param transcriptId
	 * @param seqChange
	 * @param resultsList
	 * @param useSimple
	 * @param resultsSoFar
	 * @return
	 */
	boolean anyResultMatches(String transcriptId, Variant seqChange, List<ChangeEffect> resultsList, boolean useSimple, StringBuilder resultsSoFar) {
		boolean ok = false;
		for (ChangeEffect chEff : resultsList) {
			String resStr = chEff.toStringSimple(!useSimple);
			if (verbose) System.out.println(seqChange + "\t'" + resStr + "'");

			String effTrId = chEff.getMarker().findParent(Transcript.class).getId();
			if ((transcriptId == null) || (transcriptId.equals(effTrId))) {

				if (!seqChange.getId().equals(resStr)) {
					// SNP effect does not match this result  
					if (verbose) Gpr.debug("SeqChange: " + seqChange + "\tResult: '" + chEff + "'");
					resultsSoFar.append(seqChange + "\t'" + resStr + "'\n");
				} else {
					// SNP effect matches one result
					ok = true;
					break;
				}
			}
		}
		return ok;
	}

	/**
	 * Compare a string and a file. Throw an exception if they don't match
	 * @param result
	 * @param fileName
	 */
	public void compareToFile(String result, String fileName) {
		String file = Gpr.readFile(fileName).trim();
		result = result.trim();
		if (!result.equals(file)) {
			String errMsg = "Results and file '" + fileName + "' do not match:\n--- Results ---\n" + result + "\n--- File '" + fileName + "' ---\n" + file + "\n--- End ---\n";
			System.err.println(errMsg);
			throw new RuntimeException(errMsg);
		}
	}

	/**
	 * Tests 'baseAt' method in DnaSequence
	 */
	public void DnaSequenceBaseAt(int len) {
		// Create a random sequence
		char bases[] = new char[len];
		for (int i = 0; i < bases.length; i++) {
			char base = GprSeq.BASES[(int) (Math.random() * 4)];
			bases[i] = base;
		}

		String sequence = new String(bases);
		DnaSequence DnaSequence = new DnaSequence(sequence);
		System.out.println("DnaSequence (len:" + len + ") : " + DnaSequence);

		for (int i = 0; i < bases.length; i++) {
			char base = Character.toUpperCase(DnaSequence.getBase(i));
			if (base != bases[i]) throw new RuntimeException("Bases do not match! Base:" + base + "\tOriginal sequence: " + bases[i]);
		}
	}

	void initRand() {
		rand = new Random(20100629);
	}

	public Markers intersects(Markers interval1, Marker intervals) {
		Markers ints = new Markers();
		for (Marker i : interval1)
			if (i.intersects(intervals)) ints.add(i);
		return ints;
	}

	/**
	 * Create a random set of intervals
	 * @param numIntervals
	 * @param maxStart
	 * @param maxLength
	 * @return
	 */
	public Markers randomIntervals(int numIntervals, int maxStart, int maxLength, int numChromo) {
		Markers ints = new Markers();

		for (int ch = 1; ch <= numChromo; ch++) {
			for (int i = 0; i < numIntervals; i++) {
				int start = rand.nextInt(maxStart);
				int end = Math.min(start + rand.nextInt(maxLength), maxStart - 1);
				Marker interval = new Marker(genome.getChromosome("" + ch), start, end, 1, "");
				ints.add(interval);
			}
		}

		return ints;
	}

	public void test_00() {
		initRand();

		for (int len = 1; len < 1000; len++)
			DnaSequenceBaseAt(len);
	}

	/**
	 * Read file
	 */
	public void test_01() {
		initRand();
		Markers intervals = MarkerUtil.readTxt("tests/interval_data_100.txt", genome, 0);
		compareToFile(intervals.toStringTxt(), "tests/test_01.txt");
	}

	/**
	 * Sort test
	 */
	public void test_02() {
		initRand();
		Markers intervals = randomIntervals(10, maxLen, 10, 5);
		intervals.sort(false, false);
		compareToFile(intervals.toStringTxt(), "tests/test_02.txt");
	}

	/**
	 * Sort (by end) test
	 */
	public void test_03() {
		initRand();
		Markers intervals = randomIntervals(100, maxLen, 25, 2);
		intervals.sort(true, false);
		compareToFile(intervals.toStringTxt(), "tests/test_03.txt");
	}

	/**
	 * Merge intervals
	 */
	public void test_04() {
		initRand();
		Markers intervals = randomIntervals(20, maxLen, 10, 2);
		Markers merge = intervals.merge();
		if (verbose) System.out.println("Merge :\n" + merge.toStringAsciiArt(maxLen));
		compareToFile(merge.toStringTxt(), "tests/test_04.txt");
	}

	/**
	 * Adding intervals
	 */
	public void test_05() {
		initRand();
		// Create and perform union
		Markers intervals = randomIntervals(5, maxLen, 10, 2);
		Markers intervals2 = randomIntervals(5, maxLen, 10, 2);

		Markers add = new Markers();
		add.add(intervals);
		add.add(intervals2);

		compareToFile(add.toStringTxt(), "tests/test_05.txt");

		if (verbose) {
			// Sort
			intervals.sort(false, false);
			intervals2.sort(false, false);
			add.sort(false, false);

			// Show
			System.out.println("Intervals 1:\n" + intervals.toStringAsciiArt(maxLen));
			System.out.println("Intervals 2:\n" + intervals2.toStringAsciiArt(maxLen));
			System.out.println("Union :\n" + add.toStringAsciiArt(maxLen));
		}
	}

	/**
	 * Intersect of 2 intervals
	 * We have 2 implementation (brute force and interval trees), so we can compare them.
	 */
	public void test_06() {

		for (int numInts = 10; numInts < 1000; numInts *= 2) {
			// Create 'original' intervals
			Markers intervals = randomIntervals(numInts, maxLen, 20, 2);

			// Create forest (one tree per chromosome)
			IntervalForest forest = new IntervalForest(intervals);

			if (verbose) {
				intervals.sort(false, false);
				System.out.println(intervals.toStringAsciiArt(maxLen));
				forest.build();
				System.out.println(forest);
			}

			Markers intervals2 = randomIntervals(numInts, maxLen, 20, 2);

			// We test one by one in order to compare individual results
			for (Marker i : intervals2) {
				Markers intersect = forest.query(i);
				Markers intersectBf = intersects(intervals, i);

				// Show
				if (!intersect.equals(intersectBf)) {
					intersect.sort(false, false);
					intersectBf.sort(false, false);
					String errMsg = "Interval: " + i + "\n\tIntersects  : " + intersect + "\n\tIntersectsBf: " + intersectBf;
					System.err.println(errMsg);
					throw new RuntimeException(errMsg);
				}
			}
		}
	}

	/**
	 * Minus operation for intervals
	 */
	public void test_07_01() {
		Chromosome chr = genome.getChromosome("1");

		// Two identical intervals => Result should be empty
		Markers intervals = new Markers();
		intervals.add(new Marker(chr, 10, 90, 1, ""));

		Markers intervals2 = new Markers();
		intervals2.add(new Marker(chr, 10, 90, 1, ""));

		Markers minus = intervals.minus(intervals2);
		assertEquals(0, minus.size());

		if (verbose) System.out.println(minus.toStringAsciiArt(maxLen));
	}

	/**
	 * Minus operation for intervals
	 */
	public void test_07_02() {
		Chromosome chr = genome.getChromosome("1");

		// Totally overlapping => result should be empty 
		Markers intervals = new Markers();
		intervals.add(new Marker(chr, 10, 90, 1, ""));

		Markers intervals2 = new Markers();
		intervals2.add(new Marker(chr, 9, 91, 1, ""));

		Markers minus = intervals.minus(intervals2);
		assertEquals(0, minus.size());

		if (verbose) System.out.println(minus.toStringAsciiArt(maxLen));
	}

	/**
	 * Minus operation for intervals
	 */
	public void test_07_03() {
		Chromosome chr = genome.getChromosome("1");

		// Overlapping right part => result should be left part 
		Markers intervals = new Markers();
		intervals.add(new Marker(chr, 10, 90, 1, ""));

		Markers intervals2 = new Markers();
		intervals2.add(new Marker(chr, 9, 50, 1, ""));

		Markers minus = intervals.minus(intervals2);
		assertEquals(1, minus.size());
		Marker minusInt = minus.iterator().next();
		assertEquals(51, minusInt.getStart());
		assertEquals(90, minusInt.getEnd());

		if (verbose) System.out.println(minus.toStringAsciiArt(maxLen));
	}

	/**
	 * Minus operation for intervals
	 */
	public void test_07_04() {
		Chromosome chr = genome.getChromosome("1");

		// Overlapping left part => result should be right part 
		Markers intervals = new Markers();
		intervals.add(new Marker(chr, 10, 90, 1, ""));

		Markers intervals2 = new Markers();
		intervals2.add(new Marker(chr, 51, 91, 1, ""));

		Markers minus = intervals.minus(intervals2);

		if (verbose) System.out.println(minus.toStringAsciiArt(maxLen));

		assertEquals(1, minus.size());
		Marker minusInt = minus.iterator().next();
		assertEquals(10, minusInt.getStart());
		assertEquals(50, minusInt.getEnd());
	}

	/**
	 * Minus operation for intervals
	 */
	public void test_07_05() {
		Chromosome chr = genome.getChromosome("1");

		// Overlapping middle => result should be left & right parts 
		Markers intervals = new Markers();
		intervals.add(new Marker(chr, 10, 90, 1, ""));

		Markers intervals2 = new Markers();
		intervals2.add(new Marker(chr, 40, 60, 1, ""));

		Markers minus = intervals.minus(intervals2);

		if (verbose) System.out.println(minus.toStringAsciiArt(maxLen));

		assertEquals(2, minus.size());
		Iterator<Marker> it = minus.iterator();
		Marker minusInt = it.next();
		assertEquals(10, minusInt.getStart());
		assertEquals(39, minusInt.getEnd());

		minusInt = it.next();
		assertEquals(61, minusInt.getStart());
		assertEquals(90, minusInt.getEnd());
	}

	/**
	 * Test distance (in bases) from a list of markers
	 */
	public void test_08() {
		Chromosome chr = genome.getChromosome("1");
		Marker m1 = new Marker(chr, 0, 100, 1, "");

		ArrayList<Marker> list = new ArrayList<Marker>();
		list.add(m1);

		int last = m1.getEnd() + 10;
		for (int i = m1.getStart(); i <= last; i++) {
			Marker m = new Marker(chr, i, i, 1, "");

			int dist = m.distanceBases(list, false);
			Assert.assertEquals(i, dist);
		}
	}

	/**
	 * Test distance (in bases) from a list of markers
	 */
	public void test_08_02() {
		Chromosome chr = genome.getChromosome("1");
		Marker m1 = new Marker(chr, 0, 99, 1, "");
		Marker m2 = new Marker(chr, 200, 299, 1, "");
		Marker m3 = new Marker(chr, 400, 499, 1, "");

		ArrayList<Marker> list = new ArrayList<Marker>();
		list.add(m1);
		list.add(m2);
		list.add(m3);

		int last = m3.getEnd() + 10;
		for (int i = m1.getStart(); i <= last; i++) {
			Marker m = new Marker(chr, i, i, 1, "");

			int dist = m.distanceBases(list, false);
			Assert.assertEquals(i % 100, dist % 100);
		}
	}

	/**
	 * Test distance (in bases) from a list of markers
	 */
	public void test_08_03() {
		Chromosome chr = genome.getChromosome("1");
		Marker m1 = new Marker(chr, 0, 99, 1, "");
		Marker m2 = new Marker(chr, 200, 299, 1, "");
		Marker m3 = new Marker(chr, 400, 499, 1, "");

		ArrayList<Marker> list = new ArrayList<Marker>();
		list.add(m1);
		list.add(m2);
		list.add(m3);

		int last = m3.getEnd() + 10;
		for (int i = m1.getStart(); i <= last; i++) {
			Marker m = new Marker(chr, i, i, 1, "");

			int dist = m.distanceBases(list, true);
			Assert.assertEquals((m3.getEnd() - i) % 100, dist % 100);
		}
	}

	public void test_09_chrOrder() {
		Genome genome = new Genome("test");

		Chromosome chrA = new Chromosome(genome, 0, 1, 1, "6");
		Chromosome chrB = new Chromosome(genome, 0, 1, 1, "10");
		Chromosome chrC = new Chromosome(genome, 0, 1, 1, "4_ctg9_hap1");

		// Order: A < B < C
		Assert.assertTrue(chrA.compareTo(chrB) < 0);
		Assert.assertTrue(chrA.compareTo(chrC) < 0);
		Assert.assertTrue(chrB.compareTo(chrC) < 0);
	}

	public void test_10_chrOrder() {
		Genome genome = new Genome("test");

		Chromosome chrA = new Chromosome(genome, 0, 1, 1, "chr1");
		Chromosome chrB = new Chromosome(genome, 0, 1, 1, "scaffold0001");

		// Order: A < B
		Assert.assertTrue(chrA.compareTo(chrB) < 0);
	}

}
