package ca.mcgill.mcb.pcingola.snpEffect.testCases;

import java.util.HashSet;
import java.util.Random;

import junit.framework.Assert;
import junit.framework.TestCase;
import ca.mcgill.mcb.pcingola.binseq.DnaNSequence;
import ca.mcgill.mcb.pcingola.binseq.coder.DnaCoder;
import ca.mcgill.mcb.pcingola.util.GprSeq;

public class TestCasesDnaNSequence extends TestCase {

	public static boolean verbose = false;

	/**
	 * Create random changes in a sequence
	 * @param overlap
	 * @param overlapChanges
	 * @return
	 */
	String change(String sequence, int numChanges, Random rand) {
		HashSet<Integer> changedPos = new HashSet<Integer>();
		char chars[] = sequence.toCharArray();

		for( int i = 0; i < numChanges; ) {
			int pos = rand.nextInt(chars.length);

			if( !changedPos.contains(pos) ) { // Already changed?
				char newBase = randBase(rand);

				if( chars[pos] != newBase ) { // Base is different?
					chars[pos] = newBase;
					changedPos.add(pos);
					i++;
				}
			}
		}

		return new String(chars);
	}

	/**
	 * Create a random base
	 * @param rand
	 * @return
	 */
	char randBase(Random rand) {
		int r = rand.nextInt(5);
		if( r < 4 ) return DnaCoder.get().toBase(r);
		return 'N';
	}

	/**
	 * Create random sequences and store them in a DnaNSequence. 
	 * Compare getting a few random bases from the original and DnaNSequence sequences.
	 * @param numTests
	 * @param lenMask
	 * @param seed
	 */
	public void randDnaSeqGetBasesTest(int numTests, int numTestsPerSeq, int lenMask, long seed) {
		Random rand = new Random(seed);

		for( int t = 0; t < numTests; t++ ) {
			String seq = "";
			int len = (rand.nextInt() & lenMask) + 10; // Randomly select sequence length
			seq = randSeq(len, rand); // Create a random sequence

			DnaNSequence bseq = new DnaNSequence(seq);

			// Retrieve numTestsPerSeq random bases from the sequence
			for( int i = 0; i < numTestsPerSeq; i++ ) {
				int randPos = rand.nextInt(len);
				int randLen = rand.nextInt(len - randPos);
				String basesOri = seq.substring(randPos, randPos + randLen);
				String basesBin = bseq.getBases(randPos, randLen);
				Assert.assertEquals(basesOri, basesBin);
				if( verbose ) System.out.println("randDnaSeqGetBasesTest:\tPos: " + randPos + "\t" + "Len: " + randLen + "\t'" + basesOri + "'\t=\t'" + basesBin + "'");
			}
		}
	}

	/**
	 * Create random sequences and store them in a DnaNSequence. 
	 * Compare getting a single random base from the original and DnaNSequence sequences.
	 * @param numTests
	 * @param lenMask
	 * @param seed
	 */
	public void randDnaSeqGetBaseTest(int numTests, int numTestsPerSeq, int lenMask, long seed) {
		Random rand = new Random(seed);

		for( int t = 0; t < numTests; t++ ) {
			String seq = "";
			int len = (rand.nextInt() & lenMask) + 10; // Randomly select sequence length
			seq = randSeq(len, rand); // Create a random sequence

			if( verbose ) System.out.println("DnaNSequence test:" + t + "\tlen:" + len + "\t" + seq);
			DnaNSequence bseq = new DnaNSequence(seq);

			// Retrieve numTestsPerSeq random bases from the sequence
			for( int i = 0; i < numTestsPerSeq; i++ ) {
				int randPos = rand.nextInt(len);
				char baseOri = seq.charAt(randPos);
				char baseBin = bseq.getBase(randPos);
				Assert.assertEquals(baseOri, baseBin);
			}
		}
	}

	/**
	 * Create random sequences and compare to storing them in a DnaNSequence
	 * @param numTests
	 * @param lenMask
	 * @param seed
	 */
	public void randDnaSeqTest(int numTests, int lenMask, long seed) {
		Random rand = new Random(seed);

		for( int t = 0; t < numTests; t++ ) {
			String seq = "";
			int len = (rand.nextInt() & lenMask) + 10; // Randomly select sequence length
			seq = randSeq(len, rand); // Create a random sequence

			if( verbose ) System.out.println("DnaNSequence test:" + t + "\tlen:" + len + "\t" + seq);
			DnaNSequence bseq = new DnaNSequence(seq);
			Assert.assertEquals(seq, bseq.toString());
		}
	}

	/**
	 * Create random sequences and store them in a DnaNSequence. 
	 * Compare after replacing random bases from the original and DnaNSequence sequences.
	 * @param numTests
	 * @param lenMask
	 * @param seed
	 */
	public void randReplaceBaseTest(int numTests, int numTestsPerSeq, int lenMask, long seed) {
		Random rand = new Random(seed);

		for( int t = 0; t < numTests; t++ ) {
			String seq = "";
			int len = (rand.nextInt() & lenMask) + 10; // Randomly select sequence length
			seq = randSeq(len, rand); // Create a random sequence

			DnaNSequence bseq = new DnaNSequence(seq);

			// Replace numTestsPerSeq random bases from the sequence
			if( verbose ) System.out.println("randReplaceBaseTest\nOri    :\t" + seq);
			for( int i = 0; i < numTestsPerSeq; i++ ) {
				// Random position
				int randPos = rand.nextInt(len);
				char baseOri = seq.charAt(randPos);

				// Random base (different than baseOri)
				char randBase = baseOri;
				while(randBase == baseOri)
					randBase = randBase(rand);

				// Replace base in sequence (string)
				char seqChars[] = seq.toCharArray();
				seqChars[randPos] = randBase;
				seq = new String(seqChars);

				// Replace i DnaNSequence
				bseq.setBase(randPos, randBase);
				if( verbose ) System.out.println("Changed:\t" + seq + "\tpos: " + randPos + "\trandbase: " + randBase + "\n\t\t" + bseq);

				// Compare results
				Assert.assertEquals(seq, bseq.toString());
			}
		}
	}

	/**
	 * Create a random sequence of length 'len'
	 * @param len
	 * @param rand
	 * @return
	 */
	String randSeq(int len, Random rand) {
		StringBuilder sb = new StringBuilder();
		// Create a random sequence
		for( int i = 0; i < len; i++ )
			sb.append(randBase(rand));
		return sb.toString();
	}

	public void test_01_short() {
		long seed = 20100615;
		int lenMask = 0xff;
		int numTests = 1000;
		randDnaSeqTest(numTests, lenMask, seed);
	}

	public void test_01_short_getBase() {
		long seed = 20110217;
		int lenMask = 0xff;
		int numTests = 1000;
		int numTestsPerSeq = 100;
		randDnaSeqGetBaseTest(numTests, numTestsPerSeq, lenMask, seed);
	}

	public void test_01_short_getBases() {
		long seed = 20110218;
		int lenMask = 0xff;
		int numTests = 1000;
		int numTestsPerSeq = 100;
		randDnaSeqGetBasesTest(numTests, numTestsPerSeq, lenMask, seed);
	}

	public void test_01_short_replaceBase() {
		long seed = 20110218;
		int lenMask = 0xff;
		int numTests = 1000;
		int numTestsPerSeq = 100;
		randReplaceBaseTest(numTests, numTestsPerSeq, lenMask, seed);
	}

	public void test_02_long() {
		long seed = 20100614;
		int lenMask = 0xffff;
		int numTests = 10;
		randDnaSeqTest(numTests, lenMask, seed);
	}

	public void test_02_long_getBase() {
		long seed = 20110217;
		int lenMask = 0xffff;
		int numTests = 10;
		int numTestsPerSeq = 1000;
		randDnaSeqGetBaseTest(numTests, numTestsPerSeq, lenMask, seed);
	}

	public void test_02_long_getBases() {
		long seed = 20110218;
		int lenMask = 0xffff;
		int numTests = 10;
		int numTestsPerSeq = 1000;
		randDnaSeqGetBasesTest(numTests, numTestsPerSeq, lenMask, seed);
	}

	public void test_02_long_replaceBase() {
		long seed = 20110217;
		int lenMask = 0xffff;
		int numTests = 10;
		int numTestsPerSeq = 1000;
		randReplaceBaseTest(numTests, numTestsPerSeq, lenMask, seed);
	}

	public void test_13_reverseWc() {
		long seed = 20100615;
		int lenMask = 0xfff;
		int numTests = 1000;

		Random rand = new Random(seed);

		for( int t = 0; t < numTests; t++ ) {
			int len = (rand.nextInt() & lenMask) + 10; // Randomly select sequence length
			String seq = randSeq(len, rand); // Create a random sequence
			String seqRwc = GprSeq.reverseWc(seq);

			DnaNSequence bseq = new DnaNSequence(seq);
			DnaNSequence rwc = (DnaNSequence) bseq.reverseWc();

			assertEquals(seqRwc.toUpperCase(), rwc.getSequence().toUpperCase());
		}
	}

}
