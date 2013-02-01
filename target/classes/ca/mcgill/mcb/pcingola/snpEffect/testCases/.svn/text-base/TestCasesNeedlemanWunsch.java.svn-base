package ca.mcgill.mcb.pcingola.snpEffect.testCases;

import junit.framework.Assert;
import junit.framework.TestCase;
import ca.mcgill.mcb.pcingola.fileIterator.NeedlemanWunsch;

/**
 * test cases for Needleman Wunsch algorithm
 * 
 * @author pcingola
 */
public class TestCasesNeedlemanWunsch extends TestCase {

	public void test_01() {
		String as[] = { "TTT", "TTTGTT", "GCG", "G" };
		String bs[] = { "TTTGTT", "TTT", "G", "GCG" };
		String res[] = { "-GTT", "+GTT", "+CG", "-CG" };

		for( int i = 0; i < as.length; i++ ) {
			System.out.println("---------------------------------------- " + i + " ----------------------------------------");
			String a = as[i];
			String b = bs[i];
			NeedlemanWunsch nw = new NeedlemanWunsch(a, b);
			nw.align();
			System.out.println("a    : '" + a + "'\nb    : '" + b + "'\nAlign: '" + nw.getAlignment() + "'" + "\tOffset: " + nw.getOffset() + "\n");

			Assert.assertEquals(res[i], nw.getAlignment());
		}
	}
}
