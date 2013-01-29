package ca.mcgill.mcb.pcingola.snpEffect.testCases;

import junit.framework.Assert;
import junit.framework.TestCase;
import ca.mcgill.mcb.pcingola.probablility.CochranArmitageTest;

/**
 * Cochran-Armitage test statistic test case
 * 
 * @author pcingola
 */
public class TestCochranArmitage extends TestCase {

	public void test_01() {
		int N1[] = { 20, 20, 20 };
		int N2[] = { 10, 20, 30 };

		double test = CochranArmitageTest.get().test(N1, N2, CochranArmitageTest.WEIGHT_DOMINANT);
		Assert.assertEquals(1.851, test, 0.001);
	}

	public void test_02() {
		int N1[] = { 20, 20, 20 };
		int N2[] = { 10, 20, 30 };
		double test = CochranArmitageTest.get().test(N1, N2, CochranArmitageTest.WEIGHT_RECESSIVE);
		Assert.assertEquals(-2.108, test, 0.001);
	}

	public void test_03() {
		int N1[] = { 20, 20, 20 };
		int N2[] = { 10, 20, 30 };
		double test = CochranArmitageTest.get().test(N1, N2, CochranArmitageTest.WEIGHT_CODOMINANT);
		System.out.println("Test: " + test);
		Assert.assertEquals(-2.284, test, 0.001);
	}

	public void test_04() {
		int n1[] = { 17066, 14464, 788, 126, 37 };
		int n2[] = { 48, 38, 5, 1, 1 };
		double w[] = { 1, 2, 3, 4, 5 };

		double p = CochranArmitageTest.get().p(n1, n2, w);
		Assert.assertEquals(0.088, p, 0.001);
	}

	public void test_05() {
		int n1[] = { 17066, 14464, 788, 126, 37 };
		int n2[] = { 48, 38, 5, 1, 1 };
		double w[] = { 0, 0.5, 1.5, 4, 8 };

		double p = CochranArmitageTest.get().p(n1, n2, w);
		Assert.assertEquals(0.0039, p, 0.000001);
	}

}
