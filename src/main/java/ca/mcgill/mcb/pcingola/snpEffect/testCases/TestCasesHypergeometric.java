package ca.mcgill.mcb.pcingola.snpEffect.testCases;

import java.util.Random;

import junit.framework.TestCase;
import ca.mcgill.mcb.pcingola.stats.FisherExactTest;
import ca.mcgill.mcb.pcingola.stats.Hypergeometric;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * Test for Hypergeometric distribution and Fisher exact test
 * 
 * @author pcingola
 *
 */
public class TestCasesHypergeometric extends TestCase {

	boolean verbose = true;
	double threshold = 0.01;
	int numTests = 100;
	int MAX = 1000;
	Random rand;

	public TestCasesHypergeometric() {
		super();
		initRand();
	}

	/**
	 * Compare a result form a Fisher exact test (lower tail)
	 * @param k
	 * @param N
	 * @param D
	 * @param n
	 * @param result
	 */
	void compareFisherDown(int k, int N, int D, int n, double result) {
		double p = FisherExactTest.get().pValueDown(k, N, D, n, threshold);

		double abs = Math.abs(p - result);
		double diff = abs / Math.min(p, result);

		if( (abs > 1E-300) && (diff > 0.00001) ) {
			String err = "Difference:" + diff //
					+ "\n\t\tpValue:\t" + p //
					+ "\n\tExpected:\t" + result //
					+ "\n\tR: " + FisherExactTest.get().toR(k, N, D, n, true);
			Gpr.debug(err);
			throw new RuntimeException(err);
		}
	}

	/**
	 * Compare a result form a Fisher exact test (upper tail)
	 * @param k
	 * @param N
	 * @param D
	 * @param n
	 * @param result
	 */
	void compareFisherUp(int k, int N, int D, int n, double result) {
		double p = FisherExactTest.get().pValueUp(k, N, D, n, threshold);

		double abs = Math.abs(p - result);
		double diff = abs / Math.min(p, result);
		if( (abs > 1E-300) && (diff > 0.00001) ) {
			String err = "Difference:" + diff //
					+ "\n\t\tpValue:\t" + p //
					+ "\n\tExpected:\t" + result //
					+ "\n\tR: " + FisherExactTest.get().toR(k, N, D, n, false);
			Gpr.debug(err);
			throw new RuntimeException(err);
		}
	}

	void compareHypergeometric(int k, int N, int D, int n, double result) {
		double p = Hypergeometric.get().hypergeometric(k, N, D, n);

		double abs = Math.abs(p - result);
		double diff = abs / Math.min(p, result);
		if( (abs > 1E-300) && (diff > 0.00001) ) throw new RuntimeException("Difference:" + diff + "\t\t" + p + " != " + result);
	}

	/**
	 * 
	 */
	public void generate_test() {
		boolean lowerTail = true;
		for( int i = 0; i < numTests; ) {
			int N = rand.nextInt(MAX) + 1;
			int D = rand.nextInt(N) + 1;
			int n = rand.nextInt(N) + 1;
			int k = Math.max(rand.nextInt(Math.min(n, D)), 1);

			if( Hypergeometric.get().checkHypergeometricParams(k, N, D, n) ) {
				double p = FisherExactTest.get().pValueDown(k, N, D, n, threshold);

				if( p < threshold ) {
					System.out.print("\t print ( paste( 'compareFisher( " + k + ", " + N + ", " + D + ", " + n + ", ' , " + FisherExactTest.get().toR(k, N, D, n, lowerTail) + " , ');' ) );");
					System.out.println("");
					i++;
				}
			}
		}

		// assertEquals(90, minusInt.getEnd());
	}

	void initRand() {
		rand = new Random(20110124);
	}

	public void test_0() {
		// generate_test();
	}

	public void test_01_hyper() {
		// Compare to values calculated using other programs (R)
		compareHypergeometric(1, 19, 2, 12, 0.4912281);
		compareHypergeometric(1, 70, 51, 1, 0.7285714);
		compareHypergeometric(7, 73, 9, 38, 0.0773475);
		compareHypergeometric(4, 45, 11, 5, 0.00918348);
		compareHypergeometric(4, 33, 17, 14, 0.02327653);
		compareHypergeometric(3, 28, 5, 7, 0.07478632);
		compareHypergeometric(1, 73, 31, 2, 0.4954338);
		compareHypergeometric(1, 2, 1, 2, 1);
		compareHypergeometric(1, 10, 1, 3, 0.3);
		compareHypergeometric(1, 78, 57, 11, 2.584604e-06);
		compareHypergeometric(3, 45, 17, 23, 0.0005133975);
		compareHypergeometric(4, 86, 59, 11, 0.01652058);
		compareHypergeometric(1, 28, 5, 22, 0.003357753);
		compareHypergeometric(1, 42, 4, 30, 0.05896542);
		compareHypergeometric(12, 95, 13, 49, 0.001216176);
		compareHypergeometric(1, 45, 2, 27, 0.4909091);
		compareHypergeometric(1, 38, 3, 14, 0.458037);
		compareHypergeometric(1, 43, 22, 5, 0.1367861);
		compareHypergeometric(1, 25, 2, 6, 0.38);
		compareHypergeometric(7, 80, 34, 13, 0.1598987);
		compareHypergeometric(1, 78, 34, 3, 0.4227877);
		compareHypergeometric(1, 88, 9, 70, 5.361089e-06);
		compareHypergeometric(1, 13, 11, 2, 0.2820513);
		compareHypergeometric(1, 36, 4, 22, 0.1359477);
		compareHypergeometric(1, 20, 6, 7, 0.2324303);
		compareHypergeometric(1, 21, 9, 4, 0.3308271);
		compareHypergeometric(1, 39, 6, 1, 0.1538462);
		compareHypergeometric(5, 19, 9, 12, 0.3000714);
		compareHypergeometric(1, 1, 1, 1, 1);
		compareHypergeometric(1, 10, 5, 3, 0.4166667);
		compareHypergeometric(1, 86, 2, 67, 0.34829);
		compareHypergeometric(1, 21, 4, 6, 0.4561404);
		compareHypergeometric(29, 92, 43, 55, 0.06424561);
		compareHypergeometric(1, 23, 2, 21, 0.1660079);
		compareHypergeometric(5, 72, 7, 24, 0.03254646);
		compareHypergeometric(19, 64, 29, 54, 0.0001322347);
		compareHypergeometric(1, 40, 3, 35, 0.0354251);
		compareHypergeometric(6, 18, 7, 13, 0.2696078);
		compareHypergeometric(1, 8, 2, 5, 0.5357143);
		compareHypergeometric(15, 80, 31, 21, 0.0004160716);
		compareHypergeometric(3, 32, 18, 13, 0.002351405);
		compareHypergeometric(13, 35, 22, 18, 0.1410843);
		compareHypergeometric(7, 63, 45, 15, 0.01625885);
		compareHypergeometric(1, 17, 14, 4, 0.005882353);
		compareHypergeometric(2, 20, 6, 15, 0.01354489);
		compareHypergeometric(1, 10, 8, 2, 0.3555556);
		compareHypergeometric(1, 2, 1, 2, 1);
		compareHypergeometric(1, 94, 2, 93, 0.02127660);
		compareHypergeometric(11, 49, 22, 25, 0.2238699);
		compareHypergeometric(11, 95, 20, 27, 0.003754185);
		compareHypergeometric(1, 21, 9, 2, 0.5142857);
		compareHypergeometric(1, 6, 2, 2, 0.5333333);
		compareHypergeometric(2, 18, 9, 7, 0.1425339);
		compareHypergeometric(16, 61, 40, 36, 1.500428e-05);
		compareHypergeometric(2, 54, 7, 5, 0.1076724);
		compareHypergeometric(1, 89, 1, 47, 0.5280899);
		compareHypergeometric(2, 26, 5, 13, 0.3391304);
		compareHypergeometric(1, 75, 1, 24, 0.32);
		compareHypergeometric(2, 21, 10, 12, 0.001684074);
		compareHypergeometric(1, 64, 8, 7, 0.4181274);
		compareHypergeometric(16, 89, 45, 48, 0.0003430746);
		compareHypergeometric(26, 56, 39, 37, 0.2367166);
		compareHypergeometric(10, 30, 16, 11, 0.002052307);
		compareHypergeometric(11, 92, 30, 68, 4.373004e-08);
		compareHypergeometric(17, 61, 46, 19, 0.06186141);
		compareHypergeometric(15, 70, 28, 17, 4.092152e-06);
		compareHypergeometric(6, 55, 19, 33, 0.001963553);
		compareHypergeometric(18, 71, 57, 31, 3.427275e-05);
		compareHypergeometric(48, 95, 62, 66, 0.01405199);
		compareHypergeometric(3, 46, 8, 34, 0.01816302);
		compareHypergeometric(3, 26, 8, 20, 0.004378230);
		compareHypergeometric(4, 69, 9, 47, 0.08288145);
		compareHypergeometric(13, 78, 57, 28, 0.0001087879);
		compareHypergeometric(20, 68, 24, 37, 0.000332827);
		compareHypergeometric(5, 89, 7, 23, 0.01047522);
		compareHypergeometric(3, 12, 5, 10, 0.1515152);
		compareHypergeometric(16, 65, 29, 47, 0.005136288);
		compareHypergeometric(5, 76, 12, 6, 0.0002318555);
		compareHypergeometric(2, 43, 15, 8, 0.2727957);
		compareHypergeometric(3, 46, 12, 35, 9.25133e-06);
		compareHypergeometric(3, 92, 18, 4, 0.02161083);
		compareHypergeometric(2, 87, 47, 7, 0.1217291);
		compareHypergeometric(1, 26, 6, 11, 0.1434783);
		compareHypergeometric(1, 36, 3, 28, 0.1098039);
		compareHypergeometric(1, 34, 1, 24, 0.7058824);
		compareHypergeometric(1, 88, 1, 81, 0.9204545);
		compareHypergeometric(1, 2, 2, 1, 1);
		compareHypergeometric(5, 23, 12, 7, 0.1776821);
		compareHypergeometric(1, 21, 3, 2, 0.2571429);
		compareHypergeometric(4, 83, 30, 55, 3.705065e-15);
		compareHypergeometric(1, 63, 1, 61, 0.968254);
		compareHypergeometric(1, 99, 20, 5, 0.420144);
		compareHypergeometric(1, 73, 2, 36, 0.5068493);
		compareHypergeometric(1, 65, 38, 2, 0.4932692);
		compareHypergeometric(1, 54, 1, 51, 0.9444444);
		compareHypergeometric(1, 5, 1, 4, 0.8);
		compareHypergeometric(1, 1, 1, 1, 1);
		compareHypergeometric(1, 48, 1, 10, 0.2083333);
		compareHypergeometric(1, 71, 63, 2, 0.2028169);
		compareHypergeometric(1, 46, 6, 20, 0.1404532);
		compareHypergeometric(1, 70, 1, 65, 0.9285714);
		compareHypergeometric(31, 70, 34, 60, 0.1259187);
		compareHypergeometric(1, 17, 1, 5, 0.2941176);
		compareHypergeometric(5, 30, 9, 13, 0.2140930);
		compareHypergeometric(1, 26, 4, 15, 0.1655518);
		compareHypergeometric(30, 93, 51, 57, 0.1471864);
		compareHypergeometric(11, 17, 13, 13, 0.1966387);
		compareHypergeometric(1, 25, 4, 2, 0.28);
		compareHypergeometric(1, 5, 4, 1, 0.8);
		compareHypergeometric(27, 89, 40, 28, 5.527519e-12);
		compareHypergeometric(9, 53, 25, 25, 0.06875583);
		compareHypergeometric(2, 18, 12, 6, 0.05332902);
		compareHypergeometric(1, 40, 21, 18, 3.167218e-08);
	}

	public void test_02_hyper() {
		// Compare to values calculated using other programs (R)
		compareHypergeometric(57, 470, 141, 281, 1.507456e-08);
		compareHypergeometric(152, 912, 754, 203, 0.0004373848);
		compareHypergeometric(44, 682, 324, 82, 0.04638736);
		compareHypergeometric(21, 373, 32, 294, 0.03056803);
		compareHypergeometric(59, 545, 136, 95, 7.273835e-18);
		compareHypergeometric(36, 345, 190, 41, 2.082751e-06);
		compareHypergeometric(46, 770, 68, 535, 0.1018642);
		compareHypergeometric(129, 833, 418, 241, 0.02855101);
		compareHypergeometric(3, 28, 5, 7, 0.07478632);
		compareHypergeometric(14, 873, 552, 42, 4.038633e-05);
		compareHypergeometric(97, 202, 133, 105, 1.128978e-17);
		compareHypergeometric(152, 610, 181, 333, 8.438052e-23);
		compareHypergeometric(135, 444, 142, 138, 6.925182e-101);
		compareHypergeometric(7, 486, 25, 51, 0.00781056);
		compareHypergeometric(1, 42, 4, 30, 0.05896542);
		compareHypergeometric(41, 145, 82, 72, 0.1323862);
		compareHypergeometric(152, 738, 203, 446, 2.215361e-07);
		compareHypergeometric(16, 911, 728, 164, 1.023263e-115);
		compareHypergeometric(37, 243, 54, 187, 0.03645426);
		compareHypergeometric(2, 152, 51, 41, 8.419947e-07);
		compareHypergeometric(117, 480, 274, 253, 1.791619e-07);
		compareHypergeometric(247, 578, 388, 283, 3.29638e-25);
		compareHypergeometric(429, 888, 817, 446, 1.650871e-06);
		compareHypergeometric(1, 913, 749, 4, 0.01879589);
		compareHypergeometric(4, 236, 192, 10, 0.003181798);
		compareHypergeometric(1, 620, 6, 207, 0.2627278);
		compareHypergeometric(269, 351, 311, 287, 1.610554e-08);
		compareHypergeometric(27, 521, 155, 357, 3.076794e-60);
		compareHypergeometric(7, 461, 151, 184, 7.18484e-32);
		compareHypergeometric(5, 924, 236, 23, 0.1826543);
		compareHypergeometric(30, 601, 236, 60, 0.0224787);
		compareHypergeometric(29, 510, 305, 203, 5.85835e-71);
		compareHypergeometric(29, 92, 43, 55, 0.06424561);
		compareHypergeometric(12, 634, 42, 295, 0.006764663);
		compareHypergeometric(13, 237, 88, 156, 2.39519e-40);
		compareHypergeometric(91, 523, 106, 234, 1.382574e-22);
		compareHypergeometric(89, 572, 403, 188, 6.492409e-17);
		compareHypergeometric(35, 964, 146, 102, 1.617078e-07);
		compareHypergeometric(23, 470, 58, 108, 0.001111836);
		compareHypergeometric(65, 755, 119, 162, 1.070822e-18);
		compareHypergeometric(15, 380, 151, 21, 0.002024582);
		compareHypergeometric(24, 173, 63, 55, 0.0544751);
		compareHypergeometric(29, 707, 33, 649, 0.1576592);
		compareHypergeometric(27, 335, 172, 173, 4.587644e-46);
		compareHypergeometric(8, 117, 36, 18, 0.08555023);
		compareHypergeometric(17, 206, 36, 50, 0.0005816399);
		compareHypergeometric(1, 10, 8, 2, 0.3555556);
		compareHypergeometric(33, 149, 48, 87, 0.03005922);
		compareHypergeometric(1, 402, 343, 30, 1.218651e-26);
		compareHypergeometric(109, 521, 110, 472, 6.95748e-05);
		compareHypergeometric(66, 306, 122, 206, 3.471764e-05);
		compareHypergeometric(16, 61, 40, 36, 1.500428e-05);
		compareHypergeometric(32, 554, 163, 301, 3.947078e-27);
		compareHypergeometric(30, 789, 310, 52, 0.002491050);
		compareHypergeometric(42, 726, 155, 45, 4.110391e-27);
		compareHypergeometric(293, 803, 481, 450, 0.0001804299);
		compareHypergeometric(131, 375, 151, 249, 1.094179e-12);
		compareHypergeometric(210, 321, 274, 246, 0.1475397);
		compareHypergeometric(8, 789, 388, 277, 3.078242e-96);
		compareHypergeometric(260, 556, 539, 277, 5.566594e-06);
		compareHypergeometric(82, 130, 116, 91, 0.2064143);
		compareHypergeometric(3, 461, 374, 7, 0.02295504);
		compareHypergeometric(2, 890, 10, 796, 4.312254e-07);
		compareHypergeometric(42, 770, 308, 367, 6.205572e-58);
		compareHypergeometric(6, 55, 19, 33, 0.001963553);
		compareHypergeometric(4, 171, 92, 47, 3.491095e-14);
		compareHypergeometric(310, 995, 327, 976, 3.646995e-07);
		compareHypergeometric(2, 184, 34, 28, 0.05393918);
		compareHypergeometric(82, 669, 183, 338, 0.01339627);
		compareHypergeometric(1, 278, 3, 152, 0.3379183);
		compareHypergeometric(183, 868, 584, 305, 0.0002243882);
		compareHypergeometric(13, 112, 93, 30, 2.834419e-10);
		compareHypergeometric(26, 143, 33, 81, 0.001992269);
		compareHypergeometric(19, 146, 20, 71, 3.193124e-06);
		compareHypergeometric(73, 392, 154, 304, 4.244224e-32);
		compareHypergeometric(2, 87, 47, 7, 0.1217291);
		compareHypergeometric(13, 336, 15, 256, 0.1726304);
		compareHypergeometric(1, 385, 30, 144, 7.79075e-06);
		compareHypergeometric(35, 488, 97, 57, 1.079242e-13);
		compareHypergeometric(2, 302, 168, 7, 0.1104994);
		compareHypergeometric(5, 23, 12, 7, 0.1776821);
		compareHypergeometric(119, 983, 628, 145, 1.421082e-07);
		compareHypergeometric(324, 663, 403, 541, 0.05047842);
		compareHypergeometric(1, 99, 20, 5, 0.420144);
		compareHypergeometric(1, 741, 671, 2, 0.1713171);
		compareHypergeometric(1, 73, 2, 36, 0.5068493);
		compareHypergeometric(47, 554, 101, 407, 8.9411e-11);
		compareHypergeometric(1, 152, 59, 4, 0.3582004);
		compareHypergeometric(18, 505, 326, 129, 2.489987e-44);
		compareHypergeometric(45, 148, 53, 50, 5.759724e-24);
		compareHypergeometric(19, 446, 46, 26, 1.422818e-15);
		compareHypergeometric(1, 178, 3, 142, 0.09679974);
		compareHypergeometric(607, 770, 734, 620, 2.162936e-09);
		compareHypergeometric(185, 546, 206, 407, 1.97752e-11);
		compareHypergeometric(136, 917, 169, 294, 7.691036e-48);
		compareHypergeometric(31, 321, 108, 43, 2.813547e-08);
		compareHypergeometric(16, 515, 475, 39, 3.418161e-19);
		compareHypergeometric(91, 926, 340, 335, 1.576175e-06);
		compareHypergeometric(97, 396, 167, 106, 1.446200e-35);
		compareHypergeometric(7, 305, 14, 51, 0.002821931);
	}

	public void test_03_fisher() {
		compareFisherUp(59, 545, 136, 95, 8.28958173422445e-18);
		compareFisherUp(36, 345, 190, 41, 2.40265087580901e-06);
		compareFisherUp(97, 202, 133, 105, 1.18466240918432e-17);
		compareFisherUp(152, 610, 181, 333, 9.77424907183547e-23);
		compareFisherUp(135, 444, 142, 138, 6.928747246508e-101);
		compareFisherUp(152, 738, 203, 446, 3.67088986300015e-07);
		compareFisherUp(247, 578, 388, 283, 3.79404002619577e-25);
		compareFisherUp(429, 888, 817, 446, 2.26456109642847e-06);
		compareFisherUp(269, 351, 311, 287, 1.82978544974049e-08);
		compareFisherUp(91, 523, 106, 234, 1.50955326758971e-22);
		compareFisherUp(35, 964, 146, 102, 2.21430371305547e-07);
		compareFisherUp(23, 470, 58, 108, 0.00174305424896066);
		compareFisherUp(65, 755, 119, 162, 1.25368328189431e-18);
		compareFisherUp(15, 380, 151, 21, 0.00257864863363782);
		compareFisherUp(17, 206, 36, 50, 0.000767632652092611);
		compareFisherUp(109, 521, 110, 472, 7.4260450370284e-05);
		compareFisherUp(30, 789, 310, 52, 0.00423623853376874);
		compareFisherUp(42, 726, 155, 45, 4.16785300462029e-27);
		compareFisherUp(293, 803, 481, 450, 0.000438131577741957);
		compareFisherUp(131, 375, 151, 249, 1.30995155177008e-12);
		compareFisherUp(26, 143, 33, 81, 0.00262084046420338);
		compareFisherUp(19, 146, 20, 71, 3.30381924570966e-06);
		compareFisherUp(35, 488, 97, 57, 1.20102019250494e-13);
		compareFisherUp(119, 983, 628, 145, 2.10926932575014e-07);
		compareFisherUp(45, 148, 53, 50, 5.81512027352164e-24);
		compareFisherUp(19, 446, 46, 26, 1.45759552943497e-15);
		compareFisherUp(607, 770, 734, 620, 2.4342390531863e-09);
		compareFisherUp(185, 546, 206, 407, 2.49420910249361e-11);
		compareFisherUp(136, 917, 169, 294, 8.2186630995667e-48);
		compareFisherUp(31, 321, 108, 43, 3.27189426625206e-08);
		compareFisherUp(97, 396, 167, 106, 1.48934787720076e-35);
		compareFisherUp(7, 305, 14, 51, 0.00331494622577636);
		compareFisherUp(71, 835, 75, 483, 7.2876824125444e-14);
		compareFisherUp(37, 461, 70, 181, 0.00876705910701063);
		compareFisherUp(198, 761, 318, 432, 0.00580370806914619);
		compareFisherUp(133, 936, 378, 195, 8.99332118627283e-19);
		compareFisherUp(23, 291, 65, 50, 3.31838151428692e-05);
		compareFisherUp(63, 560, 67, 190, 2.1131322965245e-28);
		compareFisherUp(44, 115, 71, 49, 3.89847129521515e-08);
		compareFisherUp(153, 800, 195, 527, 9.21709696267346e-06);
		compareFisherUp(37, 397, 67, 156, 0.00285215778284655);
		compareFisherUp(94, 433, 212, 130, 1.18080726821783e-10);
		compareFisherUp(629, 915, 873, 636, 5.48138806189002e-13);
		compareFisherUp(225, 687, 254, 282, 3.86407484833323e-91);
		compareFisherUp(187, 550, 221, 305, 2.28319162936954e-31);
		compareFisherUp(183, 812, 194, 725, 0.00473364047467412);
		compareFisherUp(75, 531, 142, 193, 1.99243069292124e-06);
		compareFisherUp(133, 897, 158, 661, 0.00044883940423697);
		compareFisherUp(168, 973, 227, 599, 5.28590254525447e-06);
		compareFisherUp(266, 492, 375, 276, 5.17872981915267e-35);
		compareFisherUp(142, 814, 161, 579, 1.1947296127053e-08);
		compareFisherUp(114, 713, 149, 177, 3.80835210298479e-54);
		compareFisherUp(49, 174, 61, 92, 5.76190646498642e-08);
		compareFisherUp(70, 361, 243, 88, 0.00300291746372154);
		compareFisherUp(11, 417, 31, 56, 0.00098764195265165);
		compareFisherUp(36, 537, 157, 41, 4.94220221304424e-16);
		compareFisherUp(115, 471, 182, 175, 1.56291175266919e-20);
		compareFisherUp(495, 881, 523, 763, 2.98889043133355e-17);
		compareFisherUp(241, 674, 307, 457, 3.14514846186333e-08);
		compareFisherUp(200, 443, 303, 242, 1.08062874506896e-12);
		compareFisherUp(353, 987, 554, 366, 3.06547860908313e-101);
		compareFisherUp(11, 300, 71, 20, 0.00172542325090638);
		compareFisherUp(42, 252, 204, 44, 0.00318184189489258);
		compareFisherUp(51, 305, 80, 53, 1.47808357334896e-34);
		compareFisherUp(34, 313, 176, 46, 0.0062078519294078);
		compareFisherUp(24, 373, 25, 268, 0.0021702822049755);
		compareFisherUp(47, 677, 141, 48, 3.36388120818521e-34);
		compareFisherUp(41, 328, 65, 89, 5.84598356867524e-12);
		compareFisherUp(65, 559, 190, 81, 2.05254736770733e-20);
		compareFisherUp(75, 984, 299, 165, 5.27613211166751e-06);
		compareFisherUp(92, 270, 95, 240, 0.00107752920684479);
		compareFisherUp(203, 969, 612, 220, 1.65237135169329e-28);
		compareFisherUp(168, 255, 194, 194, 3.59390288867635e-11);
		compareFisherUp(25, 395, 74, 31, 4.64779278697149e-15);
		compareFisherUp(154, 621, 171, 426, 3.32248893082949e-14);
		compareFisherUp(12, 327, 148, 15, 0.00569010303563112);
		compareFisherUp(80, 772, 97, 245, 1.60368435504478e-28);
		compareFisherUp(21, 654, 138, 25, 1.04540516516679e-11);
		compareFisherUp(119, 432, 121, 285, 8.46674711440474e-24);
		compareFisherUp(200, 672, 245, 202, 3.72458025046742e-123);
		compareFisherUp(13, 675, 19, 137, 5.33314581564644e-06);
		compareFisherUp(48, 299, 70, 119, 2.51027306645499e-08);
		compareFisherUp(197, 694, 250, 432, 4.58839930028976e-12);
		compareFisherUp(103, 823, 210, 295, 3.59407139405134e-06);
		compareFisherUp(73, 235, 108, 102, 3.222899986319e-12);
		compareFisherUp(565, 848, 684, 623, 2.88228387511532e-31);
		compareFisherUp(364, 582, 426, 384, 6.72059678081954e-61);
		compareFisherUp(39, 781, 42, 145, 2.07138402016061e-27);
		compareFisherUp(179, 371, 184, 302, 1.23670887007989e-16);
		compareFisherUp(159, 541, 459, 170, 4.77152383854613e-05);
		compareFisherUp(137, 675, 190, 233, 5.40690104095474e-37);
		compareFisherUp(360, 800, 406, 372, 1.17796424900086e-155);
		compareFisherUp(281, 736, 303, 620, 3.36267010147432e-08);
		compareFisherUp(504, 795, 644, 606, 0.0044212941820773);
		compareFisherUp(141, 450, 183, 209, 5.03306096804749e-28);
		compareFisherUp(247, 968, 390, 299, 2.01859193466352e-74);
		compareFisherUp(25, 172, 51, 54, 0.00131684344612187);
		compareFisherUp(54, 757, 220, 71, 6.58948512308768e-18);
		compareFisherUp(70, 301, 180, 97, 0.00172465566053411);
		compareFisherUp(194, 683, 500, 204, 1.56219065349161e-20);
	}

	public void test_04_fisher() {
		compareFisherDown(57, 470, 141, 281, 6.6866974987128e-09);
		compareFisherDown(152, 912, 754, 203, 0.000440144803784442);
		compareFisherDown(14, 873, 552, 42, 1.37118944858872e-05);
		compareFisherDown(16, 911, 728, 164, 5.41965575482104e-118);
		compareFisherDown(2, 152, 51, 41, 5.37265914345806e-08);
		compareFisherDown(117, 480, 274, 253, 1.05549713246985e-07);
		compareFisherDown(4, 236, 192, 10, 0.000393391860642947);
		compareFisherDown(27, 521, 155, 357, 7.1546179587792e-62);
		compareFisherDown(7, 461, 151, 184, 2.67218847522338e-33);
		compareFisherDown(29, 510, 305, 203, 1.10542140757286e-72);
		compareFisherDown(13, 237, 88, 156, 1.71630737398072e-42);
		compareFisherDown(89, 572, 403, 188, 1.58248960811606e-17);
		compareFisherDown(27, 335, 172, 173, 1.00022714178528e-47);
		compareFisherDown(1, 402, 343, 30, 3.55291842099265e-29);
		compareFisherDown(66, 306, 122, 206, 1.87164960370351e-05);
		compareFisherDown(16, 61, 40, 36, 4.57273382026378e-07);
		compareFisherDown(32, 554, 163, 301, 4.8225896336559e-28);
		compareFisherDown(8, 789, 388, 277, 3.18813457308328e-98);
		// compareFisherDown(260, 556, 539, 277, 0);
		compareFisherDown(2, 890, 10, 796, 1.04769809538842e-08);
		compareFisherDown(42, 770, 308, 367, 4.38024039390092e-59);
		compareFisherDown(6, 55, 19, 33, 0.000296863410300295);
		compareFisherDown(4, 171, 92, 47, 1.31758038623838e-15);
		compareFisherDown(310, 995, 327, 976, 1.92919081042295e-08);
		compareFisherDown(183, 868, 584, 305, 0.000313284475334555);
		compareFisherDown(13, 112, 93, 30, 5.09345032728386e-12);
		compareFisherDown(73, 392, 154, 304, 1.16596572422222e-33);
		compareFisherDown(1, 385, 30, 144, 3.82323826836252e-07);
		compareFisherDown(47, 554, 101, 407, 2.47909349507427e-11);
		compareFisherDown(18, 505, 326, 129, 9.10158234925755e-46);
		compareFisherDown(16, 515, 475, 39, 8.60004121044127e-21);
		compareFisherDown(91, 926, 340, 335, 1.56179857832206e-06);
		compareFisherDown(49, 953, 389, 335, 2.84436360391136e-37);
		compareFisherDown(2, 763, 61, 554, 7.7377635214529e-36);
		compareFisherDown(1, 40, 21, 18, 1.67577669149464e-10);
		compareFisherDown(12, 772, 45, 579, 2.42945701198204e-13);
		compareFisherDown(22, 347, 82, 279, 1.58907702327229e-41);
		compareFisherDown(2, 446, 273, 111, 1.52663791438916e-57);
		compareFisherDown(46, 391, 323, 91, 1.8448381164756e-18);
		compareFisherDown(52, 556, 247, 163, 3.76060089022952e-05);
		compareFisherDown(2, 282, 189, 85, 3.85835151779573e-60);
		compareFisherDown(29, 130, 62, 89, 3.17545503454762e-08);
		compareFisherDown(63, 629, 201, 264, 6.5881507974015e-05);
		compareFisherDown(200, 732, 471, 397, 5.29538922431954e-19);
		compareFisherDown(49, 517, 403, 105, 9.86612011417267e-17);
		compareFisherDown(12, 530, 38, 376, 3.04971603958179e-08);
		compareFisherDown(18, 656, 67, 587, 2.68818158736483e-44);
		compareFisherDown(16, 242, 40, 187, 1.93661820051305e-09);
		compareFisherDown(89, 698, 546, 149, 1.14195550706473e-09);
		compareFisherDown(8, 556, 440, 95, 3.14204784057095e-68);
		compareFisherDown(5, 431, 44, 152, 3.39929178807456e-05);
		//compareFisherDown(70, 251, 118, 203, 0);
		compareFisherDown(32, 915, 156, 574, 1.77058007272689e-33);
		compareFisherDown(26, 627, 461, 135, 3.93976550346058e-55);
		compareFisherDown(59, 619, 460, 131, 4.10650753885808e-17);
		compareFisherDown(10, 399, 72, 155, 7.996152169811e-08);
		compareFisherDown(77, 445, 250, 169, 0.000141929399394228);
		compareFisherDown(60, 821, 522, 197, 1.24193470803979e-28);
		compareFisherDown(15, 627, 29, 495, 0.000176997186877898);
		compareFisherDown(51, 362, 296, 85, 4.83200953355619e-09);
		compareFisherDown(65, 736, 240, 460, 1.68636888045595e-44);
		compareFisherDown(8, 551, 35, 419, 5.79043383485544e-13);
		compareFisherDown(31, 181, 52, 141, 6.371377257703e-05);
		compareFisherDown(155, 931, 546, 410, 1.57141011213516e-31);
		compareFisherDown(49, 370, 120, 200, 0.000131263672590605);
		compareFisherDown(19, 215, 168, 47, 4.78138299321588e-12);
		compareFisherDown(10, 567, 424, 29, 4.27769623511992e-07);
		compareFisherDown(51, 372, 313, 107, 1.36720312068375e-34);
		compareFisherDown(13, 393, 155, 97, 3.42880198036781e-11);
		compareFisherDown(28, 452, 192, 278, 3.07029863124219e-81);
		compareFisherDown(34, 812, 594, 112, 8.01635155599195e-26);
		compareFisherDown(100, 637, 474, 179, 1.91831046845406e-11);
		compareFisherDown(5, 425, 279, 115, 6.3490979346593e-65);
		compareFisherDown(54, 838, 243, 263, 7.04361587049744e-05);
		compareFisherDown(3, 825, 376, 55, 8.19330754927712e-13);
		compareFisherDown(15, 109, 80, 35, 1.65593752618442e-07);
		compareFisherDown(73, 674, 543, 118, 2.97539829045907e-08);
		compareFisherDown(309, 699, 344, 650, 0.000310818727006591);
		compareFisherDown(46, 469, 85, 335, 4.61477161937294e-05);
		compareFisherDown(115, 551, 515, 143, 1.62070097654381e-12);
		compareFisherDown(3, 162, 8, 144, 1.98324541560583e-05);
		compareFisherDown(72, 528, 260, 271, 1.87125928237604e-28);
		compareFisherDown(21, 589, 439, 84, 1.26356949684138e-26);
		compareFisherDown(8, 997, 545, 86, 1.79712562732956e-21);
		compareFisherDown(1, 153, 5, 120, 0.000362873598189044);
		compareFisherDown(155, 476, 237, 380, 6.80742086958635e-17);
		compareFisherDown(12, 243, 32, 173, 3.46251414041836e-06);
		compareFisherDown(2, 69, 30, 20, 1.84777897818333e-05);
		compareFisherDown(10, 621, 137, 319, 1.2359893210442e-36);
		compareFisherDown(13, 908, 418, 101, 9.79596645147574e-15);
		compareFisherDown(59, 519, 227, 260, 1.35995478816319e-23);
		compareFisherDown(81, 866, 215, 463, 2.63951613872026e-08);
		compareFisherDown(2, 669, 482, 20, 2.29978215266672e-10);
		compareFisherDown(1, 255, 28, 87, 3.61797690526272e-06);
		compareFisherDown(25, 509, 204, 94, 0.000867040874821994);
		compareFisherDown(1, 128, 40, 25, 2.40252924022306e-05);
		compareFisherDown(189, 964, 348, 621, 3.35004927065627e-07);
		compareFisherDown(11, 158, 133, 23, 7.06584015934966e-07);
		compareFisherDown(65, 358, 114, 267, 8.50412154506876e-08);
		compareFisherDown(179, 775, 414, 507, 6.06818418807189e-49);
	}

	public void test_05_fisher() {
		compareFisherDown(1, 100, 50, 0, 1);
		compareFisherDown(1, 100, 0, 20, 1);
		compareFisherDown(0, 100, 50, 0, 0);
		compareFisherDown(0, 100, 0, 20, 0);
	}
}
