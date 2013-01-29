package ca.mcgill.mcb.pcingola.snpEffect.testCases;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Invoke all test cases for SnpEff
 * 
 * @author pcingola
 */
public class TestSuiteAll {

	public static void main(String args[]) {
		junit.textui.TestRunner.run(suite());
	}

	public static Test suite() {
		TestSuite suite = new TestSuite();

		// Stats
		suite.addTestSuite(TestCasesHypergeometric.class);
		suite.addTestSuite(TestCasesIntStats.class);
		suite.addTestSuite(TestCochranArmitage.class);

		// Binary sequences
		suite.addTestSuite(TestCasesNmers.class);
		suite.addTestSuite(TestCasesDnaSequence.class);
		suite.addTestSuite(TestCasesDnaSequenceByte.class);
		suite.addTestSuite(TestCasesDnaNSequence.class);
		suite.addTestSuite(TestCaseSequenceIndexer.class);
		suite.addTestSuite(TestCaseOverlap.class);
		suite.addTestSuite(TestCasesDnaOverlap.class);

		// Alignment
		suite.addTestSuite(TestCasesNeedlemanWunsch.class);

		// Intervals
		suite.addTestSuite(TestCasesIntervals.class);
		suite.addTestSuite(TestCasesMarkerUtils.class);

		// Codon tables
		suite.addTestSuite(TestCasesCodonTable.class);

		// SeqChange
		suite.addTestSuite(TestCasesSeqChange.class);
		suite.addTestSuite(TestCasesTranscript.class);
		suite.addTestSuite(TestCasesSnp.class);
		suite.addTestSuite(TestCasesMnp.class);
		suite.addTestSuite(TestCasesIns.class);
		suite.addTestSuite(TestCasesDel.class);
		suite.addTestSuite(TestCasesIntervalSeqChange.class);

		suite.addTestSuite(TestCasesSnpEnsembl.class);
		suite.addTestSuite(TestCasesMissenseSilentRatio.class);

		// Filter transcripts
		suite.addTestSuite(TestCasesFilterTranscripts.class);

		// File format
		suite.addTestSuite(TestCasesGff3.class);
		suite.addTestSuite(TestCasesGtf22.class);
		suite.addTestSuite(TestCasesVcf.class);
		suite.addTestSuite(TestCasesRefSeq.class);

		// File 
		suite.addTestSuite(TestCasesSeekableReader.class);

		// Protein coding sequences
		suite.addTestSuite(TestCasesProtein.class);

		// Loss of function and Nonsense mediated decay
		suite.addTestSuite(TestCasesLof.class);
		suite.addTestSuite(TestCasesNmd.class);

		return suite;
	}
}
