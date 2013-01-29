package ca.mcgill.mcb.pcingola.snpEffect.testCases;

import java.util.List;
import java.util.Random;

import junit.framework.Assert;
import junit.framework.TestCase;
import ca.mcgill.mcb.pcingola.codons.CodonTable;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.snpEffect.factory.SnpEffPredictorFactoryRand;
import ca.mcgill.mcb.pcingola.util.GprSeq;

/**
 * Test random SNP changes 
 * 
 * @author pcingola
 */
public class TestCasesIns extends TestCase {

	boolean debug = false;
	Random rand;
	Config config;
	Genome genome;
	Chromosome chromosome;
	Gene gene;
	Transcript transcript;
	SnpEffectPredictor snpEffectPredictor;
	String chromoSequence = "";
	char chromoBases[];

	public TestCasesIns() {
		super();
		init();
	}

	void init() {
		initRand();
		initSnpEffPredictor();
	}

	void initRand() {
		rand = new Random(20100629);
	}

	void initSnpEffPredictor() {
		// Create a config and force out snpPredictor for hg37 chromosome Y
		config = new Config("testCase", Config.DEFAULT_CONFIG_FILE);

		// Create factory
		int maxGeneLen = 1000;
		int maxTranscripts = 1;
		int maxExons = 5;
		SnpEffPredictorFactoryRand sepf = new SnpEffPredictorFactoryRand(config, 1, rand, maxGeneLen, maxTranscripts, maxExons);

		// Chromosome sequence
		chromoSequence = sepf.getChromoSequence();
		chromoBases = chromoSequence.toCharArray();

		// Create predictor
		snpEffectPredictor = sepf.create();
		config.setSnpEffectPredictor(snpEffectPredictor);

		// No upstream or downstream
		config.getSnpEffectPredictor().setUpDownStreamLength(0);

		// Build forest
		config.getSnpEffectPredictor().buildForest();

		chromosome = sepf.getChromo();
		genome = config.getGenome();
		gene = genome.getGenes().iterator().next();
		transcript = gene.iterator().next();
	}

	public void test_01() {
		int N = 1000;
		CodonTable codonTable = genome.codonTable();

		// Test N times
		//	- Create a random gene transcript, exons
		//	- Create a random Insert at each position
		//	- Calculate effect
		for (int i = 0; i < N; i++) {
			initSnpEffPredictor();
			if (debug) System.out.println("INS Test iteration: " + i + "\n" + transcript);
			else System.out.println("INS Test iteration: " + i + "\t" + transcript.cds());

			int cdsBaseNum = 0;

			// For each exon...
			for (Exon exon : transcript.sortedStrand()) {
				int step = exon.getStrand() >= 0 ? 1 : -1;
				int beg = exon.getStrand() >= 0 ? exon.getStart() : exon.getEnd();

				// For each base in this exon...
				for (int pos = beg; (pos >= exon.getStart()) && (pos <= exon.getEnd()); pos += step, cdsBaseNum++) {
					// Get a random base different from 'refBase'
					int insLen = rand.nextInt(10) + 1;
					String insPlus = GprSeq.randSequence(rand, insLen); // Indesrtio (plus strand)
					String ins = insPlus;

					// Codon number
					int cdsCodonNum = cdsBaseNum / 3;
					int cdsCodonPos = cdsBaseNum % 3;

					int minCodonPos = cdsCodonNum * 3;
					int maxCodonPos = minCodonPos + 3;
					if (maxCodonPos < transcript.cds().length()) {
						String codonOld = transcript.cds().substring(minCodonPos, maxCodonPos);
						codonOld = codonOld.toUpperCase();
						String aaOld = codonTable.aa(codonOld);

						// Codon change
						String codonNew = "", aaNew = "";

						// Create a SeqChange
						int seqChangeStrand = rand.nextBoolean() ? +1 : -1;
						if (seqChangeStrand == -exon.getStrand()) ins = GprSeq.reverseWc(insPlus);
						SeqChange seqChange = new SeqChange(chromosome, pos, "", "+" + ins, seqChangeStrand, "", 1.0, 1);

						// Is it an insertion?
						Assert.assertEquals(true, seqChange.isIns());

						// Expected Effect
						String effectExpected = "";
						if (insLen % 3 != 0) {
							aaOld = "-";
							codonNew = insPlus;
							aaNew = codonTable.aa(codonNew);
							effectExpected = "FRAME_SHIFT(" + aaOld + "/" + aaNew + ")";
						} else {
							if (cdsCodonPos == 0) {
								codonOld = aaOld = "-";
								codonNew = insPlus;
								aaNew = codonTable.aa(codonNew);
								effectExpected = "CODON_INSERTION(" + aaOld + "/" + aaNew + ")";
							} else {
								codonNew = codonOld.substring(0, cdsCodonPos) + insPlus + codonOld.substring(cdsCodonPos);
								aaNew = codonTable.aa(codonNew);

								if (codonNew.startsWith(codonOld)) effectExpected = "CODON_INSERTION(" + aaOld + "/" + aaNew + ")";
								else effectExpected = "CODON_CHANGE_PLUS_CODON_INSERTION(" + aaOld + "/" + aaNew + ")";
							}

							if ((cdsCodonNum == 0) && codonTable.isStartFirst(codonOld) && !codonTable.isStartFirst(codonNew)) effectExpected = "START_LOST(" + aaOld + "/" + aaNew + ")";
							else if ((aaOld.indexOf('*') >= 0) && (aaNew.indexOf('*') < 0)) {
								effectExpected = "STOP_LOST(" + aaOld + "/" + aaNew + ")";
							} else if ((aaNew.indexOf('*') >= 0) && (aaOld.indexOf('*') < 0)) effectExpected = "STOP_GAINED(" + aaOld + "/" + aaNew + ")";
						}

						// Calculate effects
						List<ChangeEffect> effects = snpEffectPredictor.seqChangeEffect(seqChange);

						// There should be only one effect
						Assert.assertEquals(1, effects.size());

						// Show
						ChangeEffect effect = effects.get(0);
						String effStr = effect.effect(true, true, true);
						if (debug) System.out.println("\tPos: " + pos //
								+ "\tCDS base num: " + cdsBaseNum + " [" + cdsCodonNum + ":" + cdsCodonPos + "]" //
								+ "\t" + seqChange + "\tstrand" + (seqChange.getStrand() >= 0 ? "+" : "-") //
								+ "\tCodon: " + codonOld + " -> " + codonNew //
								+ "\tAA: " + aaOld + " -> " + aaNew //
								+ "\tEffect: " + effStr //
								+ "\tEffect expected: " + effectExpected //
						);

						// Check effect
						Assert.assertEquals(effectExpected, effStr);
					}
				}
			}
		}
	}
}
