package ca.mcgill.mcb.pcingola.snpEffect.testCases;

import java.util.Random;

import junit.framework.Assert;
import junit.framework.TestCase;
import ca.mcgill.mcb.pcingola.codons.CodonTable;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.interval.Variant;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffects;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.snpEffect.factory.SnpEffPredictorFactoryRand;
import ca.mcgill.mcb.pcingola.util.GprSeq;

/**
 * Test random SNP changes 
 * 
 * @author pcingola
 */
public class TestCasesSnp extends TestCase {

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

	public TestCasesSnp() {
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
		SnpEffPredictorFactoryRand sepf = new SnpEffPredictorFactoryRand(config, rand, maxGeneLen, maxTranscripts, maxExons);

		// Create predictor
		snpEffectPredictor = sepf.create();
		config.setSnpEffectPredictor(snpEffectPredictor);

		config.getSnpEffectPredictor().setSpliceRegionExonSize(0);
		config.getSnpEffectPredictor().setSpliceRegionIntronMin(0);
		config.getSnpEffectPredictor().setSpliceRegionIntronMax(0);

		// Chromosome sequence
		chromoSequence = sepf.getChromoSequence();
		chromoBases = chromoSequence.toCharArray();

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
		//	- Change each base in the exon
		//	- Calculate effect
		for (int i = 0; i < N; i++) {
			initSnpEffPredictor();
			if (debug) System.out.println("SNP Test iteration: " + i + "\n" + transcript);
			else System.out.println("SNP Test iteration: " + i + "\t" + (transcript.getStrand() >= 0 ? "+" : "-") + "\t" + transcript.cds());

			int cdsBaseNum = 0;

			// For each exon...
			for (Exon exon : transcript.sortedStrand()) {
				int step = exon.isStrandPlus() ? 1 : -1;
				int beg = exon.isStrandPlus() ? exon.getStart() : exon.getEnd();

				// For each base in this exon...
				for (int pos = beg; (pos >= exon.getStart()) && (pos <= exon.getEnd()); pos += step, cdsBaseNum++) {
					// Reference base
					char refBase = chromoBases[pos]; // exon.basesAt(pos - exon.getStart(), 1).charAt(0);
					refBase = Character.toUpperCase(refBase);
					// Codon number
					int cdsCodonNum = cdsBaseNum / 3;
					int cdsCodonPos = cdsBaseNum % 3;

					int minCodonPos = cdsCodonNum * 3;
					int maxCodonPos = minCodonPos + 3;
					if (maxCodonPos < transcript.cds().length()) {
						String codon = transcript.cds().substring(minCodonPos, maxCodonPos);
						codon = codon.toUpperCase();
						String aa = codonTable.aa(codon);

						// Get a random base different from 'refBase'
						char snp = refBase;
						while (snp == refBase) {
							snp = Character.toUpperCase(GprSeq.randBase(rand));
						}

						// Codon change
						String newCodon = codon.substring(0, cdsCodonPos) + snp + codon.substring(cdsCodonPos + 1);
						String newAa = codonTable.aa(newCodon);
						String effectExpected = "";

						// Effect
						if (newAa.equals(aa)) {
							if ((cdsCodonNum == 0) && (codonTable.isStart(codon))) {
								if (codonTable.isStart(newCodon)) effectExpected = "SYNONYMOUS_START(" + aa + ")";
								else effectExpected = "START_LOST(" + aa + ")";
							} else if (aa.equals("*")) effectExpected = "SYNONYMOUS_STOP(" + aa + ")";
							else effectExpected = "SYNONYMOUS_CODING(" + aa + ")";
						} else {
							if ((cdsCodonNum == 0) && (codonTable.isStart(codon))) {
								if (codonTable.isStart(newCodon)) effectExpected = "NON_SYNONYMOUS_START(" + aa + "/" + newAa + ")";
								else effectExpected = "START_LOST(" + aa + "/" + newAa + ")";
							} else if (codonTable.isStop(codon)) effectExpected = "STOP_LOST(" + aa + "/" + newAa + ")";
							else if (codonTable.isStop(newCodon)) effectExpected = "STOP_GAINED(" + aa + "/" + newAa + ")";
							else effectExpected = "NON_SYNONYMOUS_CODING(" + aa + "/" + newAa + ")";
						}

						// Create a SeqChange
						int seqChangeStrand = 1;
						if (exon.isStrandMinus()) refBase = GprSeq.wc(refBase);
						if (seqChangeStrand == -exon.getStrand()) {
							snp = GprSeq.wc(snp);
							refBase = GprSeq.wc(refBase);
						}
						Variant seqChange = new Variant(chromosome, pos, refBase + "", snp + "", "");

						if (!seqChange.isVariant()) effectExpected = "EXON";

						// Calculate effects
						ChangeEffects effects = snpEffectPredictor.seqChangeEffect(seqChange);

						// There should be only one effect
						if (debug) System.out.println(effects);
						Assert.assertEquals(true, effects.size() <= 1);

						// Show
						if (effects.size() == 1) {
							ChangeEffect effect = effects.get();
							String effStr = effect.effect(true, true, true, false);
							if (debug) System.out.println("\tPos: " + pos //
									+ "\tCDS base num: " + cdsBaseNum + " [" + cdsCodonNum + ":" + cdsCodonPos + "]" //
									+ "\t" + seqChange + (seqChange.getStrand() >= 0 ? "+" : "-") //
									+ "\tCodon: " + codon + " -> " + newCodon //
									+ "\tAA: " + aa + " -> " + newAa //
									+ "\tEffect: " + effStr);

							// Check effect
							Assert.assertEquals(effectExpected, effStr);
						}
					}
				}
			}
		}
	}
}
