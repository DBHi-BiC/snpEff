package ca.mcgill.mcb.pcingola.spliceSites;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Random;

import ca.mcgill.mcb.pcingola.fileIterator.FastaFileIterator;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Intron;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Markers;
import ca.mcgill.mcb.pcingola.interval.SpliceSiteBranch;
import ca.mcgill.mcb.pcingola.interval.SpliceSiteBranchU12;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.motif.MotifLogo;
import ca.mcgill.mcb.pcingola.motif.Pwm;
import ca.mcgill.mcb.pcingola.probablility.FisherExactTest;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.commandLine.SnpEff;
import ca.mcgill.mcb.pcingola.stats.CountByType;
import ca.mcgill.mcb.pcingola.stats.IntStats;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.mcb.pcingola.util.Tuple;

/**
 * Analyze sequences from splice sites
 * 
 * @author pcingola
 */
public class SpliceAnalysis extends SnpEff {

	/**
	 * A set of PWMs
	 * @author pablocingolani
	 *
	 */
	class PwmSet implements Comparable<PwmSet> {
		String name;
		Pwm pwmAcc, pwmDonor;
		CountByType countMotif;
		CountByType countExonTypes;
		IntStats lenStats;
		int motifMatchedBases = 0, motifMatchedStr = 0;
		int updates = 0;
		int countU12 = 0;

		public PwmSet(String name) {
			this.name = name;
			pwmAcc = new Pwm(2 * SpliceTypes.MAX_SPLICE_SIZE + 1);
			pwmDonor = new Pwm(2 * SpliceTypes.MAX_SPLICE_SIZE + 1);
			lenStats = new IntStats();
			countMotif = new CountByType();
			countExonTypes = new CountByType();
		}

		@Override
		public int compareTo(PwmSet ps) {
			int diff = ps.updates - updates;
			if (diff != 0) return diff;
			return name.compareTo(ps.name);
		}

		/**
		 * Count U12 ratio = Observed / expected
		 * @return
		 */
		public double countU12ObsExp() {
			double expected = updates * (1.0 - THRESHOLD_BRANCH_U12_PERCENTILE);
			double oe = countU12 / expected;
			return oe;
		}

		void incExonTypes(String exonTypes) {
			countExonTypes.inc(exonTypes);
		}

		void incU12() {
			countU12++; // Update count
		}

		void len(int len) {
			lenStats.sample(len);
		}

		String pExonTypes() {
			StringBuilder out = new StringBuilder();

			for (String type : countExonTypes.getTypeList())
				out.append(pExonTypes(type));

			return out.toString();
		}

		String pExonTypes(String category) {
			int countBlackDrawn = 0;
			for (String type : countExonTypes.getTypeList())
				if (!type.equals(category)) countBlackDrawn += countExonTypes.get(type);
			int countWhiteDrawn = (int) countExonTypes.get(category);

			// Get 'all' counts
			PwmSet pwmSet = getPwmSet(" ALL");
			CountByType countExonTypesAll = pwmSet.countExonTypes;
			int countBlack = 0;
			for (String type : countExonTypesAll.getTypeList())
				if (!type.equals(category)) countBlack += countExonTypesAll.get(type);
			int countWhite = (int) countExonTypesAll.get(category);

			String out = "";

			double pDown = FisherExactTest.get().fisherExactTestDown(countWhiteDrawn, countBlack + countWhite, countWhite, countBlackDrawn + countWhiteDrawn);
			if (pDown < P_VALUE_THRESHOLD) out += "p-value Down (" + category + ") : " + pDown + "\n";

			double pUp = FisherExactTest.get().fisherExactTestUp(countWhiteDrawn, countBlack + countWhite, countWhite, countBlackDrawn + countWhiteDrawn);
			if (pUp < P_VALUE_THRESHOLD) out += "p-value Up   (" + category + ") : " + pUp + "\n";

			return out;
		}

		@Override
		public String toString() {
			StringBuilder out = new StringBuilder();

			//out.append("<tr>\n");
			out.append("\t<td> <b>" + name + "</b> </td>\n");
			out.append("\t<td> " + updates + "</td>\n");

			// Donor motif
			MotifLogo mlDonor = new MotifLogo(pwmDonor);
			out.append("\t<td>\n");
			out.append(mlDonor.toStringHtml(HTML_WIDTH, HTML_HEIGHT));
			out.append("\t</td>\n");

			// U12 count
			double oe = countU12ObsExp();

			// U12 Colors
			String bg = "ffffff";
			if (oe > 5) bg = "ff0000";
			else if (oe > 2) bg = "ff8888";
			else if (oe > 1.2) bg = "ffcccc";
			out.append(String.format("\t<td bgcolor=%s> <center> %d (%1.2f)" + " </center> </td>\n", bg, countU12, oe));

			// Acceptor motif
			MotifLogo mlAcc = new MotifLogo(pwmAcc);
			out.append("\t<td>\n");
			out.append(mlAcc.toStringHtml(HTML_WIDTH, HTML_HEIGHT));
			out.append("\t</td>\n");

			// Intron length stats
			out.append("\t<td> <pre>\n");
			out.append(lenStats.toString());
			out.append("\t</pre></td>\n");

			// Count exon types
			out.append("\t<td> <pre>\n");
			out.append(countExonTypes);
			out.append("\t</pre></td>\n");

			// p-Values
			out.append("\t<td> <pre>\n");
			out.append(pExonTypes());
			out.append("\t</pre></td>\n");

			// out.append("</tr>\n");

			return out.toString();
		}

		public void update(String accStr, String donorStr) {
			updates++;
			if (accStr != null) pwmAcc.updateCounts(accStr);
			if (donorStr != null) pwmDonor.updateCounts(donorStr);

		}
	}

	public static double P_VALUE_THRESHOLD = 0.001;
	public static int SIZE_CONSENSUS_DONOR = 2;
	public static int SIZE_CONSENSUS_ACCEPTOR = 2;
	public static final double THRESHOLD_ENTROPY = 0.05;
	public static final int THRESHOLD_COUNT = 100;
	public static final double THRESHOLD_P = 0.95;
	public static final double THRESHOLD_BRANCH_U12_PERCENTILE = 0.95;
	public static final double THRESHOLD_U12_OBSERVED_EXPECTED = 5;
	public static int HTML_WIDTH = 20;
	public static int HTML_HEIGHT = 100;

	String outputDir = ".";
	String genomeVer;
	String genomeFasta;
	StringBuilder out = new StringBuilder();
	Config config;
	SpliceTypes spliceTypes;
	ArrayList<String> geneList = new ArrayList<String>();
	HashMap<String, PwmSet> pwmSetsByName = new HashMap<String, PwmSet>();
	HashMap<String, PwmSet> pwmSetsExonTypeByName = new HashMap<String, PwmSet>();
	HashMap<String, Intron> intronsByStr = new HashMap<String, Intron>();
	double thresholdPDonor;
	double thresholdEntropyDonor;
	double thresholdPAcc;
	double thresholdEntropyAcc;
	double thresholdU12Score;
	int countIntrons = 0;
	Random random = new Random();

	public SpliceAnalysis() {
		super();
	}

	/**
	 * Count how many entries that have both 'donor' and 'acceptor' 
	 * @param donor
	 * @param acceptor
	 * @return
	 */
	int countDonorAcc(String donor, String acceptor) {
		int count = 0;
		for (int i = 0; i < spliceTypes.getDonorAccPairSize(); i++) {
			String d = spliceTypes.getDonor(i);
			String a = spliceTypes.getAcceptor(i);

			if (d.startsWith(donor) && a.endsWith(acceptor)) count++;
		}
		return count;
	}

	/**
	 * Find an probability threshold using THRESHOLD_P quantile
	 * @param tree
	 * @return
	 */
	double findEntropyThreshold(AcgtTree tree) {
		List<Double> values = tree.entropyAll(THRESHOLD_COUNT);
		Collections.sort(values);
		int index = (int) (values.size() * THRESHOLD_ENTROPY);
		return values.get(index);
	}

	/**
	 * Find an probability threshold using THRESHOLD_P quantile
	 * @param tree
	 * @return
	 */
	double findPthreshold(AcgtTree tree) {
		List<Double> values = tree.pAll(THRESHOLD_COUNT);
		Collections.sort(values);
		int index = (int) (values.size() * THRESHOLD_P);
		return values.get(index);
	}

	PwmSet getPwmSet(String key) {
		PwmSet ps = pwmSetsByName.get(key);
		if (ps == null) {
			ps = new PwmSet(key);
			pwmSetsByName.put(key, ps);
		}
		return ps;
	}

	PwmSet getPwmSetExonType(String key) {
		PwmSet ps = pwmSetsExonTypeByName.get(key);
		if (ps == null) {
			ps = new PwmSet(key);
			pwmSetsExonTypeByName.put(key, ps);
		}
		return ps;
	}

	/**
	 * Initialize
	 */
	void init() {
		if (verbose) Timer.showStdErr("Initializing");
		config = new Config(genomeVer);
		genomeFasta = config.getFileNameGenomeFasta();
		if (genomeFasta == null) throw new RuntimeException("Cannot find reference genome: " + config.getFileListGenomeFasta());

		outputDir = config.getDirData() + "/spliceSites";

		// Load data
		load();
	}

	/**
	 * Lad data from files
	 */
	void load() {
		if (verbose) Timer.showStdErr("Loading: " + genomeVer);
		config.loadSnpEffectPredictor();
		if (verbose) Timer.showStdErr("done");
	}

	/**
	 * Show and append an output line
	 * @param line
	 */
	void out(Object o) {
		String s = o.toString();
		out.append(s + "\n");
		//System.out.println(s);
	}

	@Override
	public void parseArgs(String[] args) {
		if (args.length == 0) usage(null);
		int idx = 0;

		if (args[idx].equals("-v")) {
			verbose = true;
			idx++;
		}

		if ((args.length - idx) != 1) usage("Missing genome");
		genomeVer = args[idx];
	}

	@Override
	public boolean run() {
		init();

		spliceTypes = new SpliceTypes(config);
		spliceTypes.setVerbose(verbose);

		// Splice site conservation
		spliceTypes.analyzeAndCreate();
		thresholdU12Score = spliceTypes.branchU12Threshold(THRESHOLD_BRANCH_U12_PERCENTILE); // Find U12 branch points
		spliceTypes.createSpliceFasta(outputDir); // Create fasta files for splice sites

		// Splice site PWM analysis, create BranchPoints U12
		splicePwmAnalysis();

		//---
		// Save 
		//---
		String outputFile = outputDir + "/" + this.getClass().getSimpleName() + "_" + genomeVer + ".html";
		if (verbose) Timer.showStdErr("Saving output to: " + outputFile);
		Gpr.toFile(outputFile, out);

		if (verbose) Timer.showStdErr("Saving database to file: " + config.getFileSnpEffectPredictor());
		config.getSnpEffectPredictor().save(config);
		if (verbose) Timer.showStdErr("Done.");

		if (verbose) Timer.showStdErr("Finished!");
		return true;
	}

	/**
	 * Run PWM analysis 
	 */
	void splicePwmAnalysis() {
		if (verbose) Timer.showStdErr("Splice analysis (PWM). Reading fasta file: " + genomeFasta);

		//---
		// Main analysis: Iterate over all chromosomes
		//---
		FastaFileIterator ffi = new FastaFileIterator(genomeFasta);
		out("<pre>\n");
		for (String chrSeq : ffi) {
			String chrName = Chromosome.simpleName(ffi.getName());
			splicePwmAnalysis(chrName, chrSeq);
		}
		out("</pre>\n");

		//---
		// Create U12 sites
		//---
		String bedFile = outputDir + "/" + this.getClass().getSimpleName() + "_" + genomeVer + "_introns_branchSitesU12.bed";
		if (verbose) Timer.showStdErr("Writing Introns and SpliceSitesBranchU12 file to BED file: '" + bedFile + "'");

		// ArrayList<Marker> markers = new ArrayList<Marker>();
		Markers markersBed = new Markers();
		for (String donorAcc : pwmSetsByName.keySet()) {
			PwmSet pwmSet = getPwmSet(donorAcc);

			// Is it over threshold?
			if ((pwmSet.updates >= THRESHOLD_COUNT) && (pwmSet.countU12ObsExp() > THRESHOLD_U12_OBSERVED_EXPECTED)) {
				List<SpliceSiteBranchU12> ssbu12sites = spliceTypes.getBranchU12(donorAcc);

				// Add sites to transcript
				for (SpliceSiteBranchU12 bu12 : ssbu12sites) {
					Transcript tr = (Transcript) bu12.getParent();
					tr.add(bu12); // Add branch site to transcript
					markersBed.add(bu12); // Add to bed file
					Gpr.debug("Adding BranchU12 '" + bu12 + "' to transcript " + tr.getId() + "\tDonor-acceptor pair: " + donorAcc + "\tObs/Expected: " + getPwmSet(donorAcc).countU12ObsExp());
				}
			}
		}

		// Add introns
		markersBed.addAll(intronsByStr.values());
		markersBed.sort(false, false);

		// Create BED file
		StringBuilder sb = new StringBuilder();
		for (Marker i : markersBed) {
			sb.append(i.getChromosomeName() + "\t" + (i.getStart() + 1) + "\t" + (i.getEnd() + 1) + "\t" + (i instanceof SpliceSiteBranch ? i.getType().toString() : i.getId()) + "\n");
		}
		Gpr.toFile(bedFile, sb);

		//---
		// Show results
		//---
		if (verbose) Timer.showStdErr("Filter out low count splice sites. Exons: " + countIntrons + "\tThreshold: " + THRESHOLD_COUNT);
		ArrayList<PwmSet> pwmsets = new ArrayList<PwmSet>();
		pwmsets.addAll(pwmSetsByName.values());
		Collections.sort(pwmsets);
		out("<p><center><h3>Analysis by Donnor-Acceptor type</h3></center><p><table border=1>\n");
		out("<p><b>U12 PWM score threshold:</b> " + thresholdU12Score + "<p>\n");
		out("<tr> <th> Rank </th> <th> Donor-Acceptor </th>  <th> Count </th>  <th> Donor Motif </th> <th> U12 matches (Observed / Expected) </th> <th> Acceptor Motif </th> <th> Intron length </th> <th> Intron Type Count </th> <th> Intron Type p-values </th> </tr>\n");
		int count = 0;
		for (PwmSet pwmset : pwmsets)
			if (pwmset.updates >= THRESHOLD_COUNT) out("<tr> <td> " + (count++) + " </td> " + pwmset + "</tr>\n");
		out("</table>\n");

		// Show PwmSetsExonType
		pwmsets = new ArrayList<PwmSet>();
		pwmsets.addAll(pwmSetsExonTypeByName.values());
		Collections.sort(pwmsets);
		out("<p><hr><p><center><h3>Analysis by Exon-Exon types</h3></center><p><table border=1>\n");
		out("<tr> <th> Rank </th> <th> Exon_Type --- Exon_Type</th>  <th> Count </th>  <th> Donor Motif </th> <th> U12 matches (Observed / Expected) </th> <th> Acceptor Motif </th> <th> Intron length </th> <th> Intron Type Count </th> <th> Intron Type p-values </th> </tr>\n");
		count = 0;
		for (PwmSet pwmset : pwmsets)
			if (pwmset.updates >= THRESHOLD_COUNT) out("<tr> <td> " + (count++) + " </td> " + pwmset + "</tr>\n");
		out("</table>\n");
	}

	/**
	 * Run PWM analysis for one chromosome
	 * @param chrName
	 * @param chrSeq
	 */
	void splicePwmAnalysis(String chrName, String chrSeq) {
		int countEx = 0;
		HashSet<String> done = new HashSet<String>();

		//---
		// Find all exons in this chromosome
		//---
		for (Gene gene : config.getGenome().getGenes()) {
			if (gene.getChromosomeName().equals(chrName)) { // Same chromosome
				for (Transcript tr : gene) {
					Exon exPrev = null;
					for (Exon ex : tr.sortedStrand()) {
						countEx++;

						if (exPrev != null) { // Not for first exon (it has no 'previous' intron)
							int start, end;
							if (tr.isStrandPlus()) {
								start = exPrev.getEnd();
								end = ex.getStart();
							} else {
								start = ex.getEnd();
								end = exPrev.getStart();
							}

							// Get exon splice type
							String exPrevType = exPrev != null ? exPrev.getSpliceType().toString() : "";
							String exType = ex != null ? ex.getSpliceType().toString() : "";
							String exonTypes = exPrevType + "-" + exType;

							// Do not analyze this Intron if it was already analyzed
							String key = chrName + ":" + start + "-" + end;
							if (!done.contains(key)) {
								updatePwm(tr, chrSeq, start, end, exonTypes);
								done.add(key);

								// Create BED file
								Intron intron = new Intron(tr, start, end, 1, exonTypes);
								intronsByStr.put(intron.toString(), intron); // We use a hash to avoid adding the same intron multiple times
							}
						}

						exPrev = ex;
					}
				}
			}
		}

		if (verbose) Timer.showStdErr("\tChromosome: " + chrName + "\tGenes: " + config.getGenome().getGenes().size() + "\tExons: " + countEx);
	}

	/**
	 * Update PWM
	 * @param tr
	 * @param intronStart
	 * @param intronEnd
	 */
	void updatePwm(Transcript tr, String chrSeq, int intronStart, int intronEnd, String exonTypes) {
		// We don't update if the intron is too short
		int len = intronEnd - intronStart;
		if (len < (2 * SpliceTypes.MAX_SPLICE_SIZE)) return;

		String donorStr = spliceTypes.seqDonor(tr, chrSeq, intronStart, intronEnd);
		String accStr = spliceTypes.seqAcceptor(tr, chrSeq, intronStart, intronEnd);
		String intronSeqDonor = donorStr.substring(SpliceTypes.MAX_SPLICE_SIZE + 1);
		String intronSeqAcc = accStr.substring(0, SpliceTypes.MAX_SPLICE_SIZE);

		countIntrons++;

		//---
		// Group by donor type
		//---

		// Donor consensus ('GT' or 'GC'?)
		String donorConsensus = donorStr.substring(SpliceTypes.MAX_SPLICE_SIZE + 1, SpliceTypes.MAX_SPLICE_SIZE + 1 + SIZE_CONSENSUS_DONOR);
		if (donorConsensus.indexOf('N') >= 0) return; // Ignore if there is an 'N'

		// Use long consensus? U12
		String accConsensus = accStr.substring(SpliceTypes.MAX_SPLICE_SIZE - SIZE_CONSENSUS_ACCEPTOR, SpliceTypes.MAX_SPLICE_SIZE);
		if (donorConsensus.indexOf('N') >= 0) return; // Ignore if there is an 'N'

		int maxLenDa = 0;
		for (int i = 0; i < spliceTypes.getDonorAccPairSize(); i++) {
			String don = spliceTypes.getDonor(i);
			String ac = spliceTypes.getAcceptor(i);
			if (intronSeqDonor.startsWith(don) && intronSeqAcc.endsWith(ac)) {
				int lenda = don.length() + ac.length();
				if (lenda > maxLenDa) {
					maxLenDa = lenda;
					donorConsensus = don;
					accConsensus = ac;
				}
			}
		}
		String consensus = donorConsensus + "_" + accConsensus;

		//---
		// Branch U12 score
		//---
		Tuple<Double, Integer> bestU12 = spliceTypes.addBestU12Score(tr, chrSeq, consensus, intronStart, intronEnd);
		double bestU12score = bestU12.first;

		//---
		// Update PWM
		//---
		PwmSet pwmSet = getPwmSet(consensus);
		pwmSet.update(accStr, donorStr);
		pwmSet.len(len);
		pwmSet.incExonTypes(exonTypes);
		if (bestU12score >= thresholdU12Score) pwmSet.incU12();

		// Update total counts
		pwmSet = getPwmSet(" ALL");
		pwmSet.update(accStr, donorStr);
		pwmSet.incExonTypes(exonTypes);
		pwmSet.len(len);

		//---
		// Update PWM for exon type
		//---
		pwmSet = getPwmSetExonType(exonTypes);
		pwmSet.update(accStr, donorStr);
		pwmSet.len(len);
		if (bestU12score >= thresholdU12Score) pwmSet.incU12();
	}

	@Override
	public void usage(String message) {
		if (message != null) System.err.println("Error: " + message + "\n");
		System.err.println("Usage: snpEff genome_version");
		System.exit(-1);
	}

}
