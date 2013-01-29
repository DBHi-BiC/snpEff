package ca.mcgill.mcb.pcingola.snpEffect.commandLine;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import akka.actor.Actor;
import akka.actor.Props;
import akka.actor.UntypedActorFactory;
import ca.mcgill.mcb.pcingola.akka.vcf.VcfWorkQueue;
import ca.mcgill.mcb.pcingola.fileIterator.SeqChangeBedFileIterator;
import ca.mcgill.mcb.pcingola.fileIterator.SeqChangeFileIterator;
import ca.mcgill.mcb.pcingola.fileIterator.SeqChangeFilePileUp;
import ca.mcgill.mcb.pcingola.fileIterator.SeqChangeFileTxt;
import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.filter.ChangeEffectFilter;
import ca.mcgill.mcb.pcingola.filter.SeqChangeFilter;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Custom;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.interval.Regulation;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.codonChange.CodonChange;
import ca.mcgill.mcb.pcingola.interval.tree.IntervalForest;
import ca.mcgill.mcb.pcingola.outputFormatter.BedAnnotationOutputFormatter;
import ca.mcgill.mcb.pcingola.outputFormatter.BedOutputFormatter;
import ca.mcgill.mcb.pcingola.outputFormatter.OutputFormatter;
import ca.mcgill.mcb.pcingola.outputFormatter.TxtOutputFormatter;
import ca.mcgill.mcb.pcingola.outputFormatter.VcfOutputFormatter;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectImpact;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.snpEffect.commandLine.eff.MasterEff;
import ca.mcgill.mcb.pcingola.stats.ChangeEffectResutStats;
import ca.mcgill.mcb.pcingola.stats.SeqChangeStats;
import ca.mcgill.mcb.pcingola.stats.VcfStats;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.mcb.pcingola.vcf.VcfEffect;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import freemarker.template.Configuration;
import freemarker.template.DefaultObjectWrapper;
import freemarker.template.Template;
import freemarker.template.TemplateException;

/**
 * Command line program: Predict changes
 * 
 * @author Pablo Cingolani
 */
public class SnpEffCmdEff extends SnpEff {

	public static final String SUMMARY_TEMPLATE = "snpEff_summary.ftl"; // Summary template file name
	public static final String SUMMARY_GENES_TEMPLATE = "snpEff_genes.ftl"; // Genes template file name

	boolean canonical = false; // Use only canonical transcripts
	boolean supressOutput = false; // Only used for debugging purposes 
	boolean createSummary = true; // Do not create summary output file 
	boolean useLocalTemplate = false; // Use template from 'local' file instead of 'jar' (this is only used for development and debugging)
	Boolean treatAllAsProteinCoding = null; // Only use coding genes. Default is 'null' which means 'auto'
	boolean chromoPlots = true; // Create methylation by chromosome plots?
	boolean onlyRegulation = false; // Only build regulation tracks
	boolean lossOfFunction = false; // Create loss of function LOF tag?
	int upDownStreamLength = SnpEffectPredictor.DEFAULT_UP_DOWN_LENGTH; // Upstream & downstream interval length
	int totalErrs = 0;
	long countInputLines = 0, countVariants = 0, countEffects = 0, countVariantsFilteredOut = 0;
	String chrStr = "";
	String inputFile = "-"; // Input file
	String summaryFile; // Summary output file
	String summaryGenesFile; // Gene table file
	String onlyTranscriptsFile = null; // Only use the transcripts in this file (Format: One transcript ID per line)
	SeqChangeFilter seqChangeFilter; // Filter seqChanges (before prediction)
	InputFormat inputFormat = InputFormat.VCF; // Format use in input files
	OutputFormat outputFormat = OutputFormat.VCF; // Output format
	ChangeEffectFilter changeEffectResutFilter; // Filter prediction results
	ArrayList<String> filterIntervalFiles;// Files used for filter intervals
	IntervalForest filterIntervals; // Filter only seqChanges that match these intervals
	ArrayList<String> customIntervalFiles; // Custom interval files
	SeqChangeStats seqChangeStats;
	ChangeEffectResutStats changeEffectResutStats;
	VcfStats vcfStats;
	HashSet<String> regulationTracks = new HashSet<String>();
	List<VcfEntry> vcfEntriesDebug = null; // Use for debugging or testing (in some test-cases)

	public SnpEffCmdEff() {
		super();
		upDownStreamLength = SnpEffectPredictor.DEFAULT_UP_DOWN_LENGTH; // Upstream & downstream interval length
		chrStr = ""; // Default: Don't show 'chr' before chromosome
		inputFile = ""; // seqChange input file
		seqChangeFilter = new SeqChangeFilter(); // Filter seqChanges (before prediction)
		changeEffectResutFilter = new ChangeEffectFilter(); // Filter prediction results
		filterIntervalFiles = new ArrayList<String>(); // Files used for filter intervals
		filterIntervals = new IntervalForest(); // Filter only seqChanges that match these intervals
		customIntervalFiles = new ArrayList<String>(); // Custom interval files
		summaryFile = DEFAULT_SUMMARY_FILE;
		summaryGenesFile = DEFAULT_SUMMARY_GENES_FILE;
	}

	public ChangeEffectResutStats getChangeEffectResutStats() {
		return changeEffectResutStats;
	}

	public SeqChangeStats getSeqChangeStats() {
		return seqChangeStats;
	}

	/**
	 * Iterate on all inputs and calculate effects.
	 * Note: This is used for all input formats except VCF, which has a different iteration modality
	 * 
	 * @param outputFormatter
	 */
	void iterateSeqChange(OutputFormatter outputFormatter) {
		SnpEffectPredictor snpEffectPredictor = config.getSnpEffectPredictor();

		// Create an input file iterator
		SeqChangeFileIterator seqChangeFileIterator;
		if (inputFormat == InputFormat.PILEUP) seqChangeFileIterator = new SeqChangeFilePileUp(inputFile, config.getGenome(), inOffset);
		else if (inputFormat == InputFormat.BED) seqChangeFileIterator = new SeqChangeBedFileIterator(inputFile, config.getGenome(), inOffset);
		else if (inputFormat == InputFormat.TXT) seqChangeFileIterator = new SeqChangeFileTxt(inputFile, config.getGenome(), inOffset);
		else throw new RuntimeException("Cannot create SeqChange file iterator on input format '" + inputFormat + "'");

		//---
		// Iterate over input file
		//---
		for (SeqChange seqChange : seqChangeFileIterator) {
			try {
				countInputLines++;

				countVariants += seqChange.getChangeOptionCount();
				if (verbose && (countVariants % 100000 == 0)) Timer.showStdErr("\t" + countVariants + " variants");

				// Does it pass the filter? => Analyze
				if ((seqChangeFilter == null) || seqChangeFilter.filter(seqChange)) {

					// Skip if there are filter intervals and they are not matched 
					if ((filterIntervals != null) && (filterIntervals.stab(seqChange).size() <= 0)) continue;

					// Perform basic statistics about this seqChange
					if (createSummary) seqChangeStats.sample(seqChange);

					// Calculate effects
					List<ChangeEffect> changeEffects = snpEffectPredictor.seqChangeEffect(seqChange);

					// Create new 'section'
					outputFormatter.startSection(seqChange);

					// Show results
					for (ChangeEffect changeEffect : changeEffects) {
						changeEffectResutStats.sample(changeEffect); // Perform basic statistics about this result
						outputFormatter.add(changeEffect);
						countEffects++;
					}

					// Finish up this section
					outputFormatter.printSection(seqChange);

				} else countVariantsFilteredOut += seqChange.getChangeOptionCount();
			} catch (Throwable t) {
				totalErrs++;
				error(t, "Error while processing variant (line " + seqChangeFileIterator.getLineNum() + ") :\n\t" + seqChange + "\n" + t);
			}
		}

		// Close file iterator (not really needed, but just in case)
		seqChangeFileIterator.close();
	}

	/**
	 * Iterate on all inputs (VCF) and calculate effects.
	 * Note: This is used only on input format VCF, which has a different iteration modality
	 * 
	 * @param outputFormatter
	 */
	void iterateVcf(OutputFormatter outputFormatter) {
		SnpEffectPredictor snpEffectPredictor = config.getSnpEffectPredictor();

		// Open VCF file
		VcfFileIterator vcfFile = new VcfFileIterator(inputFile, config.getGenome());
		vcfFile.setInOffset(inOffset); // May be there is a special inOffset (not likely to happen).

		for (VcfEntry vcfEntry : vcfFile) {
			try {
				countInputLines++;

				// Sample vcf entry
				if (createSummary) vcfStats.sample(vcfEntry);

				// Skip if there are filter intervals and they are not matched 
				if ((filterIntervals != null) && (filterIntervals.query(vcfEntry).isEmpty())) continue;

				// Create new 'section'
				outputFormatter.startSection(vcfEntry);

				for (SeqChange seqChange : vcfEntry.seqChanges()) {
					countVariants += seqChange.getChangeOptionCount();
					if (verbose && (countVariants % 100000 == 0)) Timer.showStdErr("\t" + countVariants + " variants");

					// Does it pass the filter? => Analyze
					if ((seqChangeFilter == null) || seqChangeFilter.filter(seqChange)) {
						// Perform basic statistics about this seqChange
						if (createSummary) seqChangeStats.sample(seqChange);

						// Calculate effects
						List<ChangeEffect> changeEffects = snpEffectPredictor.seqChangeEffect(seqChange);

						// Create new 'section'
						outputFormatter.startSection(seqChange);

						// Show results
						for (ChangeEffect changeEffect : changeEffects) {
							if (createSummary) changeEffectResutStats.sample(changeEffect); // Perform basic statistics about this result

							if (changeEffect.getEffectImpact() == EffectImpact.MODIFIER) {
							}

							outputFormatter.add(changeEffect);
							countEffects++;
						}

						// Finish up this section
						outputFormatter.printSection(seqChange);

					} else countVariantsFilteredOut += seqChange.getChangeOptionCount();
				}

				// Finish up this section
				outputFormatter.printSection(vcfEntry);

			} catch (Throwable t) {
				totalErrs++;
				error(t, "Error while processing VCF entry (line " + vcfFile.getLineNum() + ") :\n\t" + vcfEntry + "\n" + t);
			}
		}

		// Close file iterator (not really needed, but just in case)
		vcfFile.close();
	}

	/**
	 * Multi-threaded iteration on VCF inputs and calculates effects.
	 * Note: This is used only on input format VCF, which has a different iteration modality
	 * 
	 * @param outputFormatter
	 */
	void iterateVcfMulti(final OutputFormatter outputFormatter) {
		if (verbose) Timer.showStdErr("Running multi-threaded mode (numThreads=" + numWorkers + ").");

		outputFormatter.setShowHeader(false); // Master process takes care of the header (instead of outputFormatter). Otherwise you get the header printed one time per worker.

		// We need final variables for the inner class
		final SnpEffectPredictor snpEffectPredictor = config.getSnpEffectPredictor();
		final VcfOutputFormatter vcfOutForm = (VcfOutputFormatter) outputFormatter;
		final SnpEffCmdEff snpEffCmdEff = this;

		// Open VCF file
		VcfFileIterator vcfFile = new VcfFileIterator(inputFile, config.getGenome());
		vcfFile.setInOffset(inOffset); // May be there is a special inOffset (not likely to happen).

		// Master factory 
		Props props = new Props(new UntypedActorFactory() {

			private static final long serialVersionUID = 1L;

			@Override
			public Actor create() {
				MasterEff master = new MasterEff(numWorkers, snpEffCmdEff, snpEffectPredictor, outputFormatter, filterIntervals, seqChangeFilter);
				master.setAddHeader(vcfOutForm.getNewHeaderLines().toArray(new String[0]));
				return master;
			}
		});

		// Create and run queue
		int batchSize = 10;
		VcfWorkQueue vcfWorkQueue = new VcfWorkQueue(inputFile, batchSize, -1, props);
		vcfWorkQueue.run(true);
	}

	/**
	 * Parse command line arguments
	 * @param args
	 */
	@Override
	public void parseArgs(String[] args) {
		this.args = args;
		for (int i = 0; i < args.length; i++) {

			// Argument starts with '-'?
			if (args[i].startsWith("-")) {

				//---
				// Generic options
				//---
				if ((args[i].equals("-c") || args[i].equalsIgnoreCase("-config"))) {
					if ((i + 1) < args.length) configFile = args[++i];
					else usage("Option '-c' without config file argument");
				} else if (args[i].equals("-v") || args[i].equalsIgnoreCase("-verbose")) {
					verbose = true;
					quiet = false;
				} else if (args[i].equals("-q") || args[i].equalsIgnoreCase("-quiet")) {
					quiet = true;
					verbose = false;
				} else if (args[i].equals("-1")) inOffset = outOffset = 1;
				else if (args[i].equals("-0")) inOffset = outOffset = 0;
				else if (args[i].equals("-h") || args[i].equalsIgnoreCase("-help")) {
					usage(null);
					System.exit(0);
				} else if (args[i].equals("-t")) {
					multiThreaded = true;
					createSummary = false; // Implies '-noStats'
				}
				//---
				// Output options
				//---
				else if (args[i].equals("-o")) {
					// Output format
					if ((i + 1) < args.length) {
						String outFor = args[++i].toUpperCase();

						if (outFor.equals("TXT")) {
							outputFormat = OutputFormat.TXT;
							outOffset = 1; // Implies '-1' since TXT coordinates are one-based
						} else if (outFor.equals("VCF")) {
							outputFormat = OutputFormat.VCF;
							outOffset = 1; // Implies '-1' since VCF coordinates are one-based
						} else if (outFor.equals("GATK")) {
							outputFormat = OutputFormat.GATK;
							outOffset = 1; // Implies '-1' since VCF coordinates are one-based
						} else if (outFor.equals("BED")) {
							outputFormat = OutputFormat.BED;
							outOffset = 0; // Implies '0' since BED coordinates are zero-based
						} else if (outFor.equals("BEDANN")) {
							outputFormat = OutputFormat.BEDANN;
							outOffset = 0; // Implies '0' since BED coordinates are zero-based
						} else usage("Unknown output file format '" + outFor + "'");
					}
				} else if ((args[i].equals("-a") || args[i].equalsIgnoreCase("-around"))) {
					if ((i + 1) < args.length) CodonChange.SHOW_CODONS_AROUND_CHANGE = Gpr.parseIntSafe(args[++i]);
					else usage("Option '-i' without config interval_file argument");
				} else if ((args[i].equals("-s") || args[i].equalsIgnoreCase("-stats"))) {
					if ((i + 1) < args.length) {
						summaryFile = args[++i];
						String dir = Gpr.dirName(summaryFile);
						summaryGenesFile = (dir != null ? dir + "/" : "") + Gpr.baseName(summaryFile, ".html") + ".genes.txt";
					}
				} else if ((args[i].equals("-of") || args[i].equalsIgnoreCase("-outOffset"))) {
					if ((i + 1) < args.length) outOffset = Gpr.parseIntSafe(args[++i]);
				} else if (args[i].equalsIgnoreCase("-chr")) chrStr = args[++i];
				else if (args[i].equalsIgnoreCase("-useLocalTemplate")) useLocalTemplate = true; // Undocumented option (only used for development & debugging)
				else if (args[i].equalsIgnoreCase("-noOut")) supressOutput = true; // Undocumented option (only used for development & debugging)
				else if (args[i].equalsIgnoreCase("-noStats")) createSummary = false; // Do not create summary file. It can be much faster (e.g. when parsing VCF files with many samples)
				else if (args[i].equalsIgnoreCase("-noChromoPlots")) chromoPlots = false;
				//---
				// Annotation options
				//---
				else if (args[i].equalsIgnoreCase("-treatAllAsProteinCoding")) {
					if ((i + 1) < args.length) {
						i++;
						if (args[i].equalsIgnoreCase("auto")) treatAllAsProteinCoding = null;
						else treatAllAsProteinCoding = Gpr.parseBoolSafe(args[i]);
					}
				} else if (args[i].equalsIgnoreCase("-canon")) canonical = true; // Use canonical transcripts
				else if (args[i].equalsIgnoreCase("-lof")) lossOfFunction = true; // Add LOF tag
				else if (args[i].equalsIgnoreCase("-onlyTr")) {
					if ((i + 1) < args.length) onlyTranscriptsFile = args[++i]; // Only use the transcripts in this file
				}
				//---
				// Input options
				//---
				else if (args[i].equalsIgnoreCase("-interval")) {
					if ((i + 1) < args.length) customIntervalFiles.add(args[++i]);
					else usage("Option '-i' without config interval_file argument");
				} else if ((args[i].equals("-fi") || args[i].equalsIgnoreCase("-filterInterval"))) {
					if ((i + 1) < args.length) filterIntervalFiles.add(args[++i]);
					else usage("Option '-fi' without config filter_interval_file argument");
				} else if (args[i].equals("-i")) {
					// Input format
					if ((i + 1) < args.length) {
						String inFor = args[++i].toUpperCase();

						if (inFor.equals("TXT")) {
							inputFormat = InputFormat.TXT;
							inOffset = 1; // Implies '-1' since TXT coordinates are one-based
						} else if (inFor.equals("PILEUP")) {
							inputFormat = InputFormat.PILEUP;
							inOffset = 1; // Implies '-1' since PILEUP coordinates are one-based
						} else if (inFor.equals("VCF")) {
							inputFormat = InputFormat.VCF;
							inOffset = 1; // Implies '-1' since VCF coordinates are one-based
						} else if (inFor.equals("BED")) {
							inputFormat = InputFormat.BED;
							inOffset = 0; // Implies '-0' since Bed coordinates are zero-based
						} else usage("Unknown input file format '" + inFor + "'");
					}
				} else if ((args[i].equals("-if") || args[i].equalsIgnoreCase("-inOffset"))) {
					if ((i + 1) < args.length) inOffset = Gpr.parseIntSafe(args[++i]);
				}
				//---
				// Regulation options
				//---
				else if (args[i].equals("-onlyReg")) onlyRegulation = true;
				else if (args[i].equals("-reg")) {
					if ((i + 1) < args.length) regulationTracks.add(args[++i]); // Add this track to the list
				}
				//---
				// Filters
				//---
				else if ((args[i].equals("-minQ") || args[i].equalsIgnoreCase("-minQuality"))) {
					if ((i + 1) < args.length) seqChangeFilter.setMinQuality(Gpr.parseIntSafe(args[++i]));
				} else if ((args[i].equals("-maxQ") || args[i].equalsIgnoreCase("-maxQuality"))) {
					if ((i + 1) < args.length) seqChangeFilter.setMaxQuality(Gpr.parseIntSafe(args[++i]));
				} else if ((args[i].equals("-minC") || args[i].equalsIgnoreCase("-minCoverage"))) {
					if ((i + 1) < args.length) seqChangeFilter.setMinCoverage(Gpr.parseIntSafe(args[++i]));
				} else if ((args[i].equals("-maxC") || args[i].equalsIgnoreCase("-maxCoverage"))) {
					if ((i + 1) < args.length) seqChangeFilter.setMaxCoverage(Gpr.parseIntSafe(args[++i]));
				} else if ((args[i].equals("-ud") || args[i].equalsIgnoreCase("-upDownStreamLen"))) {
					if ((i + 1) < args.length) upDownStreamLength = Gpr.parseIntSafe(args[++i]);
				} else if (args[i].equals("-hom")) seqChangeFilter.setHeterozygous(false);
				else if (args[i].equals("-het")) seqChangeFilter.setHeterozygous(true);
				else if (args[i].equals("-snp")) seqChangeFilter.setChangeType(SeqChange.ChangeType.SNP);
				else if (args[i].equals("-mnp")) seqChangeFilter.setChangeType(SeqChange.ChangeType.MNP);
				else if (args[i].equals("-ins")) seqChangeFilter.setChangeType(SeqChange.ChangeType.INS);
				else if (args[i].equals("-del")) seqChangeFilter.setChangeType(SeqChange.ChangeType.DEL);
				else if (args[i].equalsIgnoreCase("-no-downstream")) changeEffectResutFilter.setDownstream(true);
				else if (args[i].equalsIgnoreCase("-no-upstream")) changeEffectResutFilter.setUpstream(true);
				else if (args[i].equalsIgnoreCase("-no-intergenic")) changeEffectResutFilter.setIntergenic(true);
				else if (args[i].equalsIgnoreCase("-no-intron")) changeEffectResutFilter.setIntron(true);
				else if (args[i].equalsIgnoreCase("-no-utr")) changeEffectResutFilter.setUtr(true);
				else if (args[i].equalsIgnoreCase("-no")) {
					String filterOut = "";
					if ((i + 1) < args.length) filterOut = args[++i];

					String filterOutArray[] = filterOut.split(",");
					for (String filterStr : filterOutArray) {
						if (filterStr.equalsIgnoreCase("downstream")) changeEffectResutFilter.setDownstream(true);
						else if (filterStr.equalsIgnoreCase("upstream")) changeEffectResutFilter.setUpstream(true);
						else if (filterStr.equalsIgnoreCase("intergenic")) changeEffectResutFilter.setIntergenic(true);
						else if (filterStr.equalsIgnoreCase("intron")) changeEffectResutFilter.setIntron(true);
						else if (filterStr.equalsIgnoreCase("utr")) changeEffectResutFilter.setUtr(true);
						else if (filterStr.equalsIgnoreCase("None")) ; // OK, nothing to do
						else usage("Unknown filter option '" + filterStr + "'");
					}
				} else usage("Unknow option '" + args[i] + "'");
			} else if (genomeVer.length() <= 0) genomeVer = args[i];
			else if (inputFile.length() <= 0) inputFile = args[i];
			else usage("Unknow parameter '" + args[i] + "'");
		}

		// Check: Do we have all required parameters?
		if (genomeVer.isEmpty()) usage("Missing genomer_version parameter");
		if (inputFile.isEmpty()) inputFile = "-"; // Use STDIN

		// Sanity checks
		if ((outputFormat == OutputFormat.VCF) && (inputFormat != InputFormat.VCF)) usage("Output in VCF format is only supported when the input is also in VCF format");
		if (multiThreaded && (outputFormat != OutputFormat.VCF)) usage("Multi-threaded option is only supported when when output is in VCF format");
		if (multiThreaded && createSummary) usage("Multi-threaded option should be used with 'noStats'.");
		if (lossOfFunction && (outputFormat != OutputFormat.VCF)) usage("Loss of function annotation is only supported when when output is in VCF format");
	}

	/**
	 * Read a custom interval file
	 * @param intFile
	 */
	int readCustomIntFile(String intFile) {
		String file = readFile(intFile);
		String lines[] = file.split("\n");
		int count = 0;
		for (int lineNum = 0; lineNum < lines.length; lineNum++) {
			Custom ci = new Custom(null, 0, 0, 0, "");
			ci.readTxt(lines[lineNum], lineNum + 1, config.getGenome(), inOffset);
			config.getSnpEffectPredictor().add(ci);
			count++;
		}
		return count;
	}

	/**
	 * Read a file after checking for some common error conditions
	 * @param fileName
	 * @return
	 */
	String readFile(String fileName) {
		File file = new File(fileName);
		if (!file.exists()) fatalError("No such file '" + fileName + "'");
		if (!file.canRead()) fatalError("Cannot open file '" + fileName + "'");
		return Gpr.readFile(fileName);
	}

	/**
	 * Read a filter custom interval file
	 * @param intFile
	 */
	int readFilterIntFile(String intFile) {
		String file = readFile(intFile);
		String lines[] = file.split("\n");
		int count = 0;
		for (int lineNum = 0; lineNum < lines.length; lineNum++) {
			Custom ci = new Custom(null, 0, 0, 0, "");
			ci.readTxt(lines[lineNum], lineNum + 1, config.getGenome(), inOffset);
			filterIntervals.add(ci);
			count++;
		}
		return count;
	}

	/**
	 * Read regulation track and update SnpEffectPredictor
	 * @param regTrack
	 */
	@SuppressWarnings("unchecked")
	void readRegulationTrack(String regTrack) {
		//---
		// Read file
		//---
		if (verbose) Timer.showStdErr("Reading regulation track '" + regTrack + "'");
		String regFile = config.getDirDataVersion() + "/regulation_" + regTrack + ".bin";
		ArrayList<Regulation> regulation = (ArrayList<Regulation>) Gpr.readFileSerializedGz(regFile);

		//---
		// Are all chromosomes available?
		//---
		Genome genome = config.getGenome();
		HashMap<String, Integer> chrs = new HashMap<String, Integer>();
		for (Regulation r : regulation) {
			String chr = r.getChromosomeName();
			int max = chrs.containsKey(chr) ? chrs.get(chr) : 0;
			max = Math.max(max, r.getEnd());
			chrs.put(chr, max);
		}

		// Add all chromos
		for (String chr : chrs.keySet())
			if (genome.getChromosome(chr) == null) genome.add(new Chromosome(genome, 0, chrs.get(chr), 1, chr));

		//---
		// Add all markers to predictor
		//---
		SnpEffectPredictor snpEffectPredictor = config.getSnpEffectPredictor();
		for (Regulation r : regulation)
			snpEffectPredictor.add(r);

	}

	@Override
	public HashMap<String, String> reportValues() {
		HashMap<String, String> report = super.reportValues();
		if (seqChangeStats != null) report.put("SeqChanges", seqChangeStats.getCount() + "");
		return report;
	}

	/**
	 * Run according to command line options
	 */
	@Override
	public boolean run() {
		run(false);
		return true;
	}

	/**
	 * Run according to command line options
	 */
	public List<VcfEntry> run(boolean createList) {
		//---
		// Run predictor
		//---
		// Nothing to filter out => don't waste time
		if (!changeEffectResutFilter.anythingSet()) changeEffectResutFilter = null;

		filterIntervals = null;

		// Read config file
		if (verbose) Timer.showStdErr("Reading configuration file '" + configFile + "'");
		config = new Config(genomeVer, configFile); // Read configuration
		if (verbose) Timer.showStdErr("done");

		// Read database (or create a new one)
		if (onlyRegulation) {
			// Create predictor
			config.setSnpEffectPredictor(new SnpEffectPredictor(config.getGenome()));
			config.setOnlyRegulation(true);
			config.setErrorOnMissingChromo(false); // A chromosome might be missing (e.g. no regulation tracks available for 'MT')
			config.setErrorChromoHit(false); // A chromosome's length might be smaller than the real (it's calculated using regulation features, not real chromo data)
		} else {
			// Read
			if (verbose) Timer.showStdErr("Reading database for genome version '" + genomeVer + "' from file '" + config.getFileSnpEffectPredictor() + "' (this might take a while)");
			config.loadSnpEffectPredictor(); // Read snpEffect predictor
			if (verbose) Timer.showStdErr("done");
		}

		// Check if we can open the input file (no need to check if it is STDIN)
		if (!inputFile.equals("-") && !new File(inputFile).canRead()) usage("Cannot open input file '" + inputFile + "'");

		// Set 'treatAllAsProteinCoding'
		if (treatAllAsProteinCoding != null) config.setTreatAllAsProteinCoding(treatAllAsProteinCoding);
		else {
			// treatAllAsProteinCoding was set to 'auto'
			// I.e.: Use 'true' if there is protein coding info, otherwise use false.
			boolean tapc = !config.getGenome().hasCodingInfo();
			if (verbose) Timer.showStdErr("Setting '-treatAllAsProteinCoding' to '" + tapc + "'");
			config.setTreatAllAsProteinCoding(tapc);
		}

		// Read custom interval files
		for (String intFile : customIntervalFiles) {
			if (verbose) Timer.showStdErr("Reading interval file '" + intFile + "'");
			int count = readCustomIntFile(intFile);
			if (verbose) Timer.showStdErr("done (" + count + " intervals loaded). ");
		}

		// Read filter interval files
		for (String filterIntFile : filterIntervalFiles) {
			if (filterIntervals == null) filterIntervals = new IntervalForest();
			if (verbose) Timer.showStdErr("Reading filter interval file '" + filterIntFile + "'");
			int count = readFilterIntFile(filterIntFile);
			if (verbose) Timer.showStdErr("done (" + count + " intervals loaded). ");
		}

		// Read regulation tracks
		for (String regTrack : regulationTracks)
			readRegulationTrack(regTrack);

		// Build interval forest for filter (if any)
		if (filterIntervals != null) {
			if (verbose) Timer.showStdErr("Building filter interval forest");
			filterIntervals.build();
			if (verbose) Timer.showStdErr("done.");
		}

		// Set upstream-downstream interval length
		config.getSnpEffectPredictor().setUpDownStreamLength(upDownStreamLength);

		// Filter canonical transcripts
		if (canonical) {
			if (verbose) Timer.showStdErr("Filtering out non-canonical transcripts.");
			config.getSnpEffectPredictor().removeNonCanonical();
			if (verbose) Timer.showStdErr("done.");
		}

		// Filter canonical transcripts
		if (onlyTranscriptsFile != null) {
			// Load file
			String onlyTr = Gpr.readFile(onlyTranscriptsFile);
			HashSet<String> trIds = new HashSet<String>();
			for (String trId : onlyTr.split("\n"))
				trIds.add(trId.trim());

			// Remove transcripts
			if (verbose) Timer.showStdErr("Filtering out transcripts in file '" + onlyTranscriptsFile + "'. Total " + trIds.size() + " transcript IDs.");
			int removed = config.getSnpEffectPredictor().keepTranscripts(trIds);
			if (verbose) Timer.showStdErr("Done: " + removed + " transcripts removed.");
		}

		// Build tree
		if (verbose) Timer.showStdErr("Building interval forest");
		config.getSnpEffectPredictor().buildForest();
		if (verbose) Timer.showStdErr("done.");

		// Show some genome stats. Chromosome names are shown, a lot of people has problems with the correct chromosome names.
		if (verbose) Timer.showStdErr("Genome stats :\n" + config.getGenome());

		// Store VCF results in a list?
		if (createList) vcfEntriesDebug = new ArrayList<VcfEntry>();

		// Predict
		if (verbose) Timer.showStdErr("Predicting variants");
		runAnalysis();
		if (verbose) Timer.showStdErr("done.");

		return vcfEntriesDebug;
	}

	/**
	 * Calculate the effect of variants and show results
	 * @param snpEffFile
	 */
	public void runAnalysis() {
		seqChangeStats = new SeqChangeStats(config.getGenome());
		changeEffectResutStats = new ChangeEffectResutStats(config.getGenome());
		vcfStats = new VcfStats();

		int totalErrs = 0;

		//---
		// Create output formatter
		//---
		OutputFormatter outputFormatter = null;
		switch (outputFormat) {
		case TXT:
			outputFormatter = new TxtOutputFormatter();
			break;
		case VCF:
			outputFormatter = new VcfOutputFormatter(vcfEntriesDebug);
			((VcfOutputFormatter) outputFormatter).setLossOfFunction(lossOfFunction);
			break;
		case GATK:
			outputFormatter = new VcfOutputFormatter(config.getGenome(), VcfEffect.FormatVersion.FORMAT_SNPEFF_2);
			break;
		case BED:
			outputFormatter = new BedOutputFormatter();
			break;
		case BEDANN:
			outputFormatter = new BedAnnotationOutputFormatter();
			break;
		default:
			throw new RuntimeException("Unknown output format '" + outputFormat + "'");
		}
		outputFormatter.setVersion(VERSION);
		outputFormatter.setCommandLineStr(commandLineStr(false));
		outputFormatter.setChangeEffectResutFilter(changeEffectResutFilter);
		outputFormatter.setSupressOutput(supressOutput);
		outputFormatter.setOutOffset(outOffset);
		outputFormatter.setChrStr(chrStr);

		//---
		// Iterate over all changes
		//---
		switch (inputFormat) {
		case VCF:
			if (multiThreaded) iterateVcfMulti(outputFormatter);
			else iterateVcf(outputFormatter);
			break;
		default:
			iterateSeqChange(outputFormatter);
		}

		//---
		// Create reports
		//---
		if (createSummary && (summaryFile != null)) {
			// Creates a summary output file
			if (verbose) Timer.showStdErr("Creating summary file: " + summaryFile);
			summary(SUMMARY_TEMPLATE, summaryFile, false);

			// Creates genes output file
			if (verbose) Timer.showStdErr("Creating genes file: " + summaryGenesFile);
			summary(SUMMARY_GENES_TEMPLATE, summaryGenesFile, true);
		}

		if (totalErrs > 0) System.err.println(totalErrs + " errors.");
	}

	/**
	 * Creates a summary output file (using freeMarker and a template)
	 */
	void summary(String templateFile, String outputFile, boolean noCommas) {
		try {
			// Configure FreeMaker
			Configuration cfg = new Configuration();

			// Specify the data source where the template files come from 
			if (useLocalTemplate) cfg.setDirectoryForTemplateLoading(new File("./templates/")); // Use local 'template' directory 
			else cfg.setClassForTemplateLoading(SnpEffCmdEff.class, "/"); // Use current directory in JAR file

			cfg.setObjectWrapper(new DefaultObjectWrapper()); // Specify how templates will see the data-model. This is an advanced topic...
			cfg.setLocale(java.util.Locale.US);
			if (noCommas) cfg.setNumberFormat("0.######");

			// Create the root hash (where data objects are)
			HashMap<String, Object> root = summaryCreateHash();

			// Get the template
			Template temp = cfg.getTemplate(templateFile);

			// Process the template
			Writer out = new OutputStreamWriter(new FileOutputStream(new File(outputFile)));
			temp.process(root, out);
			out.flush();
			out.close();
		} catch (IOException e) {
			error(e, "Error creating summary: " + e.getMessage());
		} catch (TemplateException e) {
			error(e, "Error creating summary: " + e.getMessage());
		}
	}

	/**
	 * Create a hash with all variables needed for creating summary pages
	 * @return
	 */
	HashMap<String, Object> summaryCreateHash() {
		// Create the root hash (where data objects are)
		HashMap<String, Object> root = new HashMap<String, Object>();
		root.put("args", commandLineStr(true));
		root.put("changeStats", changeEffectResutStats);
		root.put("chromoPlots", chromoPlots);
		root.put("countEffects", countEffects);
		root.put("countInputLines", countInputLines);
		root.put("countVariants", countVariants);
		root.put("countVariantsFilteredOut", countVariantsFilteredOut);
		root.put("date", String.format("%1$TY-%1$Tm-%1$Td %1$TH:%1$TM", new Date()));
		root.put("genesFile", Gpr.baseName(summaryGenesFile, ""));
		root.put("genome", config.getGenome());
		root.put("genomeVersion", genomeVer);
		root.put("seqChangeFilter", seqChangeFilter);
		root.put("seqStats", seqChangeStats);
		root.put("snpEffectPredictor", config.getSnpEffectPredictor());
		root.put("vcfStats", vcfStats);
		root.put("version", SnpEff.VERSION); // Version used

		return root;
	}

	/**
	 * Show 'usage;' message and exit with an error code '-1'
	 * @param message
	 */
	@Override
	public void usage(String message) {
		if (message != null) {
			System.err.println("Error        :\t" + message);
			System.err.println("Command line :\t" + commandLineStr(false) + "\n");
		}

		System.err.println("snpEff version " + VERSION);
		System.err.println("Usage: snpEff [eff] [options] genome_version [variants_file]");
		System.err.println("\nInput file: Default is STDIN");
		System.err.println("\nOptions:");
		System.err.println("\t-a , -around            : Show N codons and amino acids around change (only in coding regions). Default is " + CodonChange.SHOW_CODONS_AROUND_CHANGE + " codons.");
		System.err.println("\t-i <format>             : Input format [ vcf, txt, pileup, bed ]. Default: VCF.");
		System.err.println("\t-o <format>             : Ouput format [ txt, vcf, gatk, bed, bedAnn ]. Default: VCF.");
		System.err.println("\t-interval               : Use a custom interval file (you may use this option many times)");
		System.err.println("\t-chr <string>           : Prepend 'string' to chromosome name (e.g. 'chr1' instead of '1'). Only on TXT output.");
		System.err.println("\t-s,  -stats             : Name of stats file (summary). Default is '" + DEFAULT_SUMMARY_FILE + "'");
		System.err.println("\t-t                      : Use multiple threads (implies '-noStats'). Default 'off'");
		System.err.println("\nSequence change filter options:");
		System.err.println("\t-del                    : Analyze deletions only");
		System.err.println("\t-ins                    : Analyze insertions only");
		System.err.println("\t-hom                    : Analyze homozygous variants only");
		System.err.println("\t-het                    : Analyze heterozygous variants only");
		System.err.println("\t-minQ X, -minQuality X  : Filter out variants with quality lower than X");
		System.err.println("\t-maxQ X, -maxQuality X  : Filter out variants with quality higher than X");
		System.err.println("\t-minC X, -minCoverage X : Filter out variants with coverage lower than X");
		System.err.println("\t-maxC X, -maxCoverage X : Filter out variants with coverage higher than X");
		System.err.println("\t-nmp                    : Only MNPs (multiple nucleotide polymorphisms)");
		System.err.println("\t-snp                    : Only SNPs (single nucleotide polymorphisms)");
		System.err.println("\nResults filter options:");
		System.err.println("\t-fi  <bedFile>                  : Only analyze changes that intersect with the intervals specified in this file (you may use this option many times)");
		System.err.println("\t-no-downstream                  : Do not show DOWNSTREAM changes");
		System.err.println("\t-no-intergenic                  : Do not show INTERGENIC changes");
		System.err.println("\t-no-intron                      : Do not show INTRON changes");
		System.err.println("\t-no-upstream                    : Do not show UPSTREAM changes");
		System.err.println("\t-no-utr                         : Do not show 5_PRIME_UTR or 3_PRIME_UTR changes");
		System.err.println("\nAnnotations filter options:");
		System.err.println("\t-canon                          : Only use canonical transcripts.");
		System.err.println("\t-lof                            : Add loss of function (LOF) and Nonsense mediated decay (NMD) tags.");
		System.err.println("\t-reg <name>                     : Regulation track to use (this option can be used add several times).");
		System.err.println("\t-onlyReg                        : Only use regulation tracks.");
		System.err.println("\t-onlyTr <file.txt>              : Only use the transcripts in this file. Format: One transcript ID per line.");
		System.err.println("\t-treatAllAsProteinCoding <bool> : If true, all transcript are treated as if they were protein conding. Default: Auto");
		System.err.println("\t-ud, -upDownStreamLen           : Set upstream downstream interval length (in bases)");
		System.err.println("\nGeneric options:");
		System.err.println("\t-0                      : File positions are zero-based (same as '-inOffset 0 -outOffset 0')");
		System.err.println("\t-1                      : File positions are one-based (same as '-inOffset 1 -outOffset 1')");
		System.err.println("\t-c , -config            : Specify config file");
		System.err.println("\t-h , -help              : Show this help and exit");
		System.err.println("\t-if, -inOffset          : Offset input by a number of bases. E.g. '-inOffset 1' for one-based input files");
		System.err.println("\t-of, -outOffset         : Offset output by a number of bases. E.g. '-outOffset 1' for one-based output files");
		System.err.println("\t-noLog                  : Do not report usage statistics to server");
		System.err.println("\t-noStats                : Do not create stats (summary) file");
		System.err.println("\t-q , -quiet             : Quiet mode (do not show any messages or errors)");
		System.err.println("\t-v , -verbose           : Verbose mode");
		System.exit(-1);
	}
}
