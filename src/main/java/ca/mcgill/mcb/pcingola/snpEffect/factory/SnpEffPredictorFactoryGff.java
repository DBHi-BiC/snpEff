package ca.mcgill.mcb.pcingola.snpEffect.factory;

import java.io.BufferedReader;
import java.util.HashMap;

import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.IntergenicConserved;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.interval.Utr3prime;
import ca.mcgill.mcb.pcingola.interval.Utr5prime;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * This class creates a SnpEffectPredictor from a GFF file.
 * This includes derived formats as GTF.
 * 
 * References: http://gmod.org/wiki/GFF3
 * 
 * @author pcingola
 */
public abstract class SnpEffPredictorFactoryGff extends SnpEffPredictorFactory {

	public static boolean debug = false;

	public static final HashMap<String, String> typeMap;
	public static final String GENE = Gene.class.getSimpleName();
	public static final String TRANSCRIPT = Transcript.class.getSimpleName();
	public static final String EXON = Exon.class.getSimpleName();
	public static final String UTR5 = Utr5prime.class.getSimpleName();
	public static final String UTR3 = Utr3prime.class.getSimpleName();
	public static final String INTERGENIC_CONSERVED = IntergenicConserved.class.getSimpleName();

	public static final String FASTA_DELIMITER = "##FASTA";

	String version = "";
	boolean mainFileHasFasta = false; // Are sequences in the GFF file or in a separate FASTA file?

	static {
		// Initialize typeMap
		typeMap = new HashMap<String, String>();

		addTypeMap("gene", GENE);
		addTypeMap("pseudogene", TRANSCRIPT);

		addTypeMap("mRNA", TRANSCRIPT);
		addTypeMap("tRNA", TRANSCRIPT);
		addTypeMap("snoRNA", TRANSCRIPT);
		addTypeMap("rRNA", TRANSCRIPT);
		addTypeMap("ncRNA", TRANSCRIPT);
		addTypeMap("miRNA", TRANSCRIPT);
		addTypeMap("snRNA", TRANSCRIPT);
		addTypeMap("pseudogenic_transcript", TRANSCRIPT);

		addTypeMap("exon", EXON);
		addTypeMap("pseudogenic_exon", EXON);
		addTypeMap("CDS", EXON);
		addTypeMap("start_codon", EXON);
		addTypeMap("stop_codon", EXON);
		addTypeMap("intron_CNS", EXON);

		addTypeMap("five_prime_UTR", UTR5);
		addTypeMap("5'-UTR", UTR5);
		addTypeMap("5UTR", UTR5);

		addTypeMap("three_prime_UTR", UTR3);
		addTypeMap("3'-UTR", UTR3);
		addTypeMap("3UTR", UTR3);

		addTypeMap("inter_CNS", INTERGENIC_CONSERVED);
	}

	static void addTypeMap(String typeAliasStr, String type) {
		typeMap.put(typeAliasStr.toUpperCase(), type);
	}

	public SnpEffPredictorFactoryGff(Config config, int inOffset) {
		super(config, inOffset);
		markersById = new HashMap<String, Marker>();
		genesById = new HashMap<String, Gene>();
		transcriptsById = new HashMap<String, Transcript>();
		fileName = config.getBaseFileNameGenes() + ".gff";
	}

	@Override
	public SnpEffectPredictor create() {
		// Read gene intervals from a file
		System.out.println("Reading " + version + " data file  : '" + fileName + "'");
		try {
			// We have to read the file a few times because we want to have all genes, then all transcripts, then all exons, etc.
			System.out.print("\tReading genes       : ");
			readGff(GENE);

			System.out.print("\tReading transcripts : ");
			readGff(TRANSCRIPT);

			System.out.print("\tReading exons       : ");
			readGff(EXON);

			// This features are not present in GFF2 and are optional in GTF 2.2
			if (!version.equals("GFF2")) {
				exonsFromCds(); // We need to create exons from CDSs before UTRs are added, since UTR require exons as parents

				System.out.print("\tReading UTRs (5)    : ");
				readGff(UTR5);

				System.out.print("\tReading UTRs (3)    : ");
				readGff(UTR3);
			}

			// Some clean-up before readng exon sequences
			beforeExonSequences();

			if (readSequences) {
				// Read chromosome sequences and set exon sequences
				System.out.print("\tReading sequences   :\n");
				if (mainFileHasFasta) readExonSequencesGff(fileName); // Read from GFF file (it has a '##FASTA' delimiter)
				else readExonSequences(); // Read them from FASTA file
			}

			System.out.println("\tTotal: " + totalSeqsAdded + " sequences added, " + totalSeqsIgnored + " sequences ignored.");

			// Finish up (fix problems, add missing info, etc.)
			finishUp(false);

			// Check that exons have sequences
			System.out.println(config.getGenome());
			boolean error = config.getGenome().isMostExonsHaveSequence();
			if (error && readSequences) throw new RuntimeException("Most Exons do not have sequences!");

		} catch (Exception e) {
			e.printStackTrace();
			throw new RuntimeException("Error reading file '" + fileName + "'\n" + e);
		}

		return snpEffectPredictor;
	}

	/**
	 * Is 'term' a 'type'?
	 */
	boolean is(String term, String type) {
		String inType = typeMap.get(term.toUpperCase());
		if (inType == null) return false;
		return inType.equals(type);
	}

	/**
	 * Is this protein coding according to the source
	 * @param source
	 * @return
	 */
	protected boolean isProteingCoding(String source) {
		return source.equals("protein_coding");

	}

	/**
	 * Parse a line
	 * @param line
	 * @return true if a line was parsed
	 */
	protected abstract boolean parse(String line, String typeToRead);

	/**
	 * Read chromosome sequence from GFF3 file and extract exons' sequences 
	 */
	protected void readExonSequencesGff(String gffFileName) {
		try {
			BufferedReader reader = Gpr.reader(gffFileName);

			// Get to fasta part of the file
			for (lineNum = 1; reader.ready(); lineNum++) {
				line = reader.readLine();
				if (line.equals(FASTA_DELIMITER)) {
					mainFileHasFasta = true;
					break;
				}
			}

			// Read fasta sequence
			String chromoName = null;
			StringBuffer chromoSb = new StringBuffer();
			for (; reader.ready(); lineNum++) {
				line = reader.readLine();
				if (line.startsWith(">")) { // New fasta sequence
					// Set chromosome sequences and length (create it if it doesn't exist)
					if (chromoName != null) {
						chromoLen(chromoName, chromoSb.length());
						addExonSequences(chromoName, chromoSb.toString()); // Add all sequences
					}

					chromoName = Chromosome.simpleName(line.substring(1).trim()); // New chromosome name
					chromoSb = new StringBuffer();
					System.out.println("\t\tReading sequence '" + chromoName + "'");
				} else chromoSb.append(line.trim());
			}

			// Last chromosome
			// Set chromosome sequneces and length (create it if it doesn't exist)
			if (chromoName != null) {
				chromoLen(chromoName, chromoSb.length());
				addExonSequences(chromoName, chromoSb.toString()); // Add all sequences
			} else System.err.println("WARNING: Ignoring sequences for '" + chromoName + "'. Cannot find chromosome"); // Chromosome not found

			reader.close();
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Read GFF file from the beginning looking for 'typeToRead' elements
	 * @param typeToRead
	 */
	protected void readGff(String typeToRead) throws Exception {
		int count = 0;
		BufferedReader reader = Gpr.reader(fileName);
		if (reader == null) return; // Error

		// Parsing GFF3 (reference: http://gmod.org/wiki/GFF3)
		try {
			for (lineNum = 1; reader.ready(); lineNum++) {
				line = reader.readLine();

				// Are we done?
				if (line.equals(FASTA_DELIMITER)) {
					mainFileHasFasta = true;
					break;
				} else if (line.startsWith("#")) {
					// Ignore this line
				} else if (parse(line, typeToRead)) {
					count++;
					Gpr.showMark(count, MARK);
				}
			}
		} catch (Exception e) {
			error("Offending line (lineNum: " + lineNum + "): '" + line + "'");
		}

		reader.close();
		System.out.println((count > 0 ? "\n" : "") + "\tTotal: " + count + " " + typeToRead + "s added.");
	}
}
