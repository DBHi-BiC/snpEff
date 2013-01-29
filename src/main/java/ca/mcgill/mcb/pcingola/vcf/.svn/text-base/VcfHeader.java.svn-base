package ca.mcgill.mcb.pcingola.vcf;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;

/**
 * Represents the header of a vcf file.
 * 
 * References: http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
 * 
 * @author pablocingolani
 *
 */
public class VcfHeader {

	StringBuffer header;
	HashMap<String, VcfInfo> vcfInfoById;

	public VcfHeader() {
		header = new StringBuffer();
	}

	/**
	 * Add a 'FORMAT' meta info
	 * @param vcfGenotypeStr
	 */
	public void addFormat(String formatName, String number, String type, String description) {
		String headerFormatInfo = "##FORMAT=<ID=" + formatName + ",Number=" + number + ",Type=" + type + ",Description=\"" + description + "\">";
		addLine(headerFormatInfo);
	}

	/**
	 * Add line to header (can add many lines)
	 * @return
	 */
	public void addLine(String newHeaderLine) {
		// Split header
		String headerLines[] = header.toString().split("\n");
		header = new StringBuffer();

		// Find "#CHROM" line in header (should always be the last one)
		boolean added = false;
		for (String line : headerLines) {
			if (line.equals(newHeaderLine)) {
				newHeaderLine = null; // Line already present? => Don't add
				added = true;
			} else if (line.startsWith("#CHROM") && (newHeaderLine != null)) {
				header.append(newHeaderLine + "\n"); // Add new header right before title line
				added = true;
			}

			if (!line.isEmpty()) header.append(line + "\n"); // Add non-empty lines
		}

		// Not added yet? => Add to the end
		if (!added) header.append(newHeaderLine + "\n"); // Add new header right before title line
	}

	/**
	 * Get sample names
	 * @return
	 */
	public List<String> getSampleNames() {
		// Split header
		String headerLines[] = header.toString().split("\n");

		// Find "#CHROM" line in header
		for (String line : headerLines) {
			if (line.startsWith("#CHROM")) {
				// This line contains all the sample names (starting on column 9)
				String titles[] = line.split("\t");

				// Create a list of names
				ArrayList<String> sampleNames = new ArrayList<String>();
				for (int i = 9; i < titles.length; i++)
					sampleNames.add(titles[i]);

				// Done
				return sampleNames;
			}
		}

		// Not found
		return null;
	}

	/**
	 * Get all VcfInfo entries
	 * @return
	 */
	public Collection<VcfInfo> getVcfInfo() {
		parseInfoLines();
		return vcfInfoById.values();
	}

	/**
	 * Get Info type for a given ID
	 * @param id
	 * @return
	 */
	public VcfInfo getVcfInfo(String id) {
		parseInfoLines();
		return vcfInfoById.get(id);
	}

	/**
	 * Parse INFO fields from header
	 */
	public void parseInfoLines() {
		if (vcfInfoById == null) {
			vcfInfoById = new HashMap<String, VcfInfo>();

			// Add standard fields
			vcfInfoById.put("CHROM", new VcfInfo("CHROM", VcfInfoType.String, "1", "Chromosome name"));
			vcfInfoById.put("POS", new VcfInfo("POS", VcfInfoType.Integer, "1", "Position in chromosome"));
			vcfInfoById.put("ID", new VcfInfo("ID", VcfInfoType.String, "1", "Variant ID"));
			vcfInfoById.put("REF", new VcfInfo("REF", VcfInfoType.String, "1", "Reference sequence"));
			vcfInfoById.put("ALT", new VcfInfo("ALT", VcfInfoType.String, "A", "Alternative sequence/s"));
			vcfInfoById.put("QUAL", new VcfInfo("QUAL", VcfInfoType.Float, "1", "Mapping quality"));
			vcfInfoById.put("FILTER", new VcfInfo("FILTER", VcfInfoType.String, "1", "Filter status"));
			vcfInfoById.put("FORMAT", new VcfInfo("FORMAT", VcfInfoType.String, "1", "Format in genotype fields"));

			// Add well known fields 
			// Reference: http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
			vcfInfoById.put("AA", new VcfInfo("AA", VcfInfoType.String, "1", "Ancestral allele"));
			vcfInfoById.put("AC", new VcfInfo("AC", VcfInfoType.Integer, "A", "Allele Frequency"));
			vcfInfoById.put("AF", new VcfInfo("AF", VcfInfoType.Float, "1", "Allele Frequency"));
			vcfInfoById.put("AN", new VcfInfo("AN", VcfInfoType.Integer, "1", "Total number of alleles"));
			vcfInfoById.put("BQ", new VcfInfo("BQ", VcfInfoType.Float, "1", "RMS base quality"));
			vcfInfoById.put("CIGAR", new VcfInfo("CIGAR", VcfInfoType.String, "1", "Cigar string describing how to align an alternate allele to the reference allele"));
			vcfInfoById.put("DB", new VcfInfo("DB", VcfInfoType.Flag, "1", "dbSNP membership"));
			vcfInfoById.put("DP", new VcfInfo("DP", VcfInfoType.Integer, "1", "Combined depth across samples"));
			vcfInfoById.put("END", new VcfInfo("END", VcfInfoType.String, "1", "End position of the variant described in this record"));
			vcfInfoById.put("H2", new VcfInfo("H2", VcfInfoType.Flag, "1", "Membership in hapmap 2"));
			vcfInfoById.put("H3", new VcfInfo("H3", VcfInfoType.Flag, "1", "Membership in hapmap 3"));
			vcfInfoById.put("MQ", new VcfInfo("MQ", VcfInfoType.Float, "1", "RMS mapping quality"));
			vcfInfoById.put("MQ0", new VcfInfo("MQ0", VcfInfoType.Integer, "1", "Number of MAPQ == 0 reads covering this record"));
			vcfInfoById.put("NS", new VcfInfo("NS", VcfInfoType.Integer, "1", "Number of samples with data"));
			vcfInfoById.put("SB", new VcfInfo("SB", VcfInfoType.Float, "1", "Strand bias at this position"));
			vcfInfoById.put("SOMATIC", new VcfInfo("SOMATIC", VcfInfoType.Flag, "1", "Indicates that the record is a somatic mutation, for cancer genomics"));
			vcfInfoById.put("VALIDATED", new VcfInfo("VALIDATED", VcfInfoType.Flag, "1", "Validated by follow-up experiment"));
			vcfInfoById.put("1000G", new VcfInfo("1000G", VcfInfoType.Flag, "1", "Membership in 1000 Genomes"));

			// Structural variants
			vcfInfoById.put("IMPRECISE", new VcfInfo("IMPRECISE", VcfInfoType.Flag, "0", "Imprecise structural variation"));
			vcfInfoById.put("NOVEL", new VcfInfo("NOVEL", VcfInfoType.Flag, "0", "Indicates a novel structural variation"));
			vcfInfoById.put("END", new VcfInfo("END", VcfInfoType.Integer, "1", "End position of the variant described in this record"));
			vcfInfoById.put("SVTYPE", new VcfInfo("SVTYPE", VcfInfoType.String, "1", "Type of structural variant"));
			vcfInfoById.put("SVLEN", new VcfInfo("SVLEN", VcfInfoType.Integer, ".", "Difference in length between REF and ALT alleles"));
			vcfInfoById.put("CIPOS", new VcfInfo("CIPOS", VcfInfoType.Integer, "2", "Confidence interval around POS for imprecise variants"));
			vcfInfoById.put("CIEND", new VcfInfo("CIEND", VcfInfoType.Integer, "2", "Confidence interval around END for imprecise variants"));
			vcfInfoById.put("HOMLEN", new VcfInfo("HOMLEN", VcfInfoType.Integer, ".", "Length of base pair identical micro-homology at event breakpoints"));
			vcfInfoById.put("HOMSEQ", new VcfInfo("HOMSEQ", VcfInfoType.String, ".", "Sequence of base pair identical micro-homology at event breakpoints"));
			vcfInfoById.put("BKPTID", new VcfInfo("BKPTID", VcfInfoType.String, ".", "ID of the assembled alternate allele in the assembly file"));
			vcfInfoById.put("MEINFO", new VcfInfo("MEINFO", VcfInfoType.String, "4", "Mobile element info of the form NAME,START,END,POLARITY"));
			vcfInfoById.put("METRANS", new VcfInfo("METRANS", VcfInfoType.String, "4", "Mobile element transduction info of the form CHR,START,END,POLARITY"));
			vcfInfoById.put("DGVID", new VcfInfo("DGVID", VcfInfoType.String, "1", "ID of this element in Database of Genomic Variation"));
			vcfInfoById.put("DBVARID", new VcfInfo("DBVARID", VcfInfoType.String, "1", "ID of this element in DBVAR"));
			vcfInfoById.put("DBRIPID", new VcfInfo("DBRIPID", VcfInfoType.String, "1", "ID of this element in DBRIP"));
			vcfInfoById.put("MATEID", new VcfInfo("MATEID", VcfInfoType.String, ".", "ID of mate breakends"));
			vcfInfoById.put("PARID", new VcfInfo("PARID", VcfInfoType.String, "1", "ID of partner breakend"));
			vcfInfoById.put("EVENT", new VcfInfo("EVENT", VcfInfoType.String, "1", "ID of event associated to breakend"));
			vcfInfoById.put("CILEN", new VcfInfo("CILEN", VcfInfoType.Integer, "2", "Confidence interval around the length of the inserted material between breakends"));
			vcfInfoById.put("DP", new VcfInfo("DP", VcfInfoType.Integer, "1", "Read Depth of segment containing breakend"));
			vcfInfoById.put("DPADJ", new VcfInfo("DPADJ", VcfInfoType.Integer, ".", "Read Depth of adjacency"));
			vcfInfoById.put("CN", new VcfInfo("CN", VcfInfoType.Integer, "1", "Copy number of segment containing breakend"));
			vcfInfoById.put("CNADJ", new VcfInfo("CNADJ", VcfInfoType.Integer, ".", "Copy number of adjacency"));
			vcfInfoById.put("CICN", new VcfInfo("CICN", VcfInfoType.Integer, "2", "Confidence interval around copy number for the segment"));
			vcfInfoById.put("CICNADJ", new VcfInfo("CICNADJ", VcfInfoType.Integer, ".", "Confidence interval around copy number for the adjacency"));

			// Add SnpEff 'EFF' fields
			vcfInfoById.put("EFF.EFFECT", new VcfInfo("EFF.EFFECT", VcfInfoType.String, ".", "SnpEff effect"));
			vcfInfoById.put("EFF.IMPACT", new VcfInfo("EFF.IMPACT", VcfInfoType.String, ".", "SnpEff impact (HIGH, MODERATE, LOW, MODIFIER)"));
			vcfInfoById.put("EFF.FUNCLASS", new VcfInfo("EFF.FUNCLASS", VcfInfoType.String, ".", "SnpEff functional class (NONE, SILENT, MISSENSE, NONSENSE)"));
			vcfInfoById.put("EFF.CODON", new VcfInfo("EFF.CODON", VcfInfoType.String, ".", "SnpEff codon change"));
			vcfInfoById.put("EFF.AA", new VcfInfo("EFF.AA", VcfInfoType.String, ".", "SnpEff amino acid change"));
			vcfInfoById.put("EFF.AA_LEN", new VcfInfo("EFF.AA_LEN", VcfInfoType.Integer, ".", "Protein length in amino acids"));
			vcfInfoById.put("EFF.GENE", new VcfInfo("EFF.GENE", VcfInfoType.String, ".", "SnpEff gene name"));
			vcfInfoById.put("EFF.BIOTYPE", new VcfInfo("EFF.BIOTYPE", VcfInfoType.String, ".", "SnpEff gene bio-type"));
			vcfInfoById.put("EFF.CODING", new VcfInfo("EFF.CODING", VcfInfoType.String, ".", "SnpEff gene coding (CODING, NON_CODING)"));
			vcfInfoById.put("EFF.TRID", new VcfInfo("EFF.TRID", VcfInfoType.String, ".", "SnpEff transcript ID"));
			vcfInfoById.put("EFF.EXID", new VcfInfo("EFF.EXID", VcfInfoType.String, ".", "SnpEff exon ID"));

			// Add SnpEff 'LOF' fields
			vcfInfoById.put("LOF.GENE", new VcfInfo("LOF.GENE", VcfInfoType.String, ".", "SnpEff LOF gene name"));
			vcfInfoById.put("LOF.GENEID", new VcfInfo("LOF.GENEID", VcfInfoType.String, ".", "SnpEff LOF gene ID"));
			vcfInfoById.put("LOF.NUMTR", new VcfInfo("LOF.NUMTR", VcfInfoType.Integer, ".", "SnpEff LOF number of transcripts in gene"));
			vcfInfoById.put("LOF.PERC", new VcfInfo("LOF.PERC", VcfInfoType.Float, ".", "SnpEff LOF percentage of transcripts in this gene that are affected"));

			// Add SnpEff 'NMD' fields
			vcfInfoById.put("NMD.GENE", new VcfInfo("NMD.GENE", VcfInfoType.String, ".", "SnpEff NMD gene name"));
			vcfInfoById.put("NMD.GENEID", new VcfInfo("NMD.GENEID", VcfInfoType.String, ".", "SnpEff NMD gene ID"));
			vcfInfoById.put("NMD.NUMTR", new VcfInfo("NMD.NUMTR", VcfInfoType.Integer, ".", "SnpEff NMD number of transcripts in gene"));
			vcfInfoById.put("NMD.PERC", new VcfInfo("NMD.PERC", VcfInfoType.Float, ".", "SnpEff NMD percentage of transcripts in this gene that are affected"));

			// Add all INFO fields from header
			String headerLines[] = header.toString().split("\n");
			for (String line : headerLines) {
				if (line.startsWith("##INFO=") || line.startsWith("##FORMAT=")) {
					VcfInfo vcfInfo = new VcfInfo(line);
					vcfInfoById.put(vcfInfo.getId(), vcfInfo);
				}
			}
		}
	}

	/**
	 * Get header information
	 * @return
	 */
	@Override
	public String toString() {
		if (header.length() <= 0) return "";

		// Delete last character, if it's a '\n' or a '\r'
		for (char c = header.charAt(header.length() - 1); (c == '\n') || (c == '\r'); c = header.charAt(header.length() - 1))
			header.deleteCharAt(header.length() - 1);

		return header.toString();
	}

}
