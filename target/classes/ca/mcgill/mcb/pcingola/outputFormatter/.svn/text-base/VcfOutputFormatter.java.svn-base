package ca.mcgill.mcb.pcingola.outputFormatter;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.interval.Intron;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Regulation;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.FunctionalClass;
import ca.mcgill.mcb.pcingola.snpEffect.LossOfFunction;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.vcf.VcfEffect;
import ca.mcgill.mcb.pcingola.vcf.VcfEffect.FormatVersion;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;

/**
 * Formats output as VCF
 * 
 * @author pcingola
 */
public class VcfOutputFormatter extends OutputFormatter {

	public static final boolean debug = false;

	public static final String VCF_INFO_EFF_NAME = "EFF";
	public static final String VCF_INFO_LOF_NAME = "LOF";
	public static final String VCF_INFO_NMD_NAME = "NMD";

	boolean needAddInfo = false;
	boolean needAddHeader = true;
	FormatVersion formatVersion = VcfEffect.FormatVersion.FORMAT_SNPEFF_3;
	List<VcfEntry> vcfEntries;
	boolean lossOfFunction;
	Genome genome;

	/**
	 * This constructor is used mostly by 'clone()' method
	 */
	protected VcfOutputFormatter() {
		super();
	}

	public VcfOutputFormatter(Genome genome, FormatVersion formatVersion) {
		super();
		this.genome = genome;
		this.formatVersion = formatVersion;
	}

	/**
	 * Add all vcf entries to a list (used only for debugging and test-cases)
	 * @param vcfEntries
	 */
	public VcfOutputFormatter(List<VcfEntry> vcfEntries) {
		super();
		this.vcfEntries = vcfEntries;
	}

	/**
	 * Add header
	 * @param vcfEntry
	 */
	protected void addHeader() {
		VcfEntry vcfEntry = (VcfEntry) section;

		// Get header
		VcfFileIterator vcfFile = vcfEntry.getVcfFileIterator();

		// Add new lines
		for (String newHeaderLine : getNewHeaderLines())
			vcfFile.getVcfHeader().addLine(newHeaderLine);

		needAddHeader = false;
	}

	/**
	 * Add effects to INFO field
	 */
	protected void addInfo(VcfEntry vcfEntry) {
		// No effects to show?
		if (changeEffects.isEmpty()) return;

		//---
		// Calculate all effects and genes
		//---
		HashSet<String> effs = new HashSet<String>();
		for (ChangeEffect changeEffect : changeEffects) {
			// If it is not filtered out by changeEffectResutFilter  => Show it
			if ((changeEffectResutFilter == null) || (!changeEffectResutFilter.filter(changeEffect))) {
				StringBuilder effBuff = new StringBuilder();

				// Add effect
				effBuff.append(changeEffect.effect(true, false, false));
				effBuff.append("(");

				// Add effect impact
				effBuff.append(changeEffect.getEffectImpact());
				effBuff.append("|");

				// Add functional class
				FunctionalClass fc = changeEffect.getFunctionalClass();
				effBuff.append(fc == FunctionalClass.NONE ? "" : fc.toString()); // Show only if it is not empty
				effBuff.append("|");

				// Codon change
				effBuff.append(changeEffect.getCodonChange());
				effBuff.append("|");

				// Add amino acid change
				effBuff.append(changeEffect.getAaChangeHgvs());
				effBuff.append("|");

				// Add amino acid length
				if (formatVersion != FormatVersion.FORMAT_SNPEFF_2) { // This field is not in format version 2
					int aalen = changeEffect.getAaLength();
					effBuff.append(aalen >= 0 ? aalen : "");
					effBuff.append("|");
				}

				// Add gene info
				Gene gene = changeEffect.getGene();
				Transcript tr = changeEffect.getTranscript();
				if (gene != null) {
					// Gene name
					effBuff.append(vcfInfoSafeString(gene.getGeneName()));
					effBuff.append("|");

					// Transcript ID
					effBuff.append(tr != null ? tr.getBioType() : "");
					effBuff.append("|");

					// Protein coding gene?
					String coding = "";
					if (gene.getGenome().hasCodingInfo()) coding = (gene.isProteinCoding() ? ChangeEffect.Coding.CODING.toString() : ChangeEffect.Coding.NON_CODING.toString());
					effBuff.append(coding);
					effBuff.append("|");
				} else if (changeEffect.isRegulation()) {
					Regulation reg = (Regulation) changeEffect.getMarker();
					effBuff.append("|" + reg.getCellType() + "||");
				} else if (changeEffect.isCustom()) {
					Marker m = changeEffect.getMarker();
					if (m != null) effBuff.append("|" + m.getId() + "||");
					else effBuff.append("|||");
				} else effBuff.append("|||");

				// Add transcript info
				if (tr != null) effBuff.append(vcfInfoSafeString(tr.getId()));
				effBuff.append("|");

				// Add exon (or intron) rank info
				Exon ex = changeEffect.getExon();
				int rank = -1;
				if (ex != null) rank = ex.getRank();
				else {
					// Do we have an intron?
					Intron intron = changeEffect.getIntron();
					if (intron != null) rank = intron.getRank();
				}

				effBuff.append(rank >= 0 ? rank : "");

				// Errors or warnings (this is the last thing in the list)
				if (!changeEffect.getWarning().isEmpty()) effBuff.append("|" + changeEffect.getWarning());
				if (!changeEffect.getError().isEmpty()) effBuff.append("|" + changeEffect.getError());

				effBuff.append(")");

				//---
				// Add effect
				//---
				if (!effs.add(effBuff.toString())) {
					if (debug) {
						// Effect has already been added? Something is wrong, the information should be unique for each effect
						StringBuilder sb = new StringBuilder();
						sb.append("--------------------------------------------------------------------------------\n");
						sb.append("VCF Entry   :\t" + vcfEntry + "\n");
						sb.append("REPEAT (VCF):\t" + effBuff + "\n");
						sb.append("REPEAT (TXT):\t" + changeEffect + "\n");
						sb.append("All    (VCF):\n");
						for (String ce : effs)
							sb.append("\t" + ce + "\n");
						sb.append("All    (TXT):\n");
						for (ChangeEffect ce : changeEffects)
							sb.append("\t" + ce + "\n");
						sb.append("--------------------------------------------------------------------------------\n");
						Gpr.debug("WARNING: Repeated effect!\n" + sb);
					}
				}
			}
		}

		//---
		// Add effects (sorted)
		//---
		ArrayList<String> listEffs = new ArrayList<String>(effs);
		Collections.sort(listEffs);
		StringBuffer sbEffs = new StringBuffer();
		for (String eff : listEffs)
			sbEffs.append(eff + ",");

		if (sbEffs.length() > 0) sbEffs.deleteCharAt(sbEffs.length() - 1); // Remove last comma

		// Add 'EFF' info field
		vcfEntry.addInfo(VCF_INFO_EFF_NAME, sbEffs.toString());

		//---
		// Add LOF info?
		//---
		if (lossOfFunction) {
			// Perform LOF analysis and add annotations
			LossOfFunction lof = new LossOfFunction(changeEffects);
			if (lof.isLof()) vcfEntry.addInfo(VCF_INFO_LOF_NAME, lof.vcfLofValue());
			if (lof.isNmd()) vcfEntry.addInfo(VCF_INFO_NMD_NAME, lof.vcfNmdValue());
		}

		needAddInfo = false; // Don't add info twice
	}

	@Override
	public OutputFormatter clone() {
		try {
			VcfOutputFormatter newOutputFormatter = (VcfOutputFormatter) super.clone();
			newOutputFormatter.formatVersion = formatVersion;
			newOutputFormatter.needAddInfo = needAddInfo;
			newOutputFormatter.needAddHeader = needAddHeader;
			newOutputFormatter.lossOfFunction = lossOfFunction;
			newOutputFormatter.genome = genome;
			return newOutputFormatter;
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Finish up section
	 * @param marker
	 */
	@Override
	public String endSection(Marker marker) {
		// Ignore other markers (e.g. seqChanges)
		if (marker instanceof VcfEntry) {
			if (vcfEntries != null) vcfEntries.add((VcfEntry) marker);
			return super.endSection(marker);
		}
		return null;
	}

	/**
	 * New lines to be added to header
	 * @return
	 */
	public List<String> getNewHeaderLines() {
		ArrayList<String> newLines = new ArrayList<String>();
		newLines.add("##SnpEffVersion=\"" + version + "\"");
		newLines.add("##SnpEffCmd=\"" + commandLineStr + "\"");
		newLines.add("##INFO=<ID=EFF,Number=.,Type=String,Description=\"Predicted effects for this variant.Format: 'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change| Amino_Acid_length | Gene_Name | Gene_BioType | Coding | Transcript | Exon [ | ERRORS | WARNINGS ] )' \">");
		if (lossOfFunction) {
			newLines.add("##INFO=<ID=LOF,Number=.,Type=String,Description=\"Predicted loss of function effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected' \">");
			newLines.add("##INFO=<ID=NMD,Number=.,Type=String,Description=\"Predicted nonsense mediated decay effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected' \">");
		}
		return newLines;
	}

	public void setLossOfFunction(boolean lossOfFunction) {
		this.lossOfFunction = lossOfFunction;
	}

	@Override
	public void startSection(Marker marker) {
		// Ignore other markers (e.g. seqChanges)
		if (marker instanceof VcfEntry) super.startSection(marker);
		needAddInfo = true;
	}

	@Override
	public String toString() {
		VcfEntry vcfEntry = (VcfEntry) section;
		if (needAddInfo) addInfo(vcfEntry);
		return vcfEntry.toString();
	}

	/**
	 * Show header
	 */
	@Override
	protected String toStringHeader() {
		if (needAddHeader) addHeader(); // Add header lines

		VcfEntry vcfEntry = (VcfEntry) section;
		VcfFileIterator vcfFile = vcfEntry.getVcfFileIterator();
		return vcfFile.getVcfHeader().toString();
	}

	/**
	 * Create a string that is safe (i.e. valid) to add in an INFO field
	 */
	public String vcfInfoSafeString(String value) {
		if (value == null) return value;
		value = value.replaceAll("[ ,;|=()]", "_");
		return value;
	}
}
