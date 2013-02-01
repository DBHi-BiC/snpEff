package ca.mcgill.mcb.pcingola.vcf;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import ca.mcgill.mcb.pcingola.fileIterator.NeedlemanWunsch;
import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.SeqChange.ChangeType;
import ca.mcgill.mcb.pcingola.outputFormatter.VcfOutputFormatter;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.vcf.VcfEffect.FormatVersion;

/**
 * A VCF entry (a line) in a VCF file

 * @author pablocingolani
 */
public class VcfEntry extends Marker implements Iterable<VcfGenotype> {

	private static final long serialVersionUID = 4226374412681243433L;

	String line; // Line from VCF file
	int lineNum; // Line number
	VcfFileIterator vcfFileIterator; // Iterator where this entry was red from
	String chromosomeName; // Original chromosome name
	String ref;
	String[] alts;
	Double quality;
	String filterPass;
	String infoStr = "";
	HashMap<String, String> info;
	String format;
	ArrayList<VcfGenotype> vcfGenotypes = null;
	ChangeType changeType;
	String genotypeFields[]; // Raw fields from VCF file
	String genotypeFieldsStr; // Raw fields from VCF file (one string, tab separated)

	public VcfEntry(VcfFileIterator vcfFileIterator, Marker parent, String chromosomeName, int start, String id, String ref, String altsStr, double quality, String filterPass, String infoStr, String format) {
		super(parent, start, start + ref.length(), 1, id);
		this.chromosomeName = chromosomeName;
		this.ref = ref;
		parseAlts(altsStr);
		this.quality = quality;
		this.filterPass = filterPass;
		this.infoStr = infoStr;
		parseInfo();
		this.format = format;
	}

	/**
	 * Create a line form a file iterator
	 * @param vcfFileIterator
	 * @param line
	 * @param lineNum
	 */
	public VcfEntry(VcfFileIterator vcfFileIterator, String line, int lineNum, boolean parseNow) {
		super(null, 0, 0, 1, "");
		this.vcfFileIterator = vcfFileIterator;
		this.lineNum = lineNum;
		this.line = line;

		if (parseNow) parse();
	}

	/**
	 * Add a 'FORMAT' field
	 * @param vcfGenotypeStr
	 */
	public void addFormat(String formatName) {
		if (format == null) format = "";
		if (format.indexOf(formatName) >= 0) throw new RuntimeException("Format field '" + formatName + "' already exists!");

		// Add to format
		format += (format.endsWith(":") ? "" : ":") + formatName;
	}

	/**
	 * Add a genotype as a string
	 * @param vcfGenotypeStr
	 */
	public void addGenotype(String vcfGenotypeStr) {
		if (vcfGenotypes == null) vcfGenotypes = new ArrayList<VcfGenotype>();
		if (format == null) format = "";
		vcfGenotypes.add(new VcfGenotype(this, format, vcfGenotypeStr));
	}

	/**
	 * Append a 'raw' INFO string (format is not checked)
	 * WARNING: Info fields are NOT added to the hash, so trying to retrieve them using 'getInfo(name)' will fail!
	 * 
	 * @param name
	 * @param value
	 */
	public void addInfo(String addInfoStr) {
		addInfo(addInfoStr, true);
	}

	/**
	 * Append a 'raw' INFO string (format is not checked)
	 * WARNING: Info fields are NOT added to the hash, so trying to retrieve them using 'getInfo(name)' will fail!
	 * 
	 * @param name
	 * @param value
	 */
	protected void addInfo(String addInfoStr, boolean invalidateCache) {
		if ((infoStr == null) || infoStr.isEmpty()) infoStr = addInfoStr;
		else {
			if (!infoStr.endsWith(";")) infoStr += ";"; // Do we need to add a semicolon?
			infoStr += addInfoStr; // Add info string
		}

		if (invalidateCache) addInfoStr = null; // Invalidate cache
	}

	/**
	 * Add a "key=value" tuple the info field
	 * 
	 * @param name
	 * @param value : Can be null if it is a boolean field.
	 */
	public void addInfo(String name, String value) {
		if ((value != null) && ((value.indexOf(' ') >= 0) || (value.indexOf(';') >= 0) || (value.indexOf('=') >= 0) || (value.indexOf('\t') >= 0) || (value.indexOf('\n') >= 0))) throw new RuntimeException("No white-space, semi-colons, or equals-signs are permitted in INFO field. Name:\"" + name + "\" Value:\"" + value + "\"");

		String addInfoStr = name + (value != null ? "=" + value : "");
		if (info != null) info.put(name, value); // Add to info hash (if available)
		addInfo(addInfoStr, false);
	}

	/**
	 * Is this entry heterozygous?
	 * 
	 * 		Infer Hom/Her if there is only one sample in the file.
	 * 		Ohtherwise the field is null.
	 * 
	 * @return
	 */
	public Boolean calcHetero() {
		// No genotyping information? => Use number of ALT fielsd
		if (genotypeFieldsStr == null) return isHeterozygous();

		Boolean isHetero = null;

		// If there is only one genotype field => parse fields
		if (genotypeFields == null) {

			// Are there more than two tabs? (i.e. more than one format field + one genotype field)  
			int countFields, fromIndex;
			for (countFields = 0, fromIndex = 0; (fromIndex >= 0) && (countFields < 1); countFields++, fromIndex++)
				fromIndex = genotypeFieldsStr.indexOf('\t', fromIndex);

			// OK only one genotype field => Parse it in order to extract homo info.
			if (countFields == 1) parseGenotypes();
		}

		// OK only one genotype field => calculate if it is heterozygous
		if ((genotypeFields != null) && (genotypeFields.length == 1)) isHetero = getVcfGenotype(0).isHeterozygous();

		return isHetero;
	}

	/**
	 * Create a seqChange 
	 * @return
	 */
	SeqChange createSeqChange(Chromosome chromo, int start, String reference, String alt, int strand, String id, double quality, int coverage) {
		// No change?
		if (alt.isEmpty()) return new SeqChange(chromo, start, reference, reference, strand, id, quality, coverage);

		alt = alt.toUpperCase();

		// Case: Structural variant
		// 2 321682    .  T   <DEL>         6     PASS    IMPRECISE;SVTYPE=DEL;END=321887;SVLEN=-105;CIPOS=-56,20;CIEND=-10,62
		if (alt.startsWith("<DEL")) {
			int end = start + reference.length() - 1;

			// If there is an 'END' tag, we should use it
			if ((getInfo("END") != null)) {
				// Get 'END' field and do some sanity check
				end = (int) getInfoInt("END");
				if (end < start) throw new RuntimeException("INFO field 'END' is before varaint's 'POS'\n\tEND : " + end + "\n\tPOS : " + start);
			}

			// Create deletion string
			// TODO: This should be changed. We should be using "imprecise" for these variants
			int size = end - start + 1;
			char change[] = new char[size];
			for (int i = 0; i < change.length; i++)
				change[i] = reference.length() > i ? reference.charAt(i) : 'N';
			String ch = "-" + new String(change);

			// Create SeqChange
			return new SeqChange(chromo, start, reference, ch, strand, id, quality, coverage);
		}

		// Case: SNP, MNP
		// 20     3 .         C      G       .   PASS  DP=100
		// 20     3 .         TC     AT      .   PASS  DP=100
		if (reference.length() == alt.length()) {
			int startDiff = Integer.MAX_VALUE;
			String ch = "";
			String ref = "";
			for (int i = 0; i < reference.length(); i++) {
				if (reference.charAt(i) != alt.charAt(i)) {
					ref += reference.charAt(i);
					ch += alt.charAt(i);
					startDiff = Math.min(startDiff, i);
				}
			}

			return new SeqChange(chromo, start + startDiff, ref, ch, strand, id, quality, coverage);
		}

		// Case: Deletion 
		// 20     2 .         TC      T      .   PASS  DP=100	
		// 20     2 .         AGAC    AAC    .   PASS  DP=100	
		if (reference.length() > alt.length()) {
			NeedlemanWunsch nw = new NeedlemanWunsch(alt, reference);
			nw.align();
			int startDiff = nw.getOffset();
			String ref = "*";
			String ch = nw.getAlignment();
			if (!ch.startsWith("-")) throw new RuntimeException("Deletion '" + ch + "' does not start with '-'. This should never happen!");

			return new SeqChange(chromo, start + startDiff, ref, ch, strand, id, quality, coverage);
		}

		// Case: Insertion of A { tC ; tCA } tC is the reference allele
		// 20     2 .         TC      TCA    .   PASS  DP=100
		if (reference.length() < alt.length()) {
			NeedlemanWunsch nw = new NeedlemanWunsch(alt, reference);
			nw.align();
			int startDiff = nw.getOffset();
			String ch = nw.getAlignment();
			String ref = "*";
			if (!ch.startsWith("+")) throw new RuntimeException("Insertion '" + ch + "' does not start with '+'. This should never happen!");

			return new SeqChange(chromo, start + startDiff, ref, ch, strand, id, quality, coverage);
		}

		// Other change type?
		throw new RuntimeException("Unsupported VCF change type '" + reference + "' => '" + alt + "'\nVcfEntry: " + this);
	}

	public String[] getAlts() {
		return alts;
	}

	/**
	 * Create a comma separated ALTS string
	 * @return
	 */
	public String getAltsStr() {
		String altsStr = "";
		for (String alt : alts)
			altsStr += alt + " ";
		return altsStr.trim().replace(' ', ',');
	}

	public ChangeType getChangeType() {
		return changeType;
	}

	/**
	 * Original chromosome name (as it appeared in the VCF file)
	 * @return
	 */
	@Override
	public String getChromosomeNameOri() {
		return chromosomeName;
	}

	public String getFilterPass() {
		return filterPass;
	}

	public String getFormat() {
		return format;
	}

	public String getInfo(String key) {
		if (info == null) parseInfo();
		return info.get(key);
	}

	/**
	 * Does the entry exists?
	 * @param key
	 * @return
	 */
	public boolean getInfoFlag(String key) {
		if (info == null) parseInfo();
		return info.containsKey(key);
	}

	/**
	 * Get info field as a 'double' number
	 * The norm specifies data type as 'FLOAT', that is why the name of this method might be not intuitive
	 * @param key
	 * @return
	 */
	public double getInfoFloat(String key) {
		if (info == null) parseInfo();
		return Gpr.parseDoubleSafe(info.get(key));
	}

	/**
	 * Get info field as an long number
	 * The norm specifies data type as 'INT', that is why the name of this method might be not intuitive
	 * @param key
	 * @return
	 */
	public long getInfoInt(String key) {
		if (info == null) parseInfo();
		return Gpr.parseLongSafe(info.get(key));
	}

	/**
	 * Get all keys available in the info field
	 * @param key
	 * @return
	 */
	public Set<String> getInfoKeys(String key) {
		if (info == null) parseInfo();
		return info.keySet();
	}

	/**
	 * Get the full (unparsed) INFO field
	 * @return
	 */
	public String getInfoStr() {
		return infoStr;
	}

	public int getLineNum() {
		return lineNum;
	}

	public double getQuality() {
		return (quality != null ? quality : 0);
	}

	public String getRef() {
		return ref;
	}

	public VcfFileIterator getVcfFileIterator() {
		return vcfFileIterator;
	}

	public VcfGenotype getVcfGenotype(int index) {
		if (vcfGenotypes == null) parseGenotypes();
		return vcfGenotypes.get(index);
	}

	public List<VcfGenotype> getVcfGenotypes() {
		if (vcfGenotypes == null) parseGenotypes();
		return vcfGenotypes;
	}

	/**
	 * Get Info type for a given ID
	 * @param id
	 * @return VcfInfoType
	 */
	public VcfInfoType getVcfInfoType(String id) {
		VcfInfo vcfInfo = vcfFileIterator.getVcfHeader().getVcfInfo(id);
		if (vcfInfo == null) return null;
		return vcfInfo.getVcfInfoType();
	}

	public boolean hasField(String filedName) {
		return vcfFileIterator.getVcfHeader().getVcfInfo(filedName) != null;
	}

	public boolean isDel() {
		return (changeType == ChangeType.DEL);
	}

	public boolean isFilterPass() {
		return filterPass.equals("PASS");
	}

	/**
	 * Is this heterozygous (based ONLY on the number of ALTs)
	 * WARINIG: You should use 'calcHetero()' method for a more precise calculation.
	 * 
	 * @return
	 */
	public boolean isHeterozygous() {
		return alts.length > 1; // More than one ALT option? => not homozygous
	}

	/**
	 * Is this homozygous (based ONLY on the number of ALTs)
	 * WARINIG: You should use 'calcHetero()' method for a more precise calculation.
	 * 
	 * @return
	 */
	public boolean isHomozygous() {
		return alts.length == 1; // Only one ALT option? => homozygous
	}

	public boolean isInDel() {
		return (changeType == ChangeType.INS) || (changeType == ChangeType.DEL);
	}

	public boolean isIns() {
		return (changeType == ChangeType.INS);
	}

	public boolean isInterval() {
		return (changeType == ChangeType.Interval);
	}

	public boolean isMixedInDel() {
		return changeType == ChangeType.MIXED;
	}

	public boolean isMnp() {
		return changeType == ChangeType.MNP;
	}

	@Override
	protected boolean isShowWarningIfParentDoesNotInclude() {
		return false;
	}

	public boolean isSnp() {
		return changeType == ChangeType.SNP;
	}

	/**
	 * Is this a change or are the ALTs actually the same as the reference
	 * @return
	 */
	public boolean isVariant() {
		for (String alt : alts)
			if (!alt.isEmpty() && !alt.equals(".") && !ref.equals(alt)) return true; // Any change option is different? => true
		return false;
	}

	@Override
	public Iterator<VcfGenotype> iterator() {
		if (vcfGenotypes == null) parseGenotypes();
		return vcfGenotypes.iterator();
	}

	/**
	 * Parse a 'line' from a 'vcfFileIterator'
	 */
	public void parse() {
		// Parse line
		String fields[] = line.split("\t", 10); // Only pare the fist 9 fields (i.e. do not parse genotypes)

		// Is line OK?
		if (fields.length >= 4) {
			// Chromosome and position. VCF files are one-base, so inOffset should be 1.
			chromosomeName = fields[0].trim();
			Chromosome chromo = vcfFileIterator.getChromosome(chromosomeName);
			parent = chromo;
			vcfFileIterator.sanityCheckChromo(chromosomeName, chromo); // Sanity check

			start = vcfFileIterator.parsePosition(vcfFileIterator.readField(fields, 1));

			// ID (e.g. might indicate dbSnp)
			id = vcfFileIterator.readField(fields, 2);

			// REF and ALT
			ref = vcfFileIterator.readField(fields, 3).toUpperCase(); // Reference and change
			end = start + ref.length();
			strand = 1;
			String altsStr = vcfFileIterator.readField(fields, 4).toUpperCase();
			parseAlts(altsStr);

			// Quality
			String qStr = vcfFileIterator.readField(fields, 5);
			if (!qStr.isEmpty()) quality = Gpr.parseDoubleSafe(qStr);
			else quality = null;

			filterPass = vcfFileIterator.readField(fields, 6); // Filter parameters

			// INFO fields
			infoStr = vcfFileIterator.readField(fields, 7);
			info = null;

			// Genotype format
			format = null;
			if (fields.length > 8) format = vcfFileIterator.readField(fields, 8); // This field is optional, So it can be null or EMPTY ('.')

			// Add genotype fields (lazy parse) 
			if (fields.length > 9) genotypeFieldsStr = fields[9];
		}
	}

	/**
	 * Parse ALT field
	 * @param altsStr
	 */
	void parseAlts(String altsStr) {
		if (altsStr.length() == 1) {
			if (altsStr.equals("A") || altsStr.equals("C") || altsStr.equals("G") || altsStr.equals("T") || altsStr.equals(".")) {
				alts = new String[1];
				alts[0] = altsStr;
			} else if (altsStr.equals("N")) { // aNy base
				alts = new String[4];
				alts[0] = "A";
				alts[1] = "C";
				alts[2] = "G";
				alts[3] = "T";
			} else if (altsStr.equals("B")) { // B: not A
				alts = new String[3];
				alts[0] = "C";
				alts[1] = "G";
				alts[2] = "T";
			} else if (altsStr.equals("D")) { // D: not C
				alts = new String[3];
				alts[0] = "A";
				alts[1] = "G";
				alts[2] = "T";
			} else if (altsStr.equals("H")) { // H: not G
				alts = new String[3];
				alts[0] = "A";
				alts[1] = "C";
				alts[2] = "T";
			} else if (altsStr.equals("V")) { // V: not T
				alts = new String[3];
				alts[0] = "A";
				alts[1] = "C";
				alts[2] = "G";
			} else if (altsStr.equals("M")) {
				alts = new String[2];
				alts[0] = "A";
				alts[1] = "C";
			} else if (altsStr.equals("R")) {
				alts = new String[2];
				alts[0] = "A";
				alts[1] = "G";
			} else if (altsStr.equals("W")) { // Weak
				alts = new String[2];
				alts[0] = "A";
				alts[1] = "T";
			} else if (altsStr.equals("S")) { // Strong
				alts = new String[2];
				alts[0] = "C";
				alts[1] = "G";
			} else if (altsStr.equals("Y")) {
				alts = new String[2];
				alts[0] = "C";
				alts[1] = "T";
			} else if (altsStr.equals("K")) {
				alts = new String[2];
				alts[0] = "G";
				alts[1] = "T";
			} else if (altsStr.equals(".")) { // No alternative (same as reference)
				alts = new String[1];
				alts[0] = ref;
			} else {
				throw new RuntimeException("WARNING: Unkown IUB code for SNP '" + altsStr + "'");
			}
		} else alts = altsStr.split(",");

		// What type of change do we have?
		int maxAltLen = Integer.MIN_VALUE, minAltLen = Integer.MAX_VALUE;
		for (int i = 0; i < alts.length; i++) {
			maxAltLen = Math.max(maxAltLen, alts[i].length());
			minAltLen = Math.min(minAltLen, alts[i].length());
		}

		// Infer change type
		if ((ref.length() == maxAltLen) && (ref.length() == minAltLen)) {
			if (ref.length() == 1) changeType = ChangeType.SNP;
			else changeType = ChangeType.MNP;
		} else if (ref.length() > minAltLen) changeType = ChangeType.DEL;
		else if (ref.length() < maxAltLen) changeType = ChangeType.INS;
		else changeType = ChangeType.MIXED;
	}

	/**
	 * Parse 'EFF' info field and get a list of effects
	 * @return
	 */
	public List<VcfEffect> parseEffects(FormatVersion formatVersion) {
		String effStr = getInfo(VcfOutputFormatter.VCF_INFO_EFF_NAME); // Get effect string from INFO field

		// Create a list of effect
		ArrayList<VcfEffect> effList = new ArrayList<VcfEffect>();
		if (effStr == null) return effList;

		// Add each effect
		String effs[] = effStr.split(",");
		for (String eff : effs) {
			VcfEffect veff = new VcfEffect(eff, formatVersion); // Create and parse this effect
			effList.add(veff);
		}
		return effList;
	}

	/**
	 * Parse GENOTPYE entries
	 */
	void parseGenotypes() {
		vcfGenotypes = new ArrayList<VcfGenotype>();

		// No genotype string? => Nothing to do
		if (genotypeFieldsStr == null) return;

		// Split genotypes and parse them
		genotypeFields = genotypeFieldsStr.split("\t");
		for (int i = 0; i < genotypeFields.length; i++) {
			String gen = genotypeFields[i];
			if (gen.equals(VcfFileIterator.MISSING)) gen = "";
			addGenotype(gen);
		}
	}

	/**
	 * Parse INFO fields
	 */
	void parseInfo() {
		// Parse info entries
		info = new HashMap<String, String>();
		for (String inf : infoStr.split(";")) {
			String vp[] = inf.split("=");

			if (vp.length > 1) info.put(vp[0], vp[1]);
			else info.put(vp[0], "true"); // A property that is present, but has no value (e.g. "INDEL")
		}
	}

	/**
	 * Create a list of seqChanges frmo this VcfEntry
	 * @return
	 */
	public List<SeqChange> seqChanges() {
		LinkedList<SeqChange> list = new LinkedList<SeqChange>();

		// Coverage
		int coverage = Gpr.parseIntSafe(getInfo("DP"));

		// Is it heterozygous, homozygous or undefined?
		Boolean isHetero = calcHetero();

		// Create one SeqChange for each alt
		for (String alt : alts) {
			Chromosome chr = (Chromosome) parent;
			SeqChange seqChange = createSeqChange(chr, start, ref, alt, strand, id, getQuality(), coverage);
			seqChange.setHeterozygous(isHetero);
			list.add(seqChange);
		}

		return list;
	}

	public void setFilterPass(String filterPass) {
		this.filterPass = filterPass;
	}

	public void setGenotypeStr(String genotypeFieldsStr) {
		this.genotypeFieldsStr = genotypeFieldsStr;
	}

	public void setLineNum(int lineNum) {
		this.lineNum = lineNum;
	}

	@Override
	public String toString() {
		boolean deleteLastTab = true;

		// Use original chromosome name or named from chromosome object
		String chr = chromosomeName != null ? chromosomeName : null;
		if (chromosomeName != null) chr = chromosomeName;
		else if (parent != null) chr = parent.getId();
		else chr = ".";

		StringBuilder sb = new StringBuilder(chr //
				+ "\t" + (start + 1) //
				+ "\t" + (id.isEmpty() ? "." : id) //
				+ "\t" + ref //
				+ "\t");

		// ALTs
		for (int i = 0; i < alts.length; i++) {
			String altStr = (alts[i].isEmpty() ? "." : alts[i]);
			sb.append(altStr + ",");
		}
		sb.deleteCharAt(sb.length() - 1); // Delete last colon

		// Quality, filter, info, format...
		sb.append("\t" + (quality != null ? quality + "" : "."));
		sb.append("\t" + (filterPass.isEmpty() ? "." : filterPass));
		sb.append("\t" + (infoStr.isEmpty() ? "." : infoStr));
		sb.append("\t");

		// Is there any 'format' field? It is optional, so it could be 'null'
		if (format != null) {
			sb.append((format.isEmpty() ? "." : format) + "\t");

			// If we have vcfGenotypes parsed, use them
			if ((vcfGenotypes != null) && !vcfGenotypes.isEmpty()) {
				for (VcfGenotype vg : vcfGenotypes)
					sb.append(vg + "\t");
			} else if (genotypeFieldsStr != null) { // If vcfGenotypes have not been parsed, use raw fields
				sb.append(genotypeFieldsStr);
				deleteLastTab = false;
			}
		}

		if (deleteLastTab) sb.deleteCharAt(sb.length() - 1); // Delete last tab

		return sb.toString();
	}
}
