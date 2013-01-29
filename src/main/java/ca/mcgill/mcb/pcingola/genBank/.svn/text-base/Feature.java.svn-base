package ca.mcgill.mcb.pcingola.genBank;

import java.util.HashMap;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * A feature in a GenBank or EMBL file
 * 
 * @author pablocingolani
 */
public class Feature {

	public enum Type {
		SOURCE, ID, CDS, GENE, MRNA, TRNA, RRNA, MISC_RNA, REPEAT_UNIT, REPEAT_REGION, MISC_FEATURE, UTR_3, UTR_5;

		/**
		 * Parse a string into a Feature.Type
		 */
		public static Feature.Type parse(String typeStr) {
			typeStr = typeStr.toUpperCase();
			typeStr = typeStr.replaceAll("[^A-Za-z0-9]", "_");

			// Some equivalences
			if (typeStr.equals("5_UTR")) return UTR_5;
			if (typeStr.equals("3_UTR")) return UTR_3;
			if (typeStr.equals("SQ")) return SOURCE;

			try {
				return Feature.Type.valueOf(typeStr);
			} catch (Exception e) {
				return null;
			}
		}

	}

	static final String FEATURE_REGEX = "/(\\S+?)=(.*)";
	static final Pattern FEATURE_PATTERN = Pattern.compile(FEATURE_REGEX);;
	public static final String COMPLEMENT_STRING = "complement";

	Type type;
	int start, end;
	HashMap<String, String> qualifiers;
	boolean complement;

	public Feature(Type type, String def) {
		this.type = type;
		qualifiers = new HashMap<String, String>();
		start = -1;
		end = -1;
		complement = false;
		parse(def);
	}

	public Feature(Type type, String def, int start, int end, boolean complement) {
		this.type = type;
		qualifiers = new HashMap<String, String>();
		this.complement = complement;

		// Assign start & end
		if (end < start) { // Order reversed? Swap them
			int tmp = end;
			end = start;
			start = tmp;
		}
		this.start = start;
		this.end = end;

		// Parse
		parse(def);

		// Sanity check
		if (start < 0) throw new RuntimeException("Feature starts with negative coordinates!\n\t" + this);
	}

	/**
	 * Get a qualifier by name
	 * 
	 * @param name
	 * @return
	 */
	public String get(String name) {
		return qualifiers.get(name);
	}

	public int getEnd() {
		return end;
	}

	public int getStart() {
		return start;
	}

	public Type getType() {
		return type;
	}

	public boolean isComplement() {
		return complement;
	}

	/**
	 * Parse definition
	 * @param def
	 */
	void parse(String def) {
		int firstLine = def.indexOf("\n");

		// Parse location (first line), if required
		if ((start < 0) && (end < 0)) {
			String loc = def.substring(0, firstLine);
			parseLocation(loc);
		}

		//---
		// Parse all other features
		//---
		def = def.substring(firstLine + 1);
		Matcher matcher = FEATURE_PATTERN.matcher(def);
		while (matcher.find()) {
			if (matcher.groupCount() >= 2) {
				String key = matcher.group(1).toLowerCase();
				String value = matcher.group(2);
				if (value.startsWith("\"") && value.endsWith("\"")) value = value.substring(1, value.length() - 1);
				qualifiers.put(key, value.trim());
			}
		}
	}

	/**
	 * Parse location
	 * @param loc
	 */
	void parseLocation(String loc) {
		loc = loc.replaceAll("[<>()]", "");
		if (loc.startsWith("complement")) {
			complement = true;
			loc = loc.substring(COMPLEMENT_STRING.length());
		}

		String se[] = loc.split("[\\.]+");
		if (se.length > 1) {
			start = Gpr.parseIntSafe(se[0]);
			end = Gpr.parseIntSafe(se[1]);
		}
	}

	/**
	 * Remove surrounding quotes from a string
	 * @param s
	 * @return
	 */
	String removeQuotes(String s) {
		if (s.startsWith("\"")) s = s.substring(1);
		if (s.endsWith("\"")) s = s.substring(0, s.length() - 1);
		return s;
	}

	@Override
	public String toString() {
		String format = "\t%-20s: \"%s\"\n";
		StringBuilder sb = new StringBuilder();

		sb.append("Feature: '" + type //
				+ "' [ " + start + ", " + end + " ]\t" //
				+ (complement ? "complement" : "") //
				+ "\n" //
		);

		for (Entry<String, String> e : qualifiers.entrySet())
			sb.append(String.format(format, e.getKey(), e.getValue()));

		return sb.toString();
	}
}
