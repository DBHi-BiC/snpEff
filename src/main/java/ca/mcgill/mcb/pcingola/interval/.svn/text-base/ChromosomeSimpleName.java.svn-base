package ca.mcgill.mcb.pcingola.interval;

import java.util.HashMap;

/**
 * Convert chromosome names to simple names
 * @author pcingola
 */
public class ChromosomeSimpleName {

	public static final String CHROMO_PREFIX[] = { "chromosome", "chromo", "chr", "group", "scaffold", "contig", "supercontig", "supercont" }; // Must be lower case (see method)
	private static ChromosomeSimpleName instance = new ChromosomeSimpleName();

	/**
	 * Get a simple name for the chromosome
	 * @param chrName
	 * @return
	 */
	public static String get(String chrName) {
		return instance.simpleNameCache(chrName);
	}

	private final HashMap<String, String> map;

	private ChromosomeSimpleName() {
		map = new HashMap<String, String>();
	}

	/**
	 * Simplify chromosome name
	 * @param chr
	 * @return
	 */
	protected String simpleName(String chr) {
		if (chr == null) return "";
		chr = chr.trim();
		String chName = chr.toLowerCase();

		// Remove any prefix string
		for (String prefix : CHROMO_PREFIX) {
			if (chName.startsWith(prefix + ":")) return chr.substring(prefix.length() + 1);
			if (chName.startsWith(prefix + "_")) return chr.substring(prefix.length() + 1);
			if (chName.startsWith(prefix + "-")) return chr.substring(prefix.length() + 1);
			if (chName.startsWith(prefix)) return chr.substring(prefix.length());
		}

		return chr;
	}

	/**
	 * Query cache before simplifying name
	 * @param chrName
	 * @return
	 */
	protected String simpleNameCache(String chrName) {
		String chr = map.get(chrName);
		if (chr == null) {
			chr = simpleName(chrName);
			map.put(chrName, chr);
		}
		return chr;
	}

}
