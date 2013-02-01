package ca.mcgill.mcb.pcingola.genBank;

import ca.mcgill.mcb.pcingola.fileIterator.LineFileIterator;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * A class representing the same data as a GenBank file (a 'GB' file)
 * 
 * References: http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord
 * 
 * @author pablocingolani
 */
public class GenBank extends Features {

	public static final int FEATURE_NAME_FIELD_LEN = 20;

	public GenBank(LineFileIterator lineFileIterator) {
		super(lineFileIterator);
	}

	/**
	 * Create a Genbank record from a 'GB' file
	 * @param fileName
	 */
	public GenBank(String fileName) {
		super(fileName);
	}

	/**
	 * Has this line a new feature?
	 * @param line
	 * @return
	 */
	@Override
	protected boolean isNewFeature(String line) {
		return !line.substring(0, FEATURE_NAME_FIELD_LEN).trim().isEmpty(); // Feature name should be within the first 20 characters
	}

	/**
	 * Parse a feature line
	 * @param name
	 * @param value
	 * @param fieldLineNum
	 */
	protected void parseFieldLine(String name, String valueOri, int fieldLineNum) {
		String value = valueOri.trim();

		if (name.equals("LOCUS")) {
			String subfields[] = value.split(" ");
			locusName = subfields[0];
			if (subfields.length > 1) sequenceLength = Gpr.parseIntSafe(subfields[1]);
			if (subfields.length > 2) moleculeType = subfields[2];
			if (subfields.length > 3) shape = subfields[3];
			if (subfields.length > 4) division = subfields[4];
			if (subfields.length > 5) date = subfields[5];
		} else if (name.equals("DEFINITION")) {
			definition += value;
		} else if (name.equals("ACCESSION")) {
			accession += value;
		} else if (name.equals("VERSION")) {
			version += value;
		} else if (name.equals("KEYWORDS")) {
			keywords += value;
		} else if (name.equals("SOURCE")) {
			source += value;
		} else if (name.equals("REFERENCE")) {
			if (fieldLineNum == 0) references.add(new StringBuffer());
			references.get(references.size() - 1).append(value + "\n");
		} else if (name.equals("FEATURES")) {
			if (fieldLineNum > 0) featuresStr.append(valueOri + "\n"); // We need all spaces preserved for this field
		} else if (name.equals("ORIGIN")) {
			String seq[] = value.split(" ", 2);

			// First line might be empty
			if (seq.length > 1) {
				String s = seq[1].replaceAll("\\s", ""); // Remove all spaces
				sequence.append(s);
			}
		} else System.err.println("Ignored feature '" + name + "'");;
	}

	/**
	 * Load and parse the contents of a data file
	 * @param fileName
	 */
	@Override
	public void readFile() {
		int fieldLineNum = 0;
		String name = null;
		String value = "";

		// Read file
		for (String line : lineFileIterator) {
			// End of current 'chromosome'?
			if (line.equals("//")) break;

			value = line;

			// Field start
			if (!line.startsWith(" ")) {
				String kv[] = line.split(" ", 2);
				if (kv.length > 1) {
					name = kv[0];
					value = kv[1];
					fieldLineNum = 0;
				}
			}

			// Parse field
			parseFieldLine(name, value, fieldLineNum);
			fieldLineNum++;
		}

		// All features are loaded. We can parse them now
		parseFeatures();
	}

}