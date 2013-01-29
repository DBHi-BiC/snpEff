package ca.mcgill.mcb.pcingola.genBank;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedList;

import ca.mcgill.mcb.pcingola.fileIterator.LineFileIterator;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * A class representing a set of features
 * 
 * References: http://www.ebi.ac.uk/embl/Documentation/FT_definitions/feature_table.html
 * 
 * @author pablocingolani
 */
public abstract class Features implements Iterable<Feature> {

	public static final int MAX_LEN_TO_SHOW = 200;
	public static final String COMPLEMENT = "complement";
	public static final String JOIN = "join";

	String locusName, moleculeType, shape, division, date;
	int sequenceLength;
	String definition = "";
	String accession = "";
	String version = "";
	String keywords = "";
	String source = "";
	String organism = "";
	StringBuffer featuresStr;
	StringBuffer sequence;
	ArrayList<Feature> features;
	ArrayList<StringBuffer> references;
	LineFileIterator lineFileIterator;

	/**
	 * Create features from a file
	 * @param fileName
	 */
	public Features(LineFileIterator lineFileIterator) {
		references = new ArrayList<StringBuffer>();
		featuresStr = new StringBuffer();
		sequence = new StringBuffer();
		features = new ArrayList<Feature>();
		this.lineFileIterator = lineFileIterator;
		readFile();
	}

	/**
	 * Create features from a file
	 * @param fileName
	 */
	public Features(String fileName) {
		references = new ArrayList<StringBuffer>();
		featuresStr = new StringBuffer();
		sequence = new StringBuffer();
		features = new ArrayList<Feature>();
		open(fileName);
		readFile();
	}

	/**
	 * Create and add a feature
	 * @param name
	 * @param values
	 */
	void addFeature(String typeStr, StringBuilder values) {
		Feature.Type type = Feature.Type.parse(typeStr);
		if (type == null) {
			Gpr.debug("WARNING: Unknown feature '" + typeStr + "', not added.");
			return;
		}

		try {
			Collection<Feature> newFeatures = featureFactory(type, values.toString()); // Create new features
			features.addAll(newFeatures); // Add all features to the list
		} catch (Exception e) {
			Gpr.debug("Error parsing feature type '" + typeStr + "' -> '" + type + "':\n" + values);
			throw new RuntimeException(e);
		}
	}

	/**
	 * Create features from a 'type' and 'values'
	 * 
	 * Note: Most of the times, only one feature is created. When the coordinates 
	 * are 'join(....)' then we must create multiple features, since each interval 
	 * is will be mapped to a different Marker object.
	 * 
	 * 
	 * @param typeStr
	 * @param values
	 * @return
	 */
	Collection<Feature> featureFactory(Feature.Type type, String def) {
		boolean complement = false;
		LinkedList<Feature> features = new LinkedList<Feature>();

		// Get first line (location)
		int firstLine = def.indexOf("\n");
		if (def.indexOf("SPAC3C7.03c") > 0) {
			Gpr.debug("def=|" + def + "|\tfirstLine: " + firstLine);
		}

		String locStr = (firstLine >= 0) ? def.substring(0, firstLine) : def;

		// Get rid of 'join' and 'complement' strings
		String locStrPrev = "";
		while (!locStr.equals(locStrPrev)) {
			locStrPrev = locStr;

			// Is it a complement?
			if (locStr.startsWith(COMPLEMENT)) {
				complement = true;
				locStr = removeStartStr(locStr, COMPLEMENT);
			}

			locStr = removeStartStr(locStr, JOIN); // Remove 'join', if any
		}

		// Split multiple locations?
		String locs[] = locStr.split(",");
		for (String loc : locs) {
			// Get rid of 'join' and 'complement' strings
			String locPrev = "";
			while (!loc.equals(locPrev)) {
				locPrev = loc;

				// Is it a complement?
				if (loc.startsWith(COMPLEMENT)) {
					complement = true;
					loc = removeStartStr(loc, COMPLEMENT);
				}

				loc = removeStartStr(loc, JOIN); // Remove 'join', if any
			}

			// Remove other characters
			loc = loc.replaceAll("[<>()]", "");

			// Calculate start & end coordinates
			String startEnd[] = loc.split("[\\.]+");

			int start, end;
			if (startEnd.length == 2) {
				start = Gpr.parseIntSafe(startEnd[0]);
				end = Gpr.parseIntSafe(startEnd[1]);
			} else if (startEnd.length == 1) {
				start = Gpr.parseIntSafe(startEnd[0]);
				end = start;
			} else throw new RuntimeException("Cannot calculate start & end coordinates: '" + loc + "'");

			// Create feature
			Feature feature = new Feature(type, def, start, end, complement);
			features.add(feature);
		}

		return features;
	}

	public String getAccession() {
		return accession;
	}

	public String getDate() {
		return date;
	}

	public String getDefinition() {
		return definition;
	}

	public String getDivision() {
		return division;
	}

	public ArrayList<Feature> getFeatures() {
		return features;
	}

	public String getKeywords() {
		return keywords;
	}

	public String getLocusName() {
		return locusName;
	}

	public String getMoleculeType() {
		return moleculeType;
	}

	public String getOrganism() {
		return organism;
	}

	public ArrayList<StringBuffer> getReferences() {
		return references;
	}

	public String getSequence() {
		return sequence.toString();
	}

	public int getSequenceLength() {
		return sequenceLength;
	}

	public String getShape() {
		return shape;
	}

	public String getSource() {
		return source;
	}

	public String getVersion() {
		return version;
	}

	public boolean isEmpty() {
		return features.isEmpty();
	}

	/**
	 * Is there a new feature in this line?
	 * @param line
	 * @return
	 */
	protected abstract boolean isNewFeature(String line);

	@Override
	public Iterator<Feature> iterator() {
		return features.iterator();
	}

	/**
	 * Open a file
	 * @param fileName
	 */
	protected void open(String fileName) {
		if (!Gpr.canRead(fileName)) throw new RuntimeException("Cannot read file '" + fileName + "'");
		if (lineFileIterator != null) lineFileIterator.close();
		lineFileIterator = new LineFileIterator(fileName);
	}

	/**
	 * Parse features
	 */
	protected void parseFeatures() {
		// Empty?
		if (featuresStr.length() <= 0) return;

		String type = null;
		String value = "";
		StringBuilder values = new StringBuilder();
		for (String line : featuresStr.toString().split("\n")) {
			// Feature start
			if (isNewFeature(line)) {
				String kv[] = line.trim().split(" ", 2);
				if (kv.length > 1) {
					// Previous feature data is available? => Add it
					if (type != null) {
						addFeature(type, values);
						values = new StringBuilder();
					}

					// Parse name
					type = kv[0];
					value = kv[1].trim();
				}
			} else value = line.trim();

			if (value.startsWith("/")) values.append("\n");
			values.append(value);
		}

		// Add last feature
		addFeature(type, values);
	}

	/**
	 * Load and parse the contents of a data file previously opened by 'open()' method.
	 */
	protected abstract void readFile();

	/**
	 * Remove start string
	 * @param str
	 * @param startStr
	 * @return
	 */
	String removeStartStr(String str, String startStr) {
		if (str.startsWith(startStr)) return str.substring(startStr.length() + 1, str.length());
		return str;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("Name            : " + locusName + "\n");
		sb.append("Sequence length : " + sequence.length() + "\n");

		// Show references
		for (StringBuffer refsb : references) {
			sb.append("Reference       :\n");
			for (String l : refsb.toString().split("\n"))
				sb.append("                 " + l + "\n");
		}

		// Show (part of) sequence
		if (sequence.length() <= MAX_LEN_TO_SHOW) sb.append("Sequence        : " + sequence + "\n");
		else sb.append("Sequence        : " + sequence.substring(0, MAX_LEN_TO_SHOW) + "..." + "\n");

		// Show all features
		for (Feature f : features)
			sb.append(f);

		return sb.toString();
	}
}
