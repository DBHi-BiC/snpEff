package ca.mcgill.mcb.pcingola.ped;

import java.util.Collection;
import java.util.HashMap;

import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * PLINK MAP file
 * 
 * References: http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml
 * 
 * @author pcingola
 */
public class PlinkMap {

	public static boolean debug = false;

	HashMap<String, Integer> genotypeNames; // TODO: This should NOT be called genotypeNames!!!
	String chrNames[];
	int positions[];
	String ids[];

	public PlinkMap() {
		genotypeNames = new HashMap<String, Integer>();
	}

	public String getChrName(int idx) {
		return chrNames[idx];
	}

	public String getId(int idx) {
		return ids[idx];
	}

	public Collection<String> getGenotypeNames() {
		return genotypeNames.keySet();
	}

	public Integer getGenotypeNames(String idStr) {
		return genotypeNames.get(idStr);
	}

	public int getPosition(int idx) {
		return positions[idx];
	}

	/**
	 * Reads MAP file
	 * 
	 * MAP file format (http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml)
	 * 
	 *  Space separated or tab columns:
	 *		chromosome (1-22, X, Y or 0 if unplaced)
	 *		rs# or snp identifier
	 *		Genetic distance (morgans)
	 *		Base-pair position (bp units)
	 *            
	 *            
	 * @param dataFileName
	 */
	public void read(String mapFileName) {
		String cols = Gpr.readFile(mapFileName);
		String lines[] = cols.split("\n");

		positions = new int[lines.length];
		chrNames = new String[lines.length];
		ids = new String[lines.length];

		int lineNum = 0;
		for (String line : lines) {
			String fields[] = line.split("\\s");

			chrNames[lineNum] = fields[0];
			ids[lineNum] = fields[1];
			positions[lineNum] = Gpr.parseIntSafe(fields[fields.length - 1]);

			if (!ids[lineNum].isEmpty()) {
				// Is it duplicate?
				if (genotypeNames.containsKey(ids)) throw new RuntimeException("Duplicate ID '" + ids[lineNum] + "'. File '" + mapFileName + "', line '" + (lineNum + 1) + "'");
				genotypeNames.put(ids[lineNum], lineNum);
				if (debug) Gpr.debug("genotypeNames.put(" + ids[lineNum] + ", " + lineNum + ")");
			}
			lineNum++;
		}
	}

	public int size() {
		return positions.length;
	}

}
