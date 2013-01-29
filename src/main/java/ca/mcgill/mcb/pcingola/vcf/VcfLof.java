package ca.mcgill.mcb.pcingola.vcf;

import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * An 'LOF' entry in a vcf line
 * 
 * @author pablocingolani
 */
public class VcfLof {

	String geneName;
	String geneId;
	int numTranscripts;
	double percentAffected;

	/**
	 * Convert from field name to field number
	 * @param name
	 * @param formatVersion
	 * @return
	 */
	public static int fieldNum(String name) {
		int fieldNum = 0;

		if (name.equals("LOF.GENE")) return fieldNum;
		fieldNum++;

		if (name.equals("LOF.GENEID")) return fieldNum;
		fieldNum++;

		if (name.equals("LOF.NUMTR")) return fieldNum;
		fieldNum++;

		if (name.equals("LOF.PERC")) return fieldNum;
		fieldNum++;

		return -1;
	}

	public VcfLof(String lofStr) {
		parse(lofStr);
	}

	void parse(String lof) {
		String lofFields[] = lof.split("|");

		try {
			Gpr.debug("PARSE: '" + lof + "'");

			// Parse each sub field
			int index = 0;

			if ((lofFields.length > index) && !lofFields[index].isEmpty()) geneName = lofFields[index];
			index++;

			if ((lofFields.length > index) && !lofFields[index].isEmpty()) geneId = lofFields[index];
			index++;

			if ((lofFields.length > index) && !lofFields[index].isEmpty()) numTranscripts = Gpr.parseIntSafe(lofFields[index]);
			index++;

			if ((lofFields.length > index) && !lofFields[index].isEmpty()) percentAffected = Gpr.parseDoubleSafe(lofFields[index]);
			index++;

		} catch (Exception e) {
			String fields = "";
			for (int i = 0; i < lofFields.length; i++)
				fields += "\t" + i + " : '" + lofFields[i] + "'\n";
			throw new RuntimeException("Error parsing: '" + lof + "'\n" + fields, e);
		}
	}

	@Override
	public String toString() {
		return String.format("(%s|%s|%d|%.2f)", geneName, geneId, numTranscripts, percentAffected);
	}
}
