package ca.mcgill.mcb.pcingola.vcf;


/**
 * An 'NMD' entry in a vcf line
 * 
 * @author pablocingolani
 */
public class VcfNmd extends VcfLof {

	/**
	 * Convert from field name to field number
	 * @param name
	 * @param formatVersion
	 * @return
	 */
	public static int fieldNum(String name) {
		int fieldNum = 0;

		if (name.equals("NMD.GENE")) return fieldNum;
		fieldNum++;

		if (name.equals("NMD.GENEID")) return fieldNum;
		fieldNum++;

		if (name.equals("NMD.NUMTR")) return fieldNum;
		fieldNum++;

		if (name.equals("NMD.PERC")) return fieldNum;
		fieldNum++;

		return -1;
	}

	public VcfNmd(String nmdStr) {
		super(nmdStr);
	}

}
