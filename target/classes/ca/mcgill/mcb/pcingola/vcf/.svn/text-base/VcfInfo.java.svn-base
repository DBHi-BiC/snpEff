package ca.mcgill.mcb.pcingola.vcf;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * Represents a info elements in a VCF file
 * 
 * References: http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
 * 
 * INFO fields should be described as follows (all keys are required):
 * 		##INFO=<ID=ID,Number=number,Type=type,Description=�description�>
 * 
 * 		Possible Types for INFO fields are: Integer, Float, Flag, Character, and String.
 * 
 * 		The Number entry is an Integer that describes the number of values that 
 * 		can be included with the INFO field. For example, if the INFO field contains 
 * 		a single number, then this value should be 1; if the INFO field describes a 
 * 		pair of numbers, then this value should be 2 and so on. If the field has one 
 * 		value per alternate allele then this value should be 'A'; if the field has 
 * 		one value for each possible genotype (more relevant to the FORMAT tags) then 
 * 		this value should be 'G'.  If the number of possible values varies, is unknown, 
 * 		or is unbounded, then this value should be '.'. The 'Flag' type indicates that 
 * 		the INFO field does not contain a Value entry, and hence the Number should be 0 in 
 * 		this case. The Description value must be surrounded by double-quotes. Double-quote 
 * 		character can be escaped with backslash (\") and backslash as \\.
 * 
 * @author pablocingolani
 *
 */
public class VcfInfo {

	String line;
	String id;
	VcfInfoType vcfInfoType;
	int number;
	boolean onePerAllele;
	boolean onePerGenotype;
	String description;

	/**
	 * Constructor using a "##INFO" line from a VCF file 
	 * @param line
	 */
	public VcfInfo(String line) {
		// Is this an Info line?
		if (line.startsWith("##INFO=") || line.startsWith("##FORMAT=")) {
			// Remove all trailing '\n'
			while (line.endsWith("\n"))
				line = line.substring(0, line.length() - 1);
			this.line = line;

			int start = line.indexOf('<');
			int end = line.lastIndexOf('>');
			String params = line.substring(start + 1, end);

			// Find ID
			Pattern pattern = Pattern.compile("ID=([^,]+),");
			Matcher matcher = pattern.matcher(params);
			if (matcher.find()) id = matcher.group(1);
			else throw new RuntimeException("Cannot find 'ID' in info line: '" + line + "'");

			// Find Number
			pattern = Pattern.compile("Number=([^,]+),");
			matcher = pattern.matcher(params);
			if (matcher.find()) parseNumber(matcher.group(1));
			else throw new RuntimeException("Cannot find 'Number' in info line: '" + line + "'");

			// Find type
			pattern = Pattern.compile("Type=([^,]+),");
			matcher = pattern.matcher(params);
			if (matcher.find()) vcfInfoType = VcfInfoType.parse(matcher.group(1).toUpperCase());
			else throw new RuntimeException("Cannot find 'Type' in info line: '" + line + "'");

			// Find description
			pattern = Pattern.compile("Description=\\\"(.+)\\\"");
			matcher = pattern.matcher(params);
			if (matcher.find()) description = matcher.group(1);
			else throw new RuntimeException("Cannot find 'Description' in info line: '" + line + "'");

		} else throw new RuntimeException("Line provided is not an INFO definition: '" + line + "'");
	}

	public VcfInfo(String id, VcfInfoType vcfInfoType, String number, String description) {
		this.id = id;
		this.vcfInfoType = vcfInfoType;
		this.description = description;
		parseNumber(number);
	}

	public String getDescription() {
		return description;
	}

	public String getId() {
		return id;
	}

	public int getNumber() {
		return number;
	}

	public VcfInfoType getVcfInfoType() {
		return vcfInfoType;
	}

	public boolean isOnePerAllele() {
		return onePerAllele;
	}

	public boolean isOnePerGenotype() {
		return onePerGenotype;
	}

	void parseNumber(String number) {
		// Parse number field
		if (number.equals("A")) onePerAllele = true;
		else if (number.equals("G")) onePerGenotype = true;
		else this.number = Gpr.parseIntSafe(number);
	}

	@Override
	public String toString() {
		if (line != null) return line;

		return "##INFO=<ID=" + id//
				+ ",Number=" + (onePerAllele ? "A" : (onePerGenotype ? "G" : number)) //
				+ ",Type=" + vcfInfoType //
				+ ",Description=\"" + description + "\"" //
				+ ">" //
		;
	}

}
