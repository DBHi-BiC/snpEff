package ca.mcgill.mcb.pcingola.ped;

import java.io.IOException;

import ca.mcgill.mcb.pcingola.fileIterator.FileIterator;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * PED file iterator (PED file from PLINK)
 * 
 * Reference: http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml
 * 
 * @author pcingola
 */
public class PedFileIterator extends FileIterator<PedEntry> {

	PlinkMap plinkMap;

	public PedFileIterator(String fileName, String mapFileName) {
		super(fileName);
		plinkMap = new PlinkMap();
		plinkMap.read(mapFileName);
	}

	public PlinkMap getPlinkMap() {
		return plinkMap;
	}

	/**
	 * Parse one line
	 * @param line
	 * @return
	 */
	PedEntry parseLine(String line) {
		String fields[] = line.split("\\s");

		int fieldNum = 0;
		String familyId = fields[fieldNum++];
		String id = fields[fieldNum++];
		String fatherId = fields[fieldNum++];
		String motherId = fields[fieldNum++];

		// Parse sex field
		int sexnum = Gpr.parseIntSafe(fields[fieldNum++]);
		Sex sex = Sex.Unknown;
		if (sexnum == 1) sex = Sex.Male;
		else if (sexnum == 2) sex = Sex.Female;

		String phenotype = fields[fieldNum++];
		if (phenotype.equals("0") || phenotype.equals("-9")) phenotype = ""; // These mean 'missing'

		// Phenotypes
		String genotypes[] = new String[fields.length - fieldNum];
		for (int j = 0, i = fieldNum; i < fields.length; i++, j++)
			genotypes[j] = fields[i];

		return new PedEntry(plinkMap, familyId, id, fatherId, motherId, sex, phenotype, genotypes);
	}

	@Override
	protected PedEntry readNext() {
		try {
			if (reader.ready()) {
				line = reader.readLine(); // Read a line (only if needed)
				if (line != null) {
					lineNum++;

					PedEntry pedEntry = parseLine(line);
					if (pedEntry != null) return pedEntry;
				}
			}
		} catch (IOException e) {
			return null;
		}

		return null;
	}
}
