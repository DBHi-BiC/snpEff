package ca.mcgill.mcb.pcingola.ped;

import java.util.Collection;
import java.util.Iterator;

/**
 * An entry in a PED table.
 * I.e. a line in a PED file (PLINK)
 * 
 * @author pcingola
 */
public class PedEntry implements Iterable<PedGenotype> {

	PlinkMap plinkMap;
	String familyId, id, fatherId, motherId;
	Sex sex;
	String genotypes[];
	String phenotype;

	public PedEntry(PlinkMap plinkMap, String familyId, String id, String fatherId, String motherId, Sex sex, String phenotype, String genotypes[]) {
		this.plinkMap = plinkMap;
		this.familyId = familyId;
		this.id = id;
		this.fatherId = fatherId;
		this.motherId = motherId;
		this.sex = sex;
		this.phenotype = phenotype;
		this.genotypes = genotypes;
	}

	/**
	 * Number of phenotypes available 
	 * @return
	 */
	public int countGenotypes() {
		int count = 0;
		for (PedGenotype gen : this)
			if (gen.isValid()) count++;
		return count;
	}

	String genoStr(String geno) {
		if (geno.equals("x") || geno.equals("0")) return "";
		return geno;
	}

	public String getFamilyId() {
		return familyId;
	}

	public String getFatherId() {
		return fatherId;
	}

	/**
	 * Get genotype
	 * WARNING: Empty string means that no genotpye is available
	 * 
	 * @param idx
	 * @return
	 */
	public PedGenotype getGenotype(int idx) {
		String geno[] = new String[2];
		geno[0] = genoStr(genotypes[idx * 2]);
		geno[1] = genoStr(genotypes[idx * 2 + 1]);
		return new PedGenotype(geno, plinkMap.getChrName(idx), plinkMap.getPosition(idx));
	}

	/**
	 * Get phenotype by String ID
	 * @param idStr
	 * @return
	 */
	public PedGenotype getGenotype(String idStr) {
		Integer idxInt = plinkMap.getGenotypeNames(idStr);
		if (idxInt == null) return null;
		return getGenotype(idxInt);
	}

	/**
	 * Get all genotype names 
	 * WARNING: the returned string collection is unsorted!
	 * @return
	 */
	public Collection<String> getGenotypeNames() {
		return plinkMap.getGenotypeNames();
	}

	public String[] getGenotypes() {
		return genotypes;
	}

	public String getId() {
		return id;
	}

	public String getMotherId() {
		return motherId;
	}

	public String getPhenotype() {
		return phenotype;
	}

	public Sex getSex() {
		return sex;
	}

	@Override
	public Iterator<PedGenotype> iterator() {
		return new Iterator<PedGenotype>() {

			int i = 0;
			int max = size();

			@Override
			public boolean hasNext() {
				return i < max;
			}

			@Override
			public PedGenotype next() {
				return getGenotype(i++);
			}

			@Override
			public void remove() {
				throw new UnsupportedOperationException("Unsupported! Go away and die!");
			}
		};
	}

	/**
	 * Number of phenotypes
	 * @return
	 */
	public int size() {
		return genotypes.length / 2;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		sb.append(familyId + "\t");
		sb.append(id + "\t");
		sb.append(fatherId + "\t");
		sb.append(motherId + "\t");
		sb.append(sex + "\t");

		int len = genotypes.length / 2;
		for (int i = 0; i < len; i += 2)
			sb.append(genotypes[i] + "/" + genotypes[i + 1] + "\t");

		return sb.toString();
	}
}
