package ca.mcgill.mcb.pcingola.ped;

/**
 * A family of PedEntries
 * 
 * @author pcingola
 */
public class PedFamily extends PedPedigree {

	String familyId = null;

	public PedFamily() {
		super();
	}

	/**
	 * Add an entry t this family
	 * @param pedEntry
	 */
	@Override
	public void add(PedEntry pedEntry) {
		if( (familyId != null) && (!familyId.equals(pedEntry.getFamilyId())) ) throw new RuntimeException("Cannot add memeber to family. Family IDs do not match: '" + familyId + "' vs '" + pedEntry.getFamilyId() + "'");
		pedById.put(pedEntry.getId(), pedEntry);
	}

}
