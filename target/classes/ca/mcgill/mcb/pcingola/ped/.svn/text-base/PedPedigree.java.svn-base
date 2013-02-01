package ca.mcgill.mcb.pcingola.ped;

import java.util.Collection;
import java.util.HashMap;
import java.util.Set;

import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * A pedigree of PedEntries
 * 
 * @author pcingola
 */
public class PedPedigree {

	boolean verbose = false;
	HashMap<String, PedEntry> pedById = new HashMap<String, PedEntry>();
	PlinkMap plinkMap;

	public PedPedigree() {
		pedById = new HashMap<String, PedEntry>();
	}

	/**
	 * Add an entry t this family
	 * @param pedEntry
	 */
	public void add(PedEntry pedEntry) {
		pedById.put(pedEntry.getId(), pedEntry);
	}

	public PedEntry get(String id) {
		return pedById.get(id);
	}

	public PlinkMap getPlinkMap() {
		return plinkMap;
	}

	public Set<String> keySet() {
		return pedById.keySet();
	}

	/**
	 * Load a pedigree from a file.
	 * Actually reads two files: PED and MAP files
	 * 
	 * @param pedFileName
	 */
	public void load(String pedFileName) {
		String pedBaseFileName = Gpr.removeExt(pedFileName);
		String mapFile = pedBaseFileName + ".map";

		PedFileIterator pedFile = new PedFileIterator(pedFileName, mapFile);

		// Load all entries for this family
		int count = 1;
		for (PedEntry pe : pedFile) {
			if (verbose) Gpr.showMarkStderr(count++, 1);
			add(pe);
		}

		plinkMap = pedFile.getPlinkMap();
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	public int size() {
		return pedById.size();
	}

	public Collection<PedEntry> values() {
		return pedById.values();
	}

}
