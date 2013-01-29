package ca.mcgill.mcb.pcingola;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import ca.mcgill.mcb.pcingola.collections.AutoHashMap;
import ca.mcgill.mcb.pcingola.fileIterator.RegulationBedFileIterator;
import ca.mcgill.mcb.pcingola.interval.Regulation;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Create a regulation consensus from multiple BED files
 * 
 * @author pcingola
 *
 */
public class RegulationConsensusMultipleBed {

	boolean verbose = false;
	AutoHashMap<String, ArrayList<String>> filesByCellType;
	HashMap<String, String> epiMarkByFile;

	String outDir;
	String inDir;
	String cellType = null;

	public RegulationConsensusMultipleBed(String inDir, String outDir) {
		this.inDir = inDir;
		this.outDir = outDir;
	}

	/**
	 * Create a consensus for each cell type
	 */
	RegulationFileConsensus consensus(RegulationFileConsensus regCons, String cellType) {
		if (verbose) Timer.showStdErr("Creating consensus for cell type: " + cellType);
		ArrayList<Regulation> regs = regCons.getRegulationList(cellType);

		if (verbose) Timer.showStdErr("Sorting: " + cellType + "\t, size: " + regs.size());
		Collections.sort(regs);

		if (verbose) Timer.showStdErr("Adding to final consensus");
		RegulationFileConsensus consCellType = new RegulationFileConsensus(false);
		for (Regulation reg : regs)
			consCellType.consensus(reg);

		regs = consCellType.getRegulationList(cellType);
		if (verbose) Timer.showStdErr("Final consensus for cell type: " + cellType + "\t, size: " + regs.size());

		return consCellType;
	}

	/**
	 * Create a consensus for each cell type
	 */
	void consensusByCellType() {
		for (String cellType : filesByCellType.keySet()) {
			ArrayList<String> bedFiles = filesByCellType.getOrCreate(cellType);
			if (verbose) Timer.showStdErr("Creating consensus for cellType '" + cellType + "', files: " + bedFiles);

			RegulationFileConsensus regCons = readBeds(bedFiles, cellType);
			regCons = consensus(regCons, cellType); // Create a new consensus from all previous consensus (yep, this is confusing...)
			regCons.save(outDir); // Save to output dir
		}
	}

	/** 
	 * Read all BED files and create a regulation consensus
	 * 
	 * @param bedFiles
	 * @return
	 */
	RegulationFileConsensus readBeds(List<String> bedFiles, String cellType) {
		RegulationFileConsensus regCons = new RegulationFileConsensus(true);

		// Read all files, creating a consensus for each
		for (String bedFile : bedFiles) {
			if (verbose) Timer.showStdErr("Reading file '" + bedFile + "'");

			String epigeneticMark = epiMarkByFile.get(bedFile);
			RegulationBedFileIterator regFile = new RegulationBedFileIterator(bedFile, epigeneticMark, cellType);
			regCons.readFile(regFile); // Add consensus 
		}

		return regCons;
	}

	/**
	 * Read all file names from a dir and looks for BED files
	 * BED file name format: regulatory.CellType.EpigeneticMark.bed[.gz]
	 * 
	 * @param dirName
	 */
	void readDir(String dirName) {
		filesByCellType = new AutoHashMap<String, ArrayList<String>>(new ArrayList<String>());
		epiMarkByFile = new HashMap<String, String>();

		// Read all files
		int countFiles = 0;
		File dir = new File(dirName);
		for (String file : dir.list()) {
			if (file.indexOf(".bed") > 0) { // Is is a 'bed' file?
				String parts[] = file.split("\\."); // Format: regulatory.CellType.EpigeneticMark.bed[.gz]

				if (parts.length > 2) {
					String cellType = parts[1];
					String epigeneticMark = parts[2];

					if ((this.cellType == null) || this.cellType.equals(cellType)) {
						// Add epigenetic mark to hash
						String filePath = dirName + "/" + file;
						epiMarkByFile.put(filePath, epigeneticMark);

						// Add file to hash
						ArrayList<String> files = filesByCellType.getOrCreate(cellType);
						files.add(filePath);

						countFiles++;
					}
				}
			}
		}

		if (verbose) Timer.showStdErr("Directory has " + countFiles + " bed files and " + filesByCellType.size() + " cell types");
	}

	/**
	 * Main method: Creates the consensus from multiple files
	 */
	public void run() {
		readDir(inDir);
		consensusByCellType();
		if (verbose) Timer.showStdErr("Done");
	}

	public void setCellType(String cellType) {
		this.cellType = cellType;
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}
}
