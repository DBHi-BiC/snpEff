package ca.mcgill.mcb.pcingola.interval;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.zip.GZIPOutputStream;

import ca.mcgill.mcb.pcingola.fileIterator.LineFileIterator;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * Serialize markers to (and from) file
 * 
 * Note: Marker's children are serialized first (e.g. a transcript get all 
 * 		exons serialized first).  
 * 
 * Note: Since Marker is a tree-like structure, we first load all the markers and then 
 * 		assign parents. Markers are assigned a fake parent object (MarkerParentId) 
 * 		which is later replaced by the 'real' parent.
 * 
 * Note: All 'IDs' used have not meaning outside this serialization process. 
 * 
 * @author pcingola
 */
public class MarkerSerializer {

	PrintStream outFile;
	int lineNum;
	String line;
	int parsedField;
	String fields[];
	int currId = 0;
	HashMap<Integer, Marker> markerById;
	HashMap<Marker, Integer> idByMarker;

	public MarkerSerializer() {
		markerById = new HashMap<Integer, Marker>();
		idByMarker = new HashMap<Marker, Integer>();
	}

	/**
	 * Add all data from 'genome' to markres
	 * @param markers
	 * @param genome
	 */
	protected void add(Markers markers, Genome genome) {
		markers.add(genome);

		for (Chromosome chr : genome)
			markers.add(chr);

		for (Gene g : genome.getGenes())
			markers.add(g);
	}

	protected int getIdByMarker(Marker m) {
		Integer id = idByMarker.get(m);
		if (id == null) { throw new RuntimeException("Marker has no numeric ID. \n" //
				+ "\tClass    : " + m.getClass().getSimpleName() + "\n" //
				+ "\tMarker ID: '" + m.getId() + "'\n" // 
				+ "\t" + m); }
		return id;
	}

	protected Marker getMarkerById(int id) {
		return markerById.get(id);
	}

	protected String getNextField() {
		if (fields.length <= parsedField) return "";
		return fields[parsedField++];
	}

	protected boolean getNextFieldBoolean() {
		return Gpr.parseBoolSafe(getNextField());
	}

	protected int getNextFieldInt() {
		return Gpr.parseIntSafe(getNextField());
	}

	protected Marker getNextFieldMarker() {
		return getMarkerById(getNextFieldInt());
	}

	protected Markers getNextFieldMarkers() {
		Markers markers = new Markers();
		String fieldIdsStr = getNextField();
		if (fieldIdsStr.isEmpty()) return markers;

		String fieldIds[] = fieldIdsStr.split(",");
		for (String idStr : fieldIds) {
			int id = Gpr.parseIntSafe(idStr);
			Marker m = getMarkerById(id);
			if (m != null) markers.add(m);
			else throw new RuntimeException("Marker '" + id + "' not found. This should never happen!");
		}
		return markers;
	}

	protected int getNextId() {
		return ++currId;
	}

	/**
	 * Load data from file
	 * 
	 * @param fileName
	 */
	public Markers load(String fileName) {
		//---
		// Load data from file
		//---
		LineFileIterator lfi = new LineFileIterator(fileName, true); // File is gzipped
		int lineNum = 0;
		for (String l : lfi) {
			line = l;
			parsedField = 0;
			fields = line.split("\t", -1);

			// Parse field type
			String typeStr = fields[0];
			EffectType type = EffectType.valueOf(typeStr);

			// Parse serialization id
			String idStr = fields[1];
			int id = Gpr.parseIntSafe(idStr);

			Marker m = null;
			switch (type) {
			case GENOME:
				m = new Genome();
				break;
			case CHROMOSOME:
				m = new Chromosome();
				break;
			case GENE:
				m = new Gene();
				break;
			case TRANSCRIPT:
				m = new Transcript();
				break;
			case CDS:
				m = new Cds();
				break;
			case EXON:
				m = new Exon();
				break;
			case UTR_3_PRIME:
				m = new Utr3prime();
				break;
			case UTR_5_PRIME:
				m = new Utr5prime();
				break;
			case RARE_AMINO_ACID:
				m = new RareAminoAcid();
				break;
			case SPLICE_SITE_ACCEPTOR:
				m = new SpliceSiteAcceptor();
				break;
			case SPLICE_SITE_BRANCH:
				m = new SpliceSiteBranch();
				break;
			case SPLICE_SITE_BRANCH_U12:
				m = new SpliceSiteBranchU12();
				break;
			case SPLICE_SITE_DONOR:
				m = new SpliceSiteDonor();
				break;
			default:
				throw new RuntimeException("Unimplemented for type '" + type + "'");
			}

			try {
				// Parse line
				m.serializeParse(this);
			} catch (Throwable t) {
				t.printStackTrace();
				throw new RuntimeException("Error parsing line " + (lineNum + 1) + " from file '" + fileName + "'\n\t" + line + "\n\tField [" + parsedField + "] : '" + (parsedField < fields.length ? fields[parsedField] : "-") + "'", t);
			}

			// Add to hash
			markerById.put(id, m);
			lineNum++;
		}

		//--- 
		// Assign parents
		//---
		Markers markers = new Markers();
		for (Marker m : markerById.values()) {
			// Find parent ID
			MarkerParentId mpid = (MarkerParentId) m.getParent();
			int parentId = mpid.getParentId();

			// Find and set parent
			Marker parent = getMarkerById(parentId);
			m.setParent(parent);

			// Add to markers
			markers.add(m);
		}

		return markers;
	}

	/**
	 * Save all markers
	 * @param markersCollection
	 * @return
	 */
	protected String save(Iterable<Marker> markersCollection) {
		StringBuilder idStr = new StringBuilder();
		for (Marker m : markersCollection) {
			int id = save(m);
			if (idStr.length() > 0) idStr.append(",");
			idStr.append(id);
		}
		return idStr.toString();
	}

	/**
	 * Save a marker
	 * @param m
	 */
	protected int save(Marker m) {
		if (m == null) return -1;
		if (idByMarker.containsKey(m)) return idByMarker.get(m); // Already done

		// Store already saved IDs
		int id = getNextId();
		idByMarker.put(m, id);

		// Print line
		String line = m.serializeSave(this);
		outFile.print(line + "\n");
		lineNum++;

		return id;
	}

	/**
	 * Save a genome
	 * @param fileName
	 * @param genome
	 */
	public void save(String fileName, Genome genome) {
		Markers markers = new Markers();
		add(markers, genome);
		save(fileName, markers); // Write
	}

	/**
	 * Save data to file
	 * @param fileName
	 */
	public void save(String fileName, Markers markers) {
		try {
			lineNum = 0;
			currId = 0;
			outFile = new PrintStream(new GZIPOutputStream(new FileOutputStream(fileName)));
			for (Marker m : markers)
				save(m);
			outFile.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Save all markers in snpEffectPredictor
	 * @param fileName
	 * @param genome
	 */
	public void save(String fileName, SnpEffectPredictor snpEffectPredictor) {
		Markers markers = new Markers();
		add(markers, snpEffectPredictor.getGenome());
		markers.add(snpEffectPredictor.getMarkers());
		save(fileName, markers);
	}
}
