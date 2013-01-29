package ca.mcgill.mcb.pcingola;

import java.util.HashMap;

import ca.mcgill.mcb.pcingola.fastq.Fastq;
import ca.mcgill.mcb.pcingola.fileIterator.FastqFileIterator;
import ca.mcgill.mcb.pcingola.fileIterator.LineFileIterator;
import ca.mcgill.mcb.pcingola.fileIterator.NeedlemanWunschOverlap;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.GprSeq;
import ca.mcgill.mcb.pcingola.util.Timer;

public class ZzzPacBIo {

	String fastq;
	String putativeOverlaps;
	HashMap<String, String> seqByName;

	public static void main(String[] args) {
		String fastq = Gpr.HOME + "/scyth/pacbio/SCYTH-8Kb_0.38nM_bind-1-2_G02_1.fastq";
		String putativeOverlaps = Gpr.HOME + "/scyth/pacbio/SCYTH-8Kb_0.38nM_bind-1-2_G02_1.subreads.100.sort.filter.overlaps.txt";

		NeedlemanWunschOverlap nwo = new NeedlemanWunschOverlap("ACGT", "AC");
		String res = nwo.align();
		System.out.println("RESULT:\n" + res);
		System.out.println(nwo);

		nwo = new NeedlemanWunschOverlap("CGGCTTATAGAAGAAGGGCTGGTGTCT", "CTGGTGTCTGCCGCGAATGAGAGCATGG");
		res = nwo.align();
		System.out.println("RESULT:\n" + res);
		System.out.println(nwo);

		nwo = new NeedlemanWunschOverlap("CTGGTGTCTGCCGCGAATGAGAGCATGG", "CGGCTTATAGAAGAAGGGCTGGTGTCT");
		res = nwo.align();
		System.out.println("RESULT:\n" + res);
		System.out.println(nwo);

		ZzzPacBIo zzz = new ZzzPacBIo(fastq, putativeOverlaps);
		zzz.readSeqs();
		zzz.readOverlaps();
	}

	public ZzzPacBIo(String fastq, String putativeOverlaps) {
		this.fastq = fastq;
		this.putativeOverlaps = putativeOverlaps;
		seqByName = new HashMap<String, String>();
	}

	/**
	 * Overlap two sequences
	 * @param fq1
	 * @param fq2
	 */
	int overlap(String readName1, String readName2, boolean rwc) {
		Timer.showStdErr("Align :" + readName1 + " , " + readName2);
		String seq1 = seqByName.get(readName1);
		String seq2 = seqByName.get(readName2);

		if (rwc) seq2 = GprSeq.reverseWc(seq2);

		NeedlemanWunschOverlap nwo = new NeedlemanWunschOverlap(seq1, seq2);
		int score = nwo.calcAlignmentScore();

		if (score > 90) System.out.println("MATCH:" + score + "\n" + nwo.calcAlignment());
		return score;
	}

	/**
	 * Read overlap file
	 */
	void readOverlaps() {
		Timer.showStdErr("Reading overlaps file: " + putativeOverlaps);
		LineFileIterator lfi = new LineFileIterator(putativeOverlaps);
		for (String line : lfi) {
			String fields[] = line.split("\t");
			String readName1 = fields[0];

			System.out.println(line);
			System.out.println(readName1);
			for (int i = 1; i < fields.length; i++) {
				String sf[] = fields[i].split(";");
				String readName2 = sf[0]; // Read name

				sf = sf[1].split("/");
				int positiveMaps = Gpr.parseIntSafe(sf[0]); // Positive strand 
				int negativeMaps = sf.length > 1 ? Gpr.parseIntSafe(sf[1]) : 0; // Negative strand

				int scorePos = 0, scoreNeg = 0;
				if (positiveMaps > 0) scorePos = overlap(readName1, readName2, false);
				if (negativeMaps > 0) scoreNeg = overlap(readName1, readName2, true);

				System.out.println("\t" + readName2 + "\t" + positiveMaps + "\t" + negativeMaps + "\t" + scorePos + "\t" + scoreNeg);
			}
		}
	}

	/**
	 * Read all sequences
	 */
	void readSeqs() {
		Timer.showStdErr("Reading seuqnece (FASTQ) file: " + fastq);

		FastqFileIterator fqfi = new FastqFileIterator(fastq);
		for (Fastq fq : fqfi) {
			String name = fq.getDescription().substring(1);
			// System.out.println("NAME:" + name);
			seqByName.put(name, fq.getSequence());
		}

		Timer.showStdErr("done. Sequences: " + seqByName.size());
	}
}
