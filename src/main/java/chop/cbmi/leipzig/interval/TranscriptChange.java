package chop.cbmi.leipzig.interval;

import java.util.ArrayList;
import java.util.List;

import ca.mcgill.mcb.pcingola.interval.*;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;

/**
 * This class handles tx-level modifications and assigns appropriate positions aka txPos
 * to allow associated HGVS descriptors to be generated
 * Unlike CodonChange.java it is also concerned with UTRs and Introns
 */
public class TranscriptChange {
    SeqChange seqChange;
    Transcript transcript;
    //Exon exon = null;
    //Intron intron = null;
    //Utr5prime utr5prime = null;
    //Utr3prime utr3prime = null;

    ChangeEffect changeEffect;
    int txPos = -1; //transcript-relative nt position
    //String netCdsChange = "";

    public TranscriptChange(SeqChange seqChange, Transcript transcript, ChangeEffect changeEffect) {
        this.seqChange = seqChange;
        this.transcript = transcript;
        this.changeEffect = changeEffect;
    }


    // Get coding start (after 5 prime UTR)
    int cdsStart = transcript.getCdsStart();

    /**
     * Calculate the transcript change at the DNA level
     * and set txPos for use in HGVS-DNA
     * @return
     */
    ChangeEffect transcriptChange() {
        ChangeEffect change = new ChangeEffect(this.seqChange);

        if (!transcript.intersects(seqChange)) return change;

        List<Exon> exons = transcript.sortedStrand();
        List<Intron> introns = transcript.introns();
        List<Utr5prime> utr5primes = transcript.get5primeUtrs();
        List<Utr3prime> utr3primes = transcript.get3primeUtrs();

        //this is all about getting txPos
        //3'UTR
        for (Utr3prime utr3prime : utr3primes){
            if(utr3prime.intersects(seqChange)){
                change.setMarker(utr3prime);
                //this is just the txPos like any exon
                //except it is off the end of the transcript
                this.txPos = -(transcript.getFirstCodingExon().distance(seqChange));
                break;
            }
        }
        //5'UTR
        for (Utr5prime utr5prime : utr5primes){
            if(utr5prime.intersects(seqChange)){
                //should be reported as a negative number
                this.txPos = -(transcript.getFirstCodingExon().distance(seqChange));
                break;
            }
        }

        // Get coding start (after 5 prime UTR)
        int cdsStart = transcript.getCdsStart();
        //---
        // Concatenate all exons
        //---
        int firstCdsBaseInExon = 0; // Where the exon maps to the CDS (i.e. which CDS base number does the first base in this exon maps to).
        for (Exon exon : exons) {
            //this.exon = exon;

            //hits an exon
            if (exon.intersects(seqChange)) {
                int cdsBaseInExon; // cdsBaseInExon: base number relative to the beginning of the coding part of this exon (i.e. excluding 5'UTRs)

                if (transcript.isStrandPlus()) {
                    int firstSeqChangeBaseInExon = Math.max(seqChange.getStart(), Math.max(exon.getStart(), cdsStart));
                    cdsBaseInExon = firstSeqChangeBaseInExon - Math.max(exon.getStart(), cdsStart);
                } else {
                    int lastSeqChangeBaseInExon = Math.min(seqChange.getEnd(), Math.min(exon.getEnd(), cdsStart));
                    cdsBaseInExon = Math.min(exon.getEnd(), cdsStart) - lastSeqChangeBaseInExon;
                }

                if (cdsBaseInExon < 0) cdsBaseInExon = 0;

                //txPos is cdsBaseinTranscript
                txPos = (firstCdsBaseInExon + cdsBaseInExon);
            }
            if (transcript.isStrandPlus()) firstCdsBaseInExon += Math.max(0, exon.getEnd() - Math.max(exon.getStart(), cdsStart) + 1);
            else firstCdsBaseInExon += Math.max(0, Math.min(cdsStart, exon.getEnd()) - exon.getStart() + 1);
        }

        return changes;
    }

    /**
     * Create a specific transcript change for a seqChange based on snp/indel/etc
     * @param seqChange
     * @param transcript
     * @param changeEffect
     * @return
     */
    TranscriptChange factory(SeqChange seqChange, Transcript transcript, ChangeEffect changeEffect) {
        if (seqChange.isSnp()) return new TranscriptChangeSnp(seqChange, transcript, changeEffect);
        if (seqChange.isIns()) return new TranscriptChangeIns(seqChange, transcript, changeEffect);
        if (seqChange.isDel()) return new TranscriptChangeDel(seqChange, transcript, changeEffect);
        if (seqChange.isMnp()) return new TranscriptChangeMnp(seqChange, transcript, changeEffect);
        if (seqChange.isInterval()) return new TranscriptChangeInterval(seqChange, transcript, changeEffect);
        throw new RuntimeException("Unimplemented factory for SeqChange: " + seqChange);
    }
}