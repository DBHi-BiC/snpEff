package chop.cbmi.leipzig.interval;

import java.util.ArrayList;
import java.util.List;

import ca.mcgill.mcb.pcingola.interval.*;
import ca.mcgill.mcb.pcingola.interval.codonChange.CodonChange;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.util.GprSeq;

/**
 * Jeremy Leipzig
 * Children's Hospital of Philadelphia
 * leipzig@gmail.com
 * 2/6/13
 *
 * This class handles tx-level modifications and assigns appropriate positions
 * to allows associated HGVS descriptors to be generated
 * Unlike CodonChange.java it is also concerned with UTRs and Introns
 */
public class TranscriptChange {
    SeqChange seqChange;
    Transcript transcript;
    ChangeEffect changeEffect;
    String txPos = "-1"; //transcript-relative nt position
    private boolean usePrevBaseIntron = true;

    public TranscriptChange(SeqChange seqChange, Transcript transcript, ChangeEffect changeEffect) {
        this.seqChange = seqChange;
        this.transcript = transcript;
        this.changeEffect = changeEffect;
    }

    /**
     * Calculate the transcript change
     *
     * @param seqChange
     * @param changeEffect
     * @return
     */
    public ChangeEffect calculate() {
        ChangeEffect change = changeEffect.clone(); // Create a copy of this result
        if(seqChange.isSnp()) setSNP();
        change = this.transcriptChange();

        return change;
    }

    /**
     * Calculate base number in a CDS where 'pos' maps, allowing negative numbers
     *
     * @returns Base number relative to cds
     */
    int cdsBaseNumberForAll(int pos) {

        // Calculate cdsStart and cdsEnd (if not already done)
        int cdsStart = transcript.getCdsStart();
        int firstCdsBaseInExon = 0; // Where the exon maps to the CDS (i.e. which CDS base number does the first base in this exon maps to).
        List<Exon> exons = transcript.sortedStrand();
        for (Exon eint : exons) {
            if (eint.intersects(pos)) {
                int cdsBaseInExon; // cdsBaseInExon: base number relative to the beginning of the coding part of this exon (i.e. excluding 5'UTRs)
                if (transcript.getStrand() >= 0) cdsBaseInExon = pos - eint.getStart(); //Math.max(eint.getStart(), cdsStart);
                else cdsBaseInExon = eint.getEnd() - pos; //Math.min(eint.getEnd(), cdsStart) - pos;
                return firstCdsBaseInExon + cdsBaseInExon;
            } else {
                // Before exon begins?
                if ((transcript.isStrandPlus() && (pos < eint.getStart())) // Before exon begins (positive strand)?
                        || (transcript.isStrandMinus() && (pos > eint.getEnd()))) // Before exon begins (negative strand)?
                    return firstCdsBaseInExon - (this.usePrevBaseIntron ? 1 : 0);
            }

            if (transcript.isStrandPlus()) firstCdsBaseInExon += eint.getEnd()-eint.getStart()+1;
            else firstCdsBaseInExon += eint.getEnd() - eint.getStart() + 1;

        }

        return firstCdsBaseInExon - 1;
    }

    /**
     * The transcript change at the DNA level
     * and set txPos for use in HGVS-DNA
     * @return
     */
    ChangeEffect transcriptChange() {
        ChangeEffect change = changeEffect.clone();
        if (!transcript.intersects(seqChange)) return change;

        this.txPos= String.valueOf(cdsBaseNumberForAll(seqChange.getStart()));
        change.setTxPos(this.txPos);
        return change;
    }


    public void setSNP(){
        if(transcript.isStrandPlus()){
            changeEffect.setTranscript(this.txPos, seqChange.reference(), seqChange.change(), null, null);
        }else{
            changeEffect.setTranscript(this.txPos, GprSeq.reverseWc(seqChange.reference()), GprSeq.reverseWc(seqChange.change()), null, null);
        }
    }
}