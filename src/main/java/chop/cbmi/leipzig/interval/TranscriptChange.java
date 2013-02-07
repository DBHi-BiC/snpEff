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
    String txPos = null; //transcript-relative nt position
    int hgvsoffset = 1;  //1-based
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
     * Calculate base number in a CDS where 'pos' maps
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
                if (transcript.getStrand() >= 0) cdsBaseInExon = pos - eint.getStart();
                else cdsBaseInExon = eint.getEnd() - pos;
                return firstCdsBaseInExon + cdsBaseInExon + hgvsoffset;
            } else {
                // Before exon begins?
                if ((transcript.isStrandPlus() && (pos < eint.getStart())) // Before exon begins (positive strand)?
                        || (transcript.isStrandMinus() && (pos > eint.getEnd()))) // Before exon begins (negative strand)?
                {
                    //System.out.println(pos + " is not in an exon but "+this.getClass()+" asked for its cds base number");
                    throw new IndexOutOfBoundsException(pos + " is not in an exon but "+this.getClass()+" asked for its cds base number");
                    //return firstCdsBaseInExon - (usePrevBaseIntron ? 1 : 0);
                }
            }

            if (transcript.isStrandPlus()) firstCdsBaseInExon += eint.getEnd()-eint.getStart()+1;
            else firstCdsBaseInExon += eint.getEnd() - eint.getStart() + 1;

        }
        return firstCdsBaseInExon - 1 + hgvsoffset;
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