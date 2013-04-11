package chop.cbmi.leipzig.interval;

import java.util.List;

import ca.mcgill.mcb.pcingola.interval.*;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
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
    final int HGVSOFFSET = 1;  //1-based
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
        if(seqChange.isIns()) setINS();
        if(seqChange.isDel()) setDEL();
        change = this.transcriptChange();

        return change;
    }


    /**
     * Calculate base number of an exon position
     * 
     */
    int cdsBaseNumberOfExonInTx(int pos){
        List<Exon> exons = transcript.sortedStrand();
        for (Exon eint : exons) {
            if (eint.intersects(pos)) {
                int cdsBaseNumber = this.transcript.cdsBaseNumber(pos,usePrevBaseIntron);
                return cdsBaseNumber+HGVSOFFSET;
            }
        }
        throw new IndexOutOfBoundsException(pos + " is not in an exon but "+this.getClass()+" asked for its cds base number");
    }

    /**
     * The transcript change at the DNA level
     * and set txPos for use in HGVS-DNA
     * @return
     */
    ChangeEffect transcriptChange() {
        ChangeEffect change = changeEffect.clone();
        if (!transcript.intersects(seqChange)) return change;

        this.txPos= String.valueOf(cdsBaseNumberOfExonInTx(seqChange.getStart()));

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
    public void setINS(){
        if(transcript.isStrandPlus()){
            changeEffect.setTranscript(this.txPos, seqChange.reference(), null, seqChange.netChange(seqChange), null);
        }else{
            changeEffect.setTranscript(this.txPos, GprSeq.reverseWc(seqChange.reference()), null, GprSeq.reverseWc(seqChange.netChange(seqChange)), null);
        }
    }
    public void setDEL(){
        if(transcript.isStrandPlus()){
            changeEffect.setTranscript(this.txPos, seqChange.reference(), null, null, seqChange.netChange(seqChange));
        }else{
            changeEffect.setTranscript(this.txPos, GprSeq.reverseWc(seqChange.reference()), null, null, GprSeq.reverseWc(seqChange.netChange(seqChange)));
        }
    }
}