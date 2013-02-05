package chop.cbmi.leipzig.interval;

import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.ErrorType;
/**
 * Set SNP details for DNA-HGVS
 */
public class TranscriptChangeSnp extends TranscriptChange{
    public TranscriptChangeSnp(SeqChange seqChange, Transcript transcript, ChangeEffect changeEffect) {
        super(seqChange, transcript, changeEffect);
        //pos,oldnt,newnt,insnt,delnt
        this.txPos = getTxPos();
        changeEffect.setTranscript(this.txPos, seqChange.reference(), seqChange.change(), null, null);
    }
}