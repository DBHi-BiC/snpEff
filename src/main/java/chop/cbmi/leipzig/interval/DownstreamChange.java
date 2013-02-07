package chop.cbmi.leipzig.interval;

import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.interval.Downstream;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;

/**
 * Jeremy Leipzig
 * Children's Hospital of Philadelphia
 * leipzig@gmail.com
 * 2/6/13
 *
 * for downstream changes
 * the nucleotide 3' of the translation stop codon is *1, the next *2, etc.
 */
public class DownstreamChange  extends TranscriptChange {
    Downstream downstream;

    public DownstreamChange(SeqChange seqChange, Downstream downstream, ChangeEffect changeEffect) {
        super(seqChange, (Transcript) downstream.getParent(), changeEffect);
        this.downstream = downstream;
    }



    @Override
    ChangeEffect transcriptChange() {
        ChangeEffect change = changeEffect.clone();


        if (downstream.intersects(seqChange)) {
            int distance = (transcript.isStrandPlus() ? seqChange.getStart() - downstream.getStart() : downstream.getEnd() - seqChange.getStart()) + 1;
            txPos = "*"+String.valueOf(distance);
            change.setTxPos(txPos);
            return change;
        }
        //}
        return change;

    }
}
