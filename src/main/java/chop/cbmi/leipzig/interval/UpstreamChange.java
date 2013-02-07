package chop.cbmi.leipzig.interval;

import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.interval.Upstream;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;


/**
 * Jeremy Leipzig
 * Children's Hospital of Philadelphia
 * leipzig@gmail.com
 * 2/6/13
 */
public class UpstreamChange extends TranscriptChange {
    Upstream upstream;

    public UpstreamChange(SeqChange seqChange, Upstream upstream, ChangeEffect changeEffect) {
        super(seqChange, (Transcript) upstream.getParent(), changeEffect);
        this.upstream = upstream;
    }



    @Override
    ChangeEffect transcriptChange() {
        ChangeEffect change = changeEffect.clone();


        if (upstream.intersects(seqChange)) {
            int distance = (transcript.isStrandPlus() ? upstream.getEnd() - seqChange.getStart() : seqChange.getStart() - upstream.getStart()) + 1;

            txPos = String.valueOf(-1*distance);
            change.setTxPos(txPos);
            return change;
        }
        //}
        return change;

    }
}
