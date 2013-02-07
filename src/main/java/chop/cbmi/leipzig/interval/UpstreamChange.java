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
            //distance to cds start
            //this may be much farther than the distance to the TSS reported by SnpEff
            //especially for alternate 5' end isoforms
            int distance = (transcript.isStrandPlus() ? transcript.getCdsStart() - seqChange.getEnd() : seqChange.getStart() - transcript.getCdsEnd());
            txPos = String.valueOf(-1*distance);
            change.setTxPos(txPos);
            return change;
        }
        //}
        return change;

    }
}
