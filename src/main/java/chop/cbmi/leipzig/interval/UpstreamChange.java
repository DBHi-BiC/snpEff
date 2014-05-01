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

    public UpstreamChange(SeqChange seqChange, Upstream upstream) {
        super(seqChange, (Transcript) upstream.getParent());
        this.upstream = upstream;
    }

    @Override
    ChangeEffect transcriptChange() {
        ChangeEffect change = changeEffect.clone();
        if (upstream.intersects(seqChange)) {
            //distance to cds start
            //this may be much farther than the distance to the TSS reported by SnpEff
            //especially for alternate 5' end isoforms
            Integer relativePosSt = (transcript.isStrandPlus() ? seqChange.getStart() - transcript.getCdsStart() : transcript.getCdsStart() - seqChange.getEnd());
            Integer relativePosEnd = (transcript.isStrandPlus() ? seqChange.getEnd() - transcript.getCdsStart() : transcript.getCdsStart() - seqChange.getStart());

            change = hgvsChangeFormatter(change, null, relativePosSt, relativePosEnd);

            //txPos = "-1"+String.valueOf(distanceToEndOfUTR);

            change.setTxPos(txPos);
            return change;
        }
        //}
        return change;

    }
}
