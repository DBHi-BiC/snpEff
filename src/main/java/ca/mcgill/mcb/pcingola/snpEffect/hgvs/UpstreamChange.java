package ca.mcgill.mcb.pcingola.snpEffect.hgvs;

import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.interval.Upstream;
import ca.mcgill.mcb.pcingola.interval.Variant;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;


/**
 * Jeremy Leipzig
 * Children's Hospital of Philadelphia
 * leipzig@gmail.com
 * 5/1/14
 */
public class UpstreamChange extends TranscriptChange {
    Upstream upstream;

    public UpstreamChange(Variant seqChange, Upstream upstream, ChangeEffect changeEffect) {
        super(seqChange, (Transcript) upstream.getParent(), changeEffect);
        this.upstream = upstream;
    }

    @Override
    boolean transcriptChange() {
        if (upstream.intersects(seqChange)) {
            //distance to cds start
            //this may be much farther than the distance to the TSS reported by SnpEff
            //especially for alternate 5' end isoforms
            Integer relativePosSt = (transcript.isStrandPlus() ? seqChange.getStart() - transcript.getCdsStart() : transcript.getCdsStart() - seqChange.getEnd());
            Integer relativePosEnd = (transcript.isStrandPlus() ? seqChange.getEnd() - transcript.getCdsStart() : transcript.getCdsStart() - seqChange.getStart());

            hgvsChangeFormatter(changeEffect, null, relativePosSt, relativePosEnd);

            //txPos = "-1"+String.valueOf(distanceToEndOfUTR);

            changeEffect.setTxPos(txPos);
            return true;
        }
        //}
        return false;

    }
}
