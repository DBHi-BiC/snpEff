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

        //a star followed by distance to cds end
        if (downstream.intersects(seqChange)) {
            //distance to cds end
            //this may be much farther than the distance to the polyA reported by SnpEff
            //especially for alternate 3' end isoforms
            int distance = (transcript.isStrandPlus() ? seqChange.getStart() - transcript.getCdsEnd() : transcript.getCdsStart() - seqChange.getEnd());
            txPos = "*"+String.valueOf(distance);

            change.setTxPos(txPos);
            return change;
        }
        //}
        return change;

    }
}
