package chop.cbmi.leipzig.interval;

import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;
import ca.mcgill.mcb.pcingola.interval.Utr5prime;

/**
 * Jeremy Leipzig
 * Children's Hospital of Philadelphia
 * leipzig@gmail.com
 * 2/6/13
 */
public class Utr5primeChange extends TranscriptChange {
    Utr5prime utr5prime;

    public Utr5primeChange(SeqChange seqChange, Utr5prime utr5prime, ChangeEffect changeEffect) {
        super(seqChange, (Transcript) utr5prime.getParent().getParent(), changeEffect);
        this.utr5prime = utr5prime;
    }

    @Override
    ChangeEffect transcriptChange() {
        ChangeEffect change = changeEffect.clone();

        //this should simply be a negative number representing the distance to the
        //coding start
        if (utr5prime.intersects(seqChange)) {
            int distance = (transcript.isStrandPlus() ? transcript.getCdsStart() - seqChange.getEnd() : seqChange.getStart() - transcript.getCdsEnd());

            txPos = String.valueOf(-1*distance);
            change.setTxPos(txPos);
            return change;

        }
        //}
        return change;

    }
}
