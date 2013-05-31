package chop.cbmi.leipzig.interval;

import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;
import ca.mcgill.mcb.pcingola.interval.Utr5prime;

import java.util.List;

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
            int distanceToCodingStart = (transcript.isStrandPlus() ? transcript.getCdsStart() - seqChange.getEnd() : seqChange.getStart() - transcript.getCdsStart());

            List<Utr5prime> utrs = transcript.get5primeUtrs();
            boolean fromEnd = !(transcript.getStrand() < 0); // We want distance from begining of transcript (TSS = End of 5'UTR)
            int distanceToEndOfUTR = seqChange.distanceFrom(utrs, fromEnd) + 1;

            txPos = String.valueOf(-1*distanceToEndOfUTR);
            change.setTxPos(txPos);
            return change;

        }
        //}
        return change;

    }
}
