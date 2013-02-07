package chop.cbmi.leipzig.interval;

import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;
import ca.mcgill.mcb.pcingola.interval.Utr5prime;

/**
 * Created with IntelliJ IDEA.
 * User: leipzigj
 * Date: 2/6/13
 * Time: 8:07 PM
 * To change this template use File | Settings | File Templates.
 */
public class Utr5primeChange extends ExonChange {
    public Utr5primeChange(SeqChange seqChange, Utr5prime utr5prime, ChangeEffect changeEffect) {
        super(seqChange, (Exon) utr5prime.getParent(), changeEffect);
    }

    @Override
    ChangeEffect transcriptChange() {
        ChangeEffect change = changeEffect.clone();

        //this should simply be a negative number representing the distance to the
        //coding star
        if (utr5prime.intersects(seqChange)) {
            txPos = String.valueOf(cdsBaseNumberForAll(seqChange.getStart()));
            change.setTxPos(txPos);
            return change;
        }
        //}
        return change;

    }
}
