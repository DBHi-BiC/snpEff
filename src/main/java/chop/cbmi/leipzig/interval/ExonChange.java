package chop.cbmi.leipzig.interval;

import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;

import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: leipzigj
 * Date: 2/6/13
 * Time: 8:07 PM
 * To change this template use File | Settings | File Templates.
 */
public class ExonChange extends TranscriptChange {
    public ExonChange(SeqChange seqChange, Exon exon, ChangeEffect changeEffect) {
        super(seqChange, (Transcript) exon.getParent(), changeEffect);
    }

    @Override
    ChangeEffect transcriptChange() {
        ChangeEffect change = changeEffect.clone();

        if (exon.intersects(seqChange)) {
            txPos = String.valueOf(cdsBaseNumberForAll(seqChange.getStart()));
            change.set(exon, EffectType.EXON, "");
            change.setTxPos(txPos);
            return change;
        }
        //}
        return change;

    }
}
