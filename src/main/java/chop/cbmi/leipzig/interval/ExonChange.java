package chop.cbmi.leipzig.interval;

import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;

/**
 * Jeremy Leipzig
 * Children's Hospital of Philadelphia
 * leipzig@gmail.com
 * 2/6/13
 */
public class ExonChange extends TranscriptChange {
    Exon exon;

    public ExonChange(SeqChange seqChange, Exon exon, ChangeEffect changeEffect) {
        super(seqChange, (Transcript) exon.getParent(), changeEffect);
        this.exon = exon;
    }

    @Override
    ChangeEffect transcriptChange() {
        ChangeEffect change = changeEffect.clone();
        try {
            txPos = String.valueOf(cdsBaseNumberOfExonInTx(seqChange.getStart()));
            change.setTxPos(txPos);
        } catch (IndexOutOfBoundsException e) {
            //sometimes a splice site will claim it belongs to an exon when it doesn't
             txPos = null;
        }

        return change;
    }
}
