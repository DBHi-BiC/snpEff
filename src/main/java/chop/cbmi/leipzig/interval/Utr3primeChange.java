package chop.cbmi.leipzig.interval;

import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;
import ca.mcgill.mcb.pcingola.interval.Utr3prime;

/**
 * Jeremy Leipzig
 * Children's Hospital of Philadelphia
 * leipzig@gmail.com
 * 2/6/13
 */
public class Utr3primeChange extends TranscriptChange {
    Utr3prime utr3prime;

    public Utr3primeChange(SeqChange seqChange, Utr3prime utr5prime, ChangeEffect changeEffect) {
        super(seqChange, (Transcript) utr5prime.getParent().getParent(), changeEffect);
        this.utr3prime = utr5prime;
    }

    @Override
    ChangeEffect transcriptChange() {
        ChangeEffect change = changeEffect.clone();

        //a star followed by distance to cds end
        //coding start
        if (utr3prime.intersects(seqChange)) {
            //distance to cds end
            //cdsStart>cdsEnd, so cdsStart is really cds start
            //this is the opposite of transcript start and end, in which tx.end > tx.start always
            Integer relativePosSt = (transcript.isStrandPlus() ? seqChange.getStart() - transcript.getCdsEnd() : transcript.getCdsEnd() - seqChange.getEnd());
            Integer relativePosEnd = (transcript.isStrandPlus() ? seqChange.getEnd() - transcript.getCdsEnd() : transcript.getCdsEnd() - seqChange.getStart());

            Exon exon = (Exon) utr3prime.findParent(Exon.class);

            change = hgvsChangeFormatter(change, exon, relativePosSt, relativePosEnd);

            //specific to 3' UTR
            txPos = "*"+txPos;
            change.setTxPos(txPos);
            return change;
        }
        return change;

    }
}
