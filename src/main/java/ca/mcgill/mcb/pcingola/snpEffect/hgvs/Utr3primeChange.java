package ca.mcgill.mcb.pcingola.snpEffect.hgvs;

import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.interval.Variant;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.interval.Utr3prime;

/**
 * Jeremy Leipzig
 * Children's Hospital of Philadelphia
 * leipzig@gmail.com
 * 5/1/14
 */
public class Utr3primeChange extends TranscriptChange {
    Utr3prime utr3prime;

    public Utr3primeChange(Variant seqChange, Utr3prime utr3prime, ChangeEffect changeEffect) {
        super(seqChange, (Transcript) utr3prime.getParent().getParent(), changeEffect);
        this.utr3prime = utr3prime;
    }

    @Override
    boolean transcriptChange() {
        ChangeEffect change = new ChangeEffect(seqChange);

        //a star followed by distance to cds end
        //coding start
        if (utr3prime.intersects(seqChange)) {
            Integer relativePosSt = (transcript.isStrandPlus() ? seqChange.getStart() - transcript.getCdsEnd() : transcript.getCdsEnd() - seqChange.getEnd());
            Integer relativePosEnd = (transcript.isStrandPlus() ? seqChange.getEnd() - transcript.getCdsEnd() : transcript.getCdsEnd() - seqChange.getStart());

            Exon exon = (Exon) utr3prime.findParent(Exon.class);

            //edge case
            //1       196716441       CD032470        TAGAA   T    is    NM_000186.3:c.3695_*2delAGAA
            if (relativePosSt < 1) {
                relativePosSt = cdsBaseNumberOfExonInTx(seqChange.getStart());
                hgvsChangeFormatter(change, exon, relativePosSt, relativePosEnd);
                txPos = txPos.replaceAll("_(\\d+)", "_\\*$1");
                change.setTxPos(txPos);
                return true;
            }


            hgvsChangeFormatter(change, exon, relativePosSt, relativePosEnd);


            //specific to 3' UTR
            //txPos = "*"+txPos;
            txPos = txPos.replaceAll("(\\d+)", "\\*$1");
            change.setTxPos(txPos);
            return true;
        }
        return false;

    }
}
