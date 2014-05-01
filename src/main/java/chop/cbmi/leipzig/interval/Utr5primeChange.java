package chop.cbmi.leipzig.interval;

import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
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
    boolean transcriptChange() {

        //this should simply be a negative number representing the distance to the
        //coding start
        if (utr5prime.intersects(seqChange)) {
            //this is the method employed by SnpEff
//            List<Utr5prime> utrs = transcript.get5primeUtrs();
//            boolean fromEnd = !(transcript.getStrand() < 0); // We want distance from beginning of transcript (TSS = End of 5'UTR)
//            int distanceToEndOfUTR = seqChange.distanceFrom(utrs, fromEnd) + 1;

            //distance to cds start
            //cdsStart>cdsEnd, so cdsStart is really cds start
            //this is the opposite of transcript start and end, in which tx.end > tx.start always
            //we are going for negative distances here since that fits in with the HGVS nomenclature
            Integer relativePosSt = (transcript.isStrandPlus() ? seqChange.getStart() - transcript.getCdsStart() : transcript.getCdsStart() - seqChange.getEnd());
            Integer relativePosEnd = (transcript.isStrandPlus() ? seqChange.getEnd() - transcript.getCdsStart() : transcript.getCdsStart() - seqChange.getStart());

            Exon exon = (Exon) utr5prime.findParent(Exon.class);

            changeEffect = hgvsChangeFormatter(changeEffect, exon, relativePosSt, relativePosEnd);

            //txPos = "-1"+String.valueOf(distanceToEndOfUTR);

            changeEffect.setTxPos(txPos);
            return true;

        }
        //}
        return false;

    }
}
