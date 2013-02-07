package chop.cbmi.leipzig.interval;

import ca.mcgill.mcb.pcingola.interval.Intron;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;

import java.util.List;

/**
 * Jeremy Leipzig
 * Children's Hospital of Philadelphia
 * leipzig@gmail.com
 * 2/6/13
 */
public class IntronChange extends TranscriptChange {
    Intron intron;

    public IntronChange(SeqChange seqChange, Intron intron, ChangeEffect changeEffect) {
        super(seqChange, (Transcript) intron.getParent(), changeEffect);
        this.intron = intron;
    }

    @Override
    ChangeEffect transcriptChange() {
        ChangeEffect change = changeEffect.clone();
        //List<Intron> introns = transcript.introns();

        //intron
        //intronic nucleotides (coding DNA reference sequence only)
        //beginning of the intron; the number of the last nucleotide of the preceding exon, a plus sign and the position in the intron, like c.77+1G, c.77+2T, ....
        //end of the intron; the number of the first nucleotide of the following exon, a minus sign and the position upstream in the intron, like ..., c.78-2A, c.78-1G.
        //in the middle of the intron, numbering changes from "c.77+.." to "c.78-.."; for introns with an uneven number of nucleotides the central nucleotide is the last described with a "+"
        //for (Intron intron : introns)
        //{
            if (intron.intersects(seqChange)) {
                change.set(intron, EffectType.INTRON, "");
                int firstAfter = transcript.isStrandPlus() ? transcript.firstExonPositionAfter(seqChange.getEnd()) : transcript.lastExonPositionBefore(seqChange.getEnd());
                int lastBefore = transcript.isStrandPlus() ? transcript.lastExonPositionBefore(seqChange.getStart()) : transcript.firstExonPositionAfter(seqChange.getStart());
                int cdsFirstAfter=cdsBaseNumberForAll(firstAfter);
                int cdsLastBefore=cdsBaseNumberForAll(lastBefore);
                int distanceToPrecedingExon=Math.abs(seqChange.getStart()-lastBefore);
                int distanceToProcedingExon=Math.abs(seqChange.getStart()-firstAfter);
                this.txPos = (distanceToPrecedingExon<distanceToProcedingExon) ? String.valueOf(cdsLastBefore)+"+"+String.valueOf(distanceToPrecedingExon) : String.valueOf(cdsFirstAfter)+"-"+String.valueOf(distanceToProcedingExon);
                change.setTxPos(txPos);
                return change;
            }
        //}
        return change;

    }
}
