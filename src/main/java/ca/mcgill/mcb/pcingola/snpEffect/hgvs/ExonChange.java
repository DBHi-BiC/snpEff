package ca.mcgill.mcb.pcingola.snpEffect.hgvs;

import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.interval.Variant;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;

/**
 * Jeremy Leipzig
 * Children's Hospital of Philadelphia
 * leipzig@gmail.com
 * 5/1/14
 */
public class ExonChange extends TranscriptChange {
    Exon exon;

    public ExonChange(Variant seqChange, Exon exon, ChangeEffect changeEffect) {
        super(seqChange, (Transcript) exon.getParent(), changeEffect);
        this.exon = exon;
    }

    @Override
    boolean transcriptChange() {
        //ChangeEffect change = changeEffect.clone();

        //exon-specific relative positions
        Integer relativePosSt;
        Integer relativePosEnd;
        try {
            relativePosSt = cdsBaseNumberOfExonInTx(seqChange.getStart());
        } catch (IndexOutOfBoundsException e) {
            int lastBefore = transcript.isStrandPlus() ? transcript.lastExonPositionBefore(seqChange.getStart()) : transcript.firstExonPositionAfter(seqChange.getStart());
            if (lastBefore == -1) {

                //off the edge of the cds, but is it the start or the end?
                if (transcript.isStrandPlus()) {
                    //utr5
                    relativePosSt = seqChange.getStart() - transcript.getCdsStart();

                    //wrap up and send off
                    relativePosEnd = cdsBaseNumberOfExonInTx(seqChange.getEnd());
                    txPos = String.valueOf(relativePosSt) + "_" + String.valueOf(relativePosEnd);
                    changeEffect.setTxPos(txPos);
                    return true;

                } else {
                    //utr3
                    relativePosSt = transcript.getCdsEnd() - seqChange.getStart();

                    relativePosEnd = cdsBaseNumberOfExonInTx(seqChange.getEnd());
                    txPos = String.valueOf(relativePosEnd) + "_" + String.valueOf(relativePosSt);
                    changeEffect.setTxPos(txPos);
                    return true;
                }

            } else {
                //intronic - assume were closer to this guy
                String relativePosEndString;
                try {
                    relativePosEnd = cdsBaseNumberOfExonInTx(seqChange.getEnd());
                    relativePosEndString = String.valueOf(relativePosEnd);
                    if (transcript.isStrandPlus()) {
                        txPos = intronFormat(seqChange.getStart()) + "_" + relativePosEndString;
                    } else {
                        txPos = relativePosEndString + "_" + intronFormat(seqChange.getStart());
                    }
                } catch (IndexOutOfBoundsException f) {
                    //totally intronic
                    allIntronTxPos();
                }
                changeEffect.setTxPos(txPos);
                return true;
            }
        }
        try {
            relativePosEnd = cdsBaseNumberOfExonInTx(seqChange.getEnd());
        } catch (IndexOutOfBoundsException e) {
            int firstAfter = transcript.isStrandPlus() ? transcript.firstExonPositionAfter(seqChange.getEnd()) : transcript.lastExonPositionBefore(seqChange.getEnd());

            if (firstAfter == -1) {

                //off the edge of the gene, but is it the start or the end?
                if (transcript.isStrandPlus()) {
                    //utr3
                    relativePosEnd = seqChange.getEnd() - transcript.getCdsEnd();
                    txPos = String.valueOf(relativePosSt) + "_*" + String.valueOf(relativePosEnd);
                    changeEffect.setTxPos(txPos);
                    return true;
                } else {
                    //utr5 - this implies neg strand
                    //the change is strand-ignorant - it's the getEnd that is utr5
                    relativePosEnd = transcript.getCdsStart() - seqChange.getEnd();
                    txPos = String.valueOf(relativePosEnd) + "_" + String.valueOf(relativePosSt);
                    changeEffect.setTxPos(txPos);
                    return true;
                }

            } else {
                //intronic - assume were closer to this guy
                if (transcript.isStrandPlus()) {
                    txPos = intronFormat(seqChange.getEnd()) + "_" + String.valueOf(relativePosSt);
                } else {
                    txPos = String.valueOf(relativePosSt) + "_" + intronFormat(seqChange.getEnd());
                }
                changeEffect.setTxPos(txPos);
                return true;
            }
        }
        hgvsChangeFormatter(changeEffect, exon, relativePosSt, relativePosEnd);
        return true;
    }


}

//http://stackoverflow.com/a/8746524/264696
class CharStack {
    StringBuilder sb;// = new StringBuilder();

    public CharStack(String inputString) {
        this.sb = new StringBuilder(inputString);
    }

    public void push(char ch) {
        sb.append(ch);
    }

    public void prepend(char ch) {
        sb.insert(0, ch);
    }

    public char pop() {
        int last = sb.length() - 1;
        char ch = sb.charAt(last);
        sb.setLength(last);
        return ch;
    }

    public int size() {
        return sb.length();
    }

    public void rollback() {
        char popped = pop();
        prepend(popped);
    }

    public String get() {
        return String.valueOf(sb);
    }
}

