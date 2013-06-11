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
        //exon-specific relative positions
        Integer relativePosSt;
        Integer relativePosEnd;
        try{
            relativePosSt=cdsBaseNumberOfExonInTx(seqChange.getStart());
        }catch (IndexOutOfBoundsException e){
            int lastBefore = transcript.isStrandPlus() ? transcript.lastExonPositionBefore(seqChange.getStart()) : transcript.firstExonPositionAfter(seqChange.getStart());
            relativePosSt= cdsBaseNumberOfExonInTx(lastBefore);
        }
        try{
            relativePosEnd=cdsBaseNumberOfExonInTx(seqChange.getEnd());
        }catch (IndexOutOfBoundsException e){
            int firstAfter = transcript.isStrandPlus() ? transcript.firstExonPositionAfter(seqChange.getEnd()) : transcript.lastExonPositionBefore(seqChange.getEnd());
            relativePosEnd= cdsBaseNumberOfExonInTx(firstAfter);
        }
        change = hgvsChangeFormatter(change, exon, relativePosSt, relativePosEnd);
        return change;
    }
}
//http://stackoverflow.com/a/8746524/264696
class CharStack {
    StringBuilder sb;// = new StringBuilder();

    public CharStack(String inputString){
        this.sb = new StringBuilder(inputString);
    }
    public void push(char ch) {
        sb.append(ch);
    }
    public void prepend(char ch) {
        sb.insert(0,ch);
    }
    public char pop() {
        int last = sb.length() -1;
        char ch= sb.charAt(last);
        sb.setLength(last);
        return ch;
    }

    public int size() {
        return sb.length();
    }

    public void rollback() {
        char popped=pop();
        prepend(popped);
    }

    public String get(){
        return String.valueOf(sb);
    }
}