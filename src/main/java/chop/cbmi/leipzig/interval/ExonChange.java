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
        Integer relativePosSt=cdsBaseNumberOfExonInTx(seqChange.getStart());
        Integer relativePosEnd=cdsBaseNumberOfExonInTx(seqChange.getEnd());
        change= hgvsChangeFormatter(change, relativePosSt, relativePosEnd);
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