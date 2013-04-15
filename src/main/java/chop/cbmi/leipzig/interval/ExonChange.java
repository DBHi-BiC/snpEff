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
            if(seqChange.isDel()){
                if(seqChange.size()==1){
                    txPos= String.valueOf(cdsBaseNumberOfExonInTx(seqChange.getStart()));
                }else{
                    if(transcript.isStrandPlus()){
                        txPos= String.valueOf(cdsBaseNumberOfExonInTx(seqChange.getStart()))+"_"+String.valueOf(cdsBaseNumberOfExonInTx(seqChange.getEnd()));
                    }else{
                        txPos= String.valueOf(cdsBaseNumberOfExonInTx(seqChange.getEnd()))+"_"+String.valueOf(cdsBaseNumberOfExonInTx(seqChange.getStart()));
                    }
                }
            }else if(seqChange.isIns()){
                txPos= String.valueOf(cdsBaseNumberOfExonInTx(seqChange.getStart()))+"_"+String.valueOf(cdsBaseNumberOfExonInTx(seqChange.getStart())+1);
                //we need to get the flanking sequence somehow
                //seems only exons get blessed with this info
                String myCodingSequence=this.exon.getSequence();

                //this is used almost nowhere else but since we get exon sequence instead of tx seq we use it
                Integer changeBaseInExon;
                if (transcript.isStrandPlus()) changeBaseInExon = seqChange.getStart() - this.exon.getStart();
                else changeBaseInExon = this.exon.getEnd() - seqChange.getStart();

                //make sure there is enough flank left to check,an insert of length 3 must be at position 4 or later
                Integer ntLen=changeEffect.getNtIns().length();
                if(changeBaseInExon-ntLen>=0){
                    String preFlank=exon.getSequence().substring(changeBaseInExon-ntLen-1,changeBaseInExon-1).toUpperCase();
                    if(preFlank.equals(changeEffect.getNtIns())){
                        seqChange.setDup(true);
                    }
                }
            }else{
                txPos= String.valueOf(cdsBaseNumberOfExonInTx(seqChange.getStart()));
            }


            change.setTxPos(txPos);
        } catch (IndexOutOfBoundsException e) {
            //sometimes a splice site will claim it belongs to an exon when it doesn't
             txPos = null;
        }

        return change;
    }
}
