package chop.cbmi.leipzig.interval;

import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.util.GprSeq;

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

    private int repeatWalker(int changeBaseInExon,ChangeEffect change,int ntLen){
        //make sure there is enough flank left to check,an insert of length 3 must be at position 4 or later
        //ntLen is the size of the insert
        //dupOffset is the adjustment for placing ins 3' of repeats if they are dups
        Integer dupOffset=0;
        //we need to know whether to use insertion or deletion string
        String flank="";
        if(seqChange.isDel()){
            flank=changeEffect.getNtDel();
        }
        if(seqChange.isIns()){
            flank=changeEffect.getNtIns();
        }
        if (transcript.isStrandPlus()){
            if(changeBaseInExon-ntLen>=0){
                boolean still_walking = true;
                while(still_walking){
                    //walk the duplication
                        String postFlank=exon.getSequence().substring(changeBaseInExon+dupOffset,changeBaseInExon+ntLen+dupOffset).toUpperCase();
                        if(postFlank.equals(flank) & seqChange.isIns()){
                            change.setDup(true);
                            dupOffset=dupOffset+ntLen;
                        }else{
                            still_walking=false;
                        }

                }
            }
        }else{
           //for negative strand inserts we need to look behind
            String preFlank=exon.getSequence().substring(changeBaseInExon-2,changeBaseInExon+ntLen-2).toUpperCase();
            if(preFlank.equals(flank) & seqChange.isIns()){
                change.setDup(true);
            }
        }
        return dupOffset;
    }
    @Override
    ChangeEffect transcriptChange() {
        ChangeEffect change = changeEffect.clone();
        try {
            if(seqChange.isDel()){
                if(seqChange.size()==1){
                    txPos= String.valueOf(cdsBaseNumberOfExonInTx(seqChange.getStart()));
                }else{
                    //end>start because a strand is given
                    if(transcript.isStrandPlus()){
                        Integer changeBaseInExon;
                        changeBaseInExon = seqChange.getEnd() - this.exon.getStart();
                        Integer ntLen=changeEffect.getNtDel().length();
                        int dupOffset=repeatWalker(changeBaseInExon,change,ntLen);
                        txPos= String.valueOf(cdsBaseNumberOfExonInTx(seqChange.getStart())+dupOffset)+"_"+String.valueOf(cdsBaseNumberOfExonInTx(seqChange.getEnd())+dupOffset);
                    }else{
                        //you don't have to walk?
                        txPos= String.valueOf(cdsBaseNumberOfExonInTx(seqChange.getEnd()))+"_"+String.valueOf(cdsBaseNumberOfExonInTx(seqChange.getStart()));
                    }
                }
            }else if(seqChange.isIns()){
                Integer hgvs_ins_offset;
                //TODO:i do not understand this at all
                if(transcript.isStrandPlus()){
                    hgvs_ins_offset=-1;
                }else{
                    hgvs_ins_offset=0;
                }

                //we need to get the flanking sequence somehow
                //seems only exons get blessed with this info
                String myCodingSequence=this.exon.getSequence();

                //this is used almost nowhere else but since we get exon sequence instead of tx seq we use it
                Integer changeBaseInExon;
                if (transcript.isStrandPlus()) changeBaseInExon = seqChange.getStart() - this.exon.getStart();
                else changeBaseInExon = this.exon.getEnd() - seqChange.getStart();
                Integer ntLen=changeEffect.getNtIns().length();
                int dupOffset=repeatWalker(changeBaseInExon,change,ntLen);

                if(change.isDup()){
                    //when the dup loop ends you are change+offset is sitting on the last base of the repeat
                    txPos= String.valueOf(cdsBaseNumberOfExonInTx(seqChange.getStart())+hgvs_ins_offset+dupOffset-ntLen+1)+"_"+String.valueOf(cdsBaseNumberOfExonInTx(seqChange.getStart())+hgvs_ins_offset+dupOffset);
                }else{
                    txPos= String.valueOf(cdsBaseNumberOfExonInTx(seqChange.getStart())+hgvs_ins_offset)+"_"+String.valueOf(cdsBaseNumberOfExonInTx(seqChange.getStart())+1+hgvs_ins_offset);
                }
                seqChange.setStart(seqChange.getStart()+dupOffset);
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
