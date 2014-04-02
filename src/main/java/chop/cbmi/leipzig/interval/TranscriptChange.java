package chop.cbmi.leipzig.interval;

import java.util.List;

import ca.mcgill.mcb.pcingola.interval.*;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.util.GprSeq;

/**
 * Jeremy Leipzig
 * Children's Hospital of Philadelphia
 * leipzig@gmail.com
 * 2/6/13
 *
 * This class handles tx-level modifications and assigns appropriate positions
 * to allows associated HGVS descriptors to be generated
 * Unlike CodonChange.java it is also concerned with UTRs and Introns
 */
public class TranscriptChange {
    SeqChange seqChange;
    Transcript transcript;
    ChangeEffect changeEffect;
    String txPos = null; //transcript-relative nt position
    final int HGVSOFFSET = 1;  //1-based
    private boolean usePrevBaseIntron = true;


    public TranscriptChange(SeqChange seqChange, Transcript transcript, ChangeEffect changeEffect) {
        this.seqChange = seqChange;
        this.transcript = transcript;
        this.changeEffect = changeEffect;
    }

    //common hgvs formatter for transcripts
    //txPos is for exons by default
    ChangeEffect hgvsChangeFormatter(ChangeEffect change, Exon exon, Integer relativePosSt, Integer relativePosEnd) {
        int dupOffset;
        Integer stPos;
        Integer endPos;
        try {
            if(seqChange.isDel()){
                //end>start because a strand is given
                if(transcript.isStrandPlus()){
                    dupOffset=repeatWalker(exon,change);
                    stPos=relativePosSt+dupOffset;
                    endPos=relativePosEnd+dupOffset;
                }else{
                    //you don't have to walk?
                    dupOffset=0;
                    stPos=relativePosEnd+dupOffset;
                    endPos=relativePosSt+dupOffset;
                }
                //sanity check - do the nucleotides specified in the deletion actually exist
                String deletedSequence;
                deletedSequence = transcript.cds().substring(stPos - 1, endPos).toUpperCase();
                //deletedSequence = exon.getSequence().substring(stPos-exon.getStart(),endPos-exon.getStart()).toUpperCase();
                if(!change.getNtDel().equals(deletedSequence))
                    throw new RuntimeException("the nucleotides at " + stPos + "-" + endPos + "(" + deletedSequence + ") are not the same as what ntDel as set to:" + change.getNtDel());

                if(seqChange.size()==1){
                    assert(stPos==endPos);
                    txPos= String.valueOf(stPos);
                }else{
                    if(stPos>endPos){
                        txPos=String.valueOf(endPos)+"_"+String.valueOf(stPos);
                    }else{
                        txPos= String.valueOf(stPos)+"_"+String.valueOf(endPos);
                    }
                }
            }else if(seqChange.isIns()){
                //we only use the startpos for insertions
                Integer hgvs_ins_offset;
                if(transcript.isStrandPlus()){
                    hgvs_ins_offset=-1;
                }else{
                    hgvs_ins_offset=0;
                }
                dupOffset=repeatWalker(exon,change);

                if(change.isDup()){
                    //when the dup loop ends you are change+offset
                    //the position is the first nucleotide of the proceeding sequence
                    //the same way a regular insert would,
                    //so the duplicated part precedes that nucleotide

                    //this logic seems one off for
                    //1	235967924	CI021900	A	AT	.	.	CLASS=DM;MUT=ALT;GENE=LYST;STRAND=-;DNA=NM_000081.2:c.3434dupA;DB=rs80338653;PHEN="Chediak-Higashi syndrome";
                    Integer ntLen=changeEffect.getNtIns().length();
                    if(ntLen==1){
                        //we need to break here
                        //61dupC
                        if(transcript.isStrandPlus()){
                            stPos = relativePosSt + dupOffset - 1;
                        }else{
                            stPos=relativePosSt+dupOffset;
                        }
                        txPos = String.valueOf(stPos);
                        change.setTxPos(txPos);
                        return change;
                    }
                    if(transcript.isStrandPlus()){
                        stPos=relativePosSt+dupOffset-ntLen;
                        endPos=relativePosSt+dupOffset-1;
                    }else{
                        //for negative strand dups the position is the end nt of the repeat
                        //and the last repeat
                        stPos=relativePosSt+dupOffset-ntLen+1;
                        endPos=relativePosSt+dupOffset;
                    }

                    //sanity check - do the nucleotides duplicated actually exist
                    String insertedSequence;
                    insertedSequence = transcript.cds().substring(stPos-1,endPos).toUpperCase();
                    if(!change.getNtIns().equals(insertedSequence))
                        throw new RuntimeException("the nucleotides at " + stPos + "-" + endPos + "(" + insertedSequence + ") are not the same as what ntIns as set to:" + change.getNtIns());
                }else{
                    stPos=relativePosSt+hgvs_ins_offset;
                    endPos=relativePosSt+hgvs_ins_offset+1;
                }


                if(stPos>endPos){
                    txPos=String.valueOf(endPos)+"_"+String.valueOf(stPos);
                }else{
                    txPos= String.valueOf(stPos)+"_"+String.valueOf(endPos);
                }



                seqChange.setStart(seqChange.getStart()+dupOffset);
            }else{
                //a SNP, start is fine
                txPos= String.valueOf(relativePosSt);
            }
            change.setTxPos(txPos);
        } catch (IndexOutOfBoundsException e) {
            //sometimes a splice site will claim it belongs to an exon when it doesn't
            txPos = null;
        }

        return change;
    }
    /**
     * Calculate the transcript change
     *
     * @param seqChange
     * @param changeEffect
     * @return
     */
    public ChangeEffect calculate() {
        ChangeEffect change = changeEffect.clone(); // Create a copy of this result
        if(seqChange.isSnp()) setSNP();
        if(seqChange.isIns()) setINS();
        if(seqChange.isDel()) setDEL();
        change = this.transcriptChange();

        return change;
    }


    /**
     * Calculate base number of an exon position
     *
     */
    int cdsBaseNumberOfExonInTx(int pos){
        List<Exon> exons = transcript.sortedStrand();
        for (Exon eint : exons) {
            if (eint.intersects(pos)) {
                int cdsBaseNumber = this.transcript.cdsBaseNumber(pos,usePrevBaseIntron);
                return cdsBaseNumber+HGVSOFFSET;
            }
        }
        //sometimes this happens with complex deletions
        throw new IndexOutOfBoundsException(pos + " is not in an exon but "+this.getClass()+" asked for its cds base number");
    }

    /**
     * The transcript change at the DNA level
     * and set txPos for use in HGVS-DNA
     * @return
     */
    ChangeEffect transcriptChange() {
        ChangeEffect change = changeEffect.clone();
        if (!transcript.intersects(seqChange)) return change;

        this.txPos= String.valueOf(cdsBaseNumberOfExonInTx(seqChange.getStart()));

        change.setTxPos(this.txPos);
        return change;
    }


    public void setSNP(){
        if(transcript.isStrandPlus()){
            changeEffect.setTranscript(this.txPos, seqChange.reference(), seqChange.change(), null, null);
        }else{
            changeEffect.setTranscript(this.txPos, GprSeq.reverseWc(seqChange.reference()), GprSeq.reverseWc(seqChange.change()), null, null);
        }
    }
    public void setINS(){
        if(transcript.isStrandPlus()){
            changeEffect.setTranscript(this.txPos, seqChange.reference(), null, seqChange.netChange(seqChange), null);
        }else{
            changeEffect.setTranscript(this.txPos, GprSeq.reverseWc(seqChange.reference()), null, GprSeq.reverseWc(seqChange.netChange(seqChange)), null);
        }
    }
    public void setDEL(){
        if(transcript.isStrandPlus()){
            changeEffect.setTranscript(this.txPos, seqChange.reference(), null, null, seqChange.netChange(seqChange));
        }else{
            changeEffect.setTranscript(this.txPos, GprSeq.reverseWc(seqChange.reference()), null, null, GprSeq.reverseWc(seqChange.netChange(seqChange)));
        }
    }

    //walk an indel looking for optimal position for hgvs
    public int repeatWalker(Exon exon,ChangeEffect change){
        //make sure there is enough flank left to check,an insert of length 3 must be at position 4 or later
        //ntLen is the size of the indel
        //dupOffset is the adjustment for placing ins 3' of repeats if they are dups
        Integer changeBaseInExon;
        Integer ntLen;
        Integer dupOffset=0;
        Integer rollOffset=0;
        //we need to know whether to use insertion or deletion string
        CharStack flank;
        if(exon == null){
            return 0;
        }
        if(seqChange.isDel()){

            changeBaseInExon = seqChange.getStart() - exon.getStart();
            ntLen=changeEffect.getNtDel().length();

            flank=new CharStack(changeEffect.getNtDel());

            if (transcript.isStrandPlus()){
                //this must have been a walk-off the end warning but it isn't written correctly
                //if(changeBaseInExon-ntLen>=0){
                    boolean walking = true;
                    boolean rolling = false;
                    while(walking || rolling){
                        String postFlank="";
                        //walk the duplication

                        //String testFlank=exon.getSequence().substring(2361,2370);
                        int sPos=changeBaseInExon+dupOffset-rollOffset;
                        int ePos=changeBaseInExon+ntLen+dupOffset-rollOffset;
                        try{
                            postFlank=exon.getSequence().substring(sPos,ePos).toUpperCase();
                        }catch(StringIndexOutOfBoundsException e){
                            //you have no room to check
                            //this occurs with
                            //1       76199309        CD052077        TG      T       .       .       CLASS=DM;MUT=ALT;GENE=ACADM;STRAND=+;DNA=NM_000016.4:c.387+1delG
                            walking=false;
                            if(change.isDup()){
                                rolling=true;
                                rollOffset+=1;
                            }
                            continue;
                        }
                        if(postFlank.equals(flank.get())){
                            change.setDup(true);
                            if(walking){
                                //try to walk further, might fail
                                dupOffset=dupOffset+ntLen;
                            }
                            if(rolling){
                                //retract ntlen later if a failure
                                //not sure why I added +ntLen;
                                dupOffset=dupOffset-rollOffset;

                                //rolling was a success
                                if(seqChange.isDel()){
                                    change.setNtDel(flank.get());
                                }else{
                                    change.setNtIns(flank.get());
                                }
                                rolling=false;
                                //can we walk again?
                                //walking=true;
                                rollOffset=0;
                            }
                        }
                        else{
                            walking=false;
                            if(change.isDup()){
                                //you don't get to roll unless you have walked at least one dup
                                rolling=true;
                                //here we need to rollback because these deletions will produce the same sequence but we want the latter
                                //torollbackrollbacsome
                                //  rollback
                                //         krollbac
                                rollOffset+=1;
                                //are you back where you started?
                                if(rollOffset==ntLen){
                                    rolling=false;
                                    //rolling was a failure
                                    dupOffset=dupOffset-ntLen;
                                }
                                else{
                                    flank.rollback();
                                }
                            }
                        }
                    //}
                }
            }
        }else{
            if(seqChange.isIns()){
                //this is used almost nowhere else but since we get exon sequence instead of tx seq we use it
                if (transcript.isStrandPlus()) changeBaseInExon = seqChange.getStart() - exon.getStart();
                else changeBaseInExon = exon.getEnd() - seqChange.getStart();
                ntLen=changeEffect.getNtIns().length();

                flank=new CharStack(changeEffect.getNtIns());
                if (transcript.isStrandPlus()){
                    //if(changeBaseInExon-ntLen>=0){
                        boolean walking = true;
                        boolean rolling = false;
                        boolean continue_flag=true;
                        while(continue_flag & (walking || rolling)){
                            String postFlank="";
                            String preFlank="";
                            //walk the duplication

                            //String testFlank=exon.getSequence().substring(2361,2370);
                            int sPos=changeBaseInExon+dupOffset-rollOffset;
                            int ePos=changeBaseInExon+ntLen+dupOffset-rollOffset;
                            //remember substring is exclusive
                            //ePos is the position after the insert
                            try{
                                postFlank=exon.getSequence().substring(sPos,ePos).toUpperCase();
                                preFlank=exon.getSequence().substring(sPos-ntLen,sPos).toUpperCase();
                            }catch(StringIndexOutOfBoundsException e){
                                //you have no room to check
                                continue_flag=false;
                                break;
                            }
                            if(walking & postFlank.equals(flank.get())){
                                change.setDup(true);
                                change.setNtIns(flank.get());
                                //try to walk further, might fail
                                dupOffset=dupOffset+ntLen;
                            }else{
                                if(!change.isDup()){
                                    continue_flag=false;
                                }
                                if(rolling){
                                    //if the frame is correct
                                    //dupOffset=dupOffset-rollOffset+ntLen;
                                    //rolling was a success
                                    if(postFlank.equals(flank.get())){
                                        change.setNtIns(flank.get());
                                        rolling=false;
                                        //the insertion is at the end of teh postflank
                                        dupOffset=dupOffset-rollOffset+ntLen;
                                        rollOffset=0;
                                        //TODO: not elegant
                                        continue_flag=false;
                                    }
                                }
                                if(walking){
                                    walking=false;
                                    rolling=true;
                                    //for a insertion of 3
                                    //we know to try 2 away then 1 away
                                    //there are only 3 permutations

                                    //I dont' get this
                                    //dupOffset=dupOffset+ntLen-1;
                                }
                                if(continue_flag){
                                    rollOffset+=1;
                                    //are you back where you started?
                                    if(rollOffset==ntLen){
                                        rolling=false;
                                        continue_flag=false;
                                        //rolling was a failure
                                        //dupOffset=dupOffset-ntLen;
                                    }
                                    else{
                                        flank.rollback();
                                        //dupOffset=dupOffset-1;
                                    }
                                }
                            }

                        }
                    //}
                }else{
                    //for negative strand inserts we need to look behind
                    try{
                        //String preFlank=exon.getSequence().substring(changeBaseInExon-2,changeBaseInExon+ntLen-2).toUpperCase();
                        String preFlank=exon.getSequence().substring(changeBaseInExon-ntLen-2,changeBaseInExon-2).toUpperCase();

                        if(preFlank.equals(flank.get()) & seqChange.isIns()){
                            change.setDup(true);
                        }
                    }catch(StringIndexOutOfBoundsException e){
                        //no room to check
                        //4	39455843	.	T	TA	100.0	PASS	DP=100
                        return 0;
                    }

                }
            }else{
                //not a del or insert
                return 0;
            }
        }
        return dupOffset;
    }

    public String intronFormat(int position){
        int firstAfter;
        int lastBefore;
        String beforeString;
        String afterString;
        //cdsStart can be > cdsEnd, so cdsStart is really cds start

        //check to see if this intron is in the utr
        //the rule says find the nearest utr nucleotide
        //http://www.hgvs.org/mutnomen/refseq.html#IVSin5

        lastBefore = transcript.isStrandPlus() ? transcript.lastExonPositionBefore(position) : transcript.firstExonPositionAfter(position);
        firstAfter = transcript.isStrandPlus() ? transcript.firstExonPositionAfter(position) : transcript.lastExonPositionBefore(position);

        if(transcript.isStrandPlus() & position>transcript.getCdsEnd() || !transcript.isStrandPlus() &position<transcript.getCdsEnd()){
            //3'utr
            int lastBefore3Utr = (transcript.isStrandPlus() ? lastBefore - transcript.getCdsEnd() : transcript.getCdsEnd() - lastBefore);
            int firstAfter3Utr = (transcript.isStrandPlus() ? firstAfter - transcript.getCdsEnd() : transcript.getCdsEnd() - firstAfter);
            beforeString="*"+String.valueOf(lastBefore3Utr);
            afterString="*"+String.valueOf(firstAfter3Utr);
        }else{
            if(transcript.isStrandPlus() & position<transcript.getCdsStart() || !transcript.isStrandPlus() & position>transcript.getCdsStart()){
                //5utr
                int lastBefore5Utr = (transcript.isStrandPlus() ? lastBefore - transcript.getCdsStart() : transcript.getCdsStart() - lastBefore);
                int firstAfter5Utr = (transcript.isStrandPlus() ? firstAfter - transcript.getCdsStart() : transcript.getCdsStart() - firstAfter);
                beforeString=String.valueOf(lastBefore5Utr);
                afterString=String.valueOf(firstAfter5Utr);
            }else{
                //normal intron
                int cdsFirstAfter= cdsBaseNumberOfExonInTx(firstAfter);
                int cdsLastBefore= cdsBaseNumberOfExonInTx(lastBefore);
                beforeString=String.valueOf(cdsLastBefore);
                afterString=String.valueOf(cdsFirstAfter);
            }
        }


        int toProceeding=Math.abs(position-firstAfter);
        int fromPreceeding=Math.abs(position-lastBefore);

        String fromPreceedingString = ((fromPreceeding == 0)   ? "" : "+"+String.valueOf(fromPreceeding));
        String toProceedingString = ((toProceeding == 0)   ? "" : "-"+String.valueOf(toProceeding));

        String intronFormat = (fromPreceeding<toProceeding) ? beforeString + fromPreceedingString : afterString + toProceedingString;
        return intronFormat;
    }

    //am I closer to the preceeding exon than the proceeding one?
    public Boolean comingFrom(int position){
        int firstAfter = transcript.isStrandPlus() ? transcript.firstExonPositionAfter(position) : transcript.lastExonPositionBefore(position);
        int lastBefore = transcript.isStrandPlus() ? transcript.lastExonPositionBefore(position) : transcript.firstExonPositionAfter(position);
        int toProceeding=Math.abs(position-firstAfter);
        int fromPreceeding=Math.abs(position-lastBefore);
        return fromPreceeding<toProceeding;
    }

    public void allIntronTxPos(){
        String txPosStringSt=intronFormat(seqChange.getStart());
        String txPosStringEnd=intronFormat(seqChange.getEnd());

        if(seqChange.isDel()){
            if(seqChange.size()==1){
                txPos= txPosStringSt;
            }else{
                if(transcript.isStrandPlus()){
                    txPos=txPosStringSt+"_"+txPosStringEnd;
                }else{
                    txPos=txPosStringEnd+"_"+txPosStringSt;
                }
            }
        }else if(seqChange.isIns()){
            Integer hgvs_ins_offset;
            if(transcript.isStrandPlus()){
                hgvs_ins_offset=-1;
            }else{
                hgvs_ins_offset=0;
            }

            int stPos=seqChange.getStart()+hgvs_ins_offset;
            int endPos=seqChange.getStart()+hgvs_ins_offset+1;
            String insSt=intronFormat(stPos);
            String endSt=intronFormat(endPos);

            //666-3-666-2
            //if(comingFrom(stPos)){
            txPos = transcript.isStrandPlus() ? insSt+"_"+endSt : endSt+"_"+insSt;
//            }else{
//                txPos= transcript.isStrandPlus() ? insSt+"_"+endSt : endSt+"_"+insSt;
//            }
            //txPos=txPosStringSt+intronFormat(seqChange.getStart()+1);
        }else{
            txPos=intronFormat(seqChange.getStart());
        }
    }
}