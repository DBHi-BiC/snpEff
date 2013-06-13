package chop.cbmi.leipzig.interval;

import ca.mcgill.mcb.pcingola.interval.Intron;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;

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
        //beginning of the intron; the number of the last nucleotide of the preceeding exon, a plus sign and the position in the intron, like c.77+1G, c.77+2T, ....
        //end of the intron; the number of the first nucleotide of the following exon, a minus sign and the position upstream in the intron, like ..., c.78-2A, c.78-1G.
        //in the middle of the intron, numbering changes from "c.77+.." to "c.78-.."; for introns with an uneven number of nucleotides the central nucleotide is the last described with a "+"
        //for (Intron intron : introns)
        //{
            if (intron.intersects(seqChange)) {
                change.set(intron, EffectType.INTRON, "");
                int firstAfter = transcript.isStrandPlus() ? transcript.firstExonPositionAfter(seqChange.getEnd()) : transcript.lastExonPositionBefore(seqChange.getEnd());
                int lastBefore = transcript.isStrandPlus() ? transcript.lastExonPositionBefore(seqChange.getStart()) : transcript.firstExonPositionAfter(seqChange.getStart());
                int cdsFirstAfter= cdsBaseNumberOfExonInTx(firstAfter);
                int cdsLastBefore= cdsBaseNumberOfExonInTx(lastBefore);
                int distanceToPreceedingExonStart=Math.abs(seqChange.getStart()-lastBefore);
                int distanceToProceedingExonStart=Math.abs(seqChange.getStart()-firstAfter);
                int distanceToPreceedingExonEnd=Math.abs(seqChange.getEnd()-lastBefore);
                int distanceToProceedingExonEnd=Math.abs(seqChange.getEnd()-firstAfter);
                //666-3-666-2
                //bigger is from
                int fromProceeding = (distanceToProceedingExonStart <= distanceToProceedingExonEnd) ? distanceToProceedingExonEnd : distanceToProceedingExonStart;
                int toProceeding =  (distanceToProceedingExonStart <= distanceToProceedingExonEnd) ? distanceToProceedingExonStart : distanceToProceedingExonEnd;

                //665+1-665+2
                int fromPreceeding = (distanceToPreceedingExonStart <= distanceToPreceedingExonEnd) ? distanceToPreceedingExonStart :distanceToPreceedingExonEnd;
                int toPreceeding =  (distanceToPreceedingExonStart <= distanceToPreceedingExonEnd) ? distanceToPreceedingExonEnd : distanceToPreceedingExonStart;

                String fromPreceedingString = ((fromPreceeding == 0) ? "" : "+"+String.valueOf(fromPreceeding));
                String fromProceedingString = ((fromProceeding == 0) ? "" : "-"+String.valueOf(fromProceeding));
                String toPreceedingString =   ((toPreceeding == 0)   ? "" : "+"+String.valueOf(toPreceeding));
                String toProceedingString =   ((toProceeding == 0)   ? "" : "-"+String.valueOf(toProceeding));

                if(seqChange.isDel()){
                    String txPosStart = (distanceToPreceedingExonStart<distanceToProceedingExonStart) ? String.valueOf(cdsLastBefore)+ fromPreceedingString : String.valueOf(cdsFirstAfter)+ fromProceedingString;
                    String txPosEnd = (distanceToPreceedingExonEnd<distanceToProceedingExonEnd) ? String.valueOf(cdsLastBefore)+toPreceedingString : String.valueOf(cdsFirstAfter)+toProceedingString;

                    if(seqChange.size()==1){
                        txPos= txPosStart;
                    }else{
                        //if(transcript.isStrandPlus()){
                             txPos=txPosStart+"_"+txPosEnd;
                        //}else{
                        //    txPos=txPosEnd+"_"+txPosStart;
                        //}
                    }
                }else if(seqChange.isIns()){
                    txPos = (distanceToPreceedingExonStart<distanceToProceedingExonStart) ? String.valueOf(cdsLastBefore)+"+"+String.valueOf(distanceToPreceedingExonStart+1)+"_"+String.valueOf(cdsLastBefore)+"+"+String.valueOf(distanceToPreceedingExonStart+2) : String.valueOf(cdsFirstAfter)+"-"+String.valueOf(distanceToProceedingExonStart)+"_"+String.valueOf(cdsFirstAfter)+"-"+String.valueOf(distanceToProceedingExonStart+1);
                }else{
                    txPos = (distanceToPreceedingExonStart<distanceToProceedingExonStart) ? String.valueOf(cdsLastBefore)+"+"+String.valueOf(distanceToPreceedingExonStart) : String.valueOf(cdsFirstAfter)+"-"+String.valueOf(distanceToProceedingExonStart);
                }


                change.setTxPos(txPos);
                return change;
            }
        //}
        return change;

    }
}
