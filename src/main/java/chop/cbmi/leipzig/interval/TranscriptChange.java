package chop.cbmi.leipzig.interval;

import java.util.ArrayList;
import java.util.List;

import ca.mcgill.mcb.pcingola.interval.*;
import ca.mcgill.mcb.pcingola.interval.codonChange.CodonChange;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.Config;

/**
 * Created with IntelliJ IDEA.
 * User: leipzig
 * Date: 2/1/13
 * Time: 2:06 PM
 * This class handles tx-level modifications and assigns appropriate positions
 * to allows associated HGVS descriptors to be generated
 * Unlike CodonChange.java it is also concerned with UTRs and Introns
 */
public class TranscriptChange {
    SeqChange seqChange;
    Transcript transcript;
    Exon exon = null;
    Intron intron = null;
    Utr5prime utr5prime = null;
    Utr3prime utr3prime = null;

    ChangeEffect changeEffect;
    String txPos = "-1"; //transcript-relative nt position
    String netCdsChange = "";
    private boolean usePrevBaseIntron = true;

    public TranscriptChange(SeqChange seqChange, Transcript transcript, ChangeEffect changeEffect) {
        this.seqChange = seqChange;
        this.transcript = transcript;
        this.changeEffect = changeEffect;
    }


   /* *//**
     * Get some details about the effect on this transcript
     * @param seqChange
     * @return
     *//*
    @Override
    public seqChangeEffect(SeqChange seqChange, ChangeEffect changeEffect) {


        // Create a list of changes
        ArrayList<ChangeEffect> changeEffectList = new ArrayList<ChangeEffect>();

        //---
        // Hits a UTR region?
        //---
        boolean included = false;
        for (Utr utr : utrs)
            if (utr.intersects(seqChange)) {
                // Calculate the effect
                List<ChangeEffect> chEffList = utr.seqChangeEffect(seqChange, changeEffect.clone());
                if (!chEffList.isEmpty()) changeEffectList.addAll(chEffList);
                included |= utr.includes(seqChange); // Is this seqChange fully included in the UTR?
            }
        if (included) return changeEffectList; // SeqChange fully included in the UTR? => We are done.


        // Does it hit an intron?
        for (Intron intron : introns())
            if (intron.intersects(seqChange)) {
                ChangeEffect cheff = changeEffect.clone();
                cheff.set(intron, ChangeEffect.EffectType.INTRON, "");
                changeEffectList.add(cheff);
            }

        //---
        // Analyze non-coding transcripts (or 'interval' seqChanges)
        //---
        if ((!Config.get().isTreatAllAsProteinCoding() && !isProteinCoding()) || seqChange.isInterval()) {
            // Do we have exon information for this transcript?
            if (!subintervals().isEmpty()) {
                // Add all exons
                for (Exon exon : this) {
                    if (exon.intersects(seqChange)) {
                        ChangeEffect cheff = changeEffect.clone();
                        cheff.set(exon, ChangeEffect.EffectType.EXON, "");
                        changeEffectList.add(cheff);
                    }
                }
            } else {
                // No exons annotated? Just mark it as hitting a transcript
                ChangeEffect cheff = changeEffect.clone();
                cheff.set(this, ChangeEffect.EffectType.TRANSCRIPT, "");
                changeEffectList.add(cheff);
            }

            return changeEffectList;
        }

        //---
        // This is a protein coding transcript.
        // We analyze codon replacement effect
        //---
        if (isCds(seqChange)) {
            // Get codon change effect
            CodonChange codonChange = new CodonChange(seqChange, this, changeEffect);
            changeEffectList.addAll(codonChange.calculate());
        }

        return changeEffectList;
    }
*/

    /**
     * Calculate the transcript change
     *
     * @param seqChange
     * @param changeEffect
     * @return
     */
    public ChangeEffect calculate() {
        ChangeEffect change = changeEffect.clone(); // Create a copy of this result
        //change.set(this.transcript, ChangeEffect.EffectType.HGVS,"HGVS");
        TranscriptChange transcriptChange = factory(seqChange, transcript, change);
        change = transcriptChange.transcriptChange();

        return change;
    }

    /**
     * Calculate base number in a CDS where 'pos' maps
     *
     * @returns Base number relative to cds
     */
     int cdsBaseNumberForAll(int pos) {

        // Calculate cdsStart and cdsEnd (if not already done)
        //calcCdsStartEnd();
         int cdsStart = transcript.getCdsStart();
        // All exons..
        int firstCdsBaseInExon = 0; // Where the exon maps to the CDS (i.e. which CDS base number does the first base in this exon maps to).
         List<Exon> exons = transcript.sortedStrand();
         for (Exon eint : exons) {
            if (eint.intersects(pos)) {
                int cdsBaseInExon; // cdsBaseInExon: base number relative to the beginning of the coding part of this exon (i.e. excluding 5'UTRs)
                if (transcript.getStrand() >= 0) cdsBaseInExon = pos - eint.getStart(); //Math.max(eint.getStart(), cdsStart);
                else cdsBaseInExon = eint.getEnd() - pos; //Math.min(eint.getEnd(), cdsStart) - pos;

                //cdsBaseInExon = Math.max(0, cdsBaseInExon);

                return firstCdsBaseInExon + cdsBaseInExon;
            } else {
                // Before exon begins?
                if ((transcript.isStrandPlus() && (pos < eint.getStart())) // Before exon begins (positive strand)?
                        || (transcript.isStrandMinus() && (pos > eint.getEnd()))) // Before exon begins (negative strand)?
                    return firstCdsBaseInExon - (this.usePrevBaseIntron ? 1 : 0);
            }

            if (transcript.isStrandPlus()){firstCdsBaseInExon += eint.getEnd()-eint.getStart()+1;
                //Math.max(0, eint.getEnd() - Math.max(eint.getStart(), cdsStart) + 1);
            }
            else {firstCdsBaseInExon += eint.getEnd() - eint.getStart() + 1;//Math.max(0, Math.min(cdsStart, eint.getEnd()) - eint.getStart() + 1);
            }
        }

        return firstCdsBaseInExon - 1;
    }

    /**
     * The transcript change at the DNA level
     * and set txPos for use in HGVS-DNA
     * @return
     */
    ChangeEffect transcriptChange() {
        ChangeEffect change = changeEffect.clone();
        if (!transcript.intersects(seqChange)) return change;

        List<Exon> exons = transcript.sortedStrand();
        List<Intron> introns = transcript.introns();
        List<Utr5prime> utr5primes = transcript.get5primeUtrs();
        List<Utr3prime> utr3primes = transcript.get3primeUtrs();

        //3'UTR
        for (Utr3prime utr3prime : utr3primes){
            if(utr3prime.intersects(seqChange)){
                //this is just the txPos like any exon
                //except it is off the end of the transcript
                txPos = String.valueOf(-(transcript.getFirstCodingExon().distance(seqChange)));
                change.setTxPos(txPos);
                return change;
            }
        }
        //5'UTR
        for (Utr5prime utr5prime : utr5primes){
            if(utr5prime.intersects(seqChange)){
                //should be reported as a negative number
                txPos = String.valueOf(-(transcript.getFirstCodingExon().distance(seqChange)));
                change.setTxPos(txPos);
                return change;
            }
        }

        // Get coding start (after 5 prime UTR)
        int cdsStart = transcript.getCdsStart();

        //intron
        //intronic nucleotides (coding DNA reference sequence only)
        //beginning of the intron; the number of the last nucleotide of the preceding exon, a plus sign and the position in the intron, like c.77+1G, c.77+2T, ....
        //end of the intron; the number of the first nucleotide of the following exon, a minus sign and the position upstream in the intron, like ..., c.78-2A, c.78-1G.
        //in the middle of the intron, numbering changes from "c.77+.." to "c.78-.."; for introns with an uneven number of nucleotides the central nucleotide is the last described with a "+"
        for (Intron intron : introns)
        {
            if (intron.intersects(seqChange)) {
                int firstAfter=transcript.firstExonPositionAfter(seqChange.getStart());
                int lastBefore=transcript.lastExonPositionBefore(seqChange.getStart());
                int cdsFirstAfter=cdsBaseNumberForAll(firstAfter);
                int cdsLastBefore=cdsBaseNumberForAll(lastBefore);
                int distanceToPrecedingExon=Math.abs(seqChange.getStart()-lastBefore);
                int distanceToProcedingExon=Math.abs(seqChange.getStart()-firstAfter);
                this.txPos = (distanceToPrecedingExon<distanceToProcedingExon) ? String.valueOf(cdsLastBefore)+"+"+String.valueOf(distanceToPrecedingExon) : String.valueOf(cdsFirstAfter)+"-"+String.valueOf(distanceToProcedingExon);
//                if(distanceToPrecedingExon<distanceToProcedingExon){
//                    this.txPos=String.valueOf(transcript.cdsBaseNumber(lastBefore))+"+"+String.valueOf(distanceToPrecedingExon);
//                }else{
//                    this.txPos=String.valueOf(transcript.cdsBaseNumber(firstAfter,true))+"-"+String.valueOf(distanceToProcedingExon);
//                }
                change.setTxPos(txPos);
                return change;
            }
        }


        //---
        // Concatenate all exons
        //---
        int firstCdsBaseInExon = 0; // Where the exon maps to the CDS (i.e. which CDS base number does the first base in this exon maps to).
        for (Exon exon : exons) {
            this.exon = exon;

            //hits an exon
            if (exon.intersects(seqChange)) {
                int cdsBaseInExon; // cdsBaseInExon: base number relative to the beginning of the coding part of this exon (i.e. excluding 5'UTRs)

                if (transcript.isStrandPlus()) {
                    int firstSeqChangeBaseInExon = Math.max(seqChange.getStart(), Math.max(exon.getStart(), cdsStart));
                    cdsBaseInExon = firstSeqChangeBaseInExon - Math.max(exon.getStart(), cdsStart);
                } else {
                    int lastSeqChangeBaseInExon = Math.min(seqChange.getEnd(), Math.min(exon.getEnd(), cdsStart));
                    cdsBaseInExon = Math.min(exon.getEnd(), cdsStart) - lastSeqChangeBaseInExon;
                }

                if (cdsBaseInExon < 0) cdsBaseInExon = 0;

                // Get codon number and index within codon (where seqChage is pointing)
                //codonNum = (firstCdsBaseInExon + cdsBaseInExon) / CODON_SIZE;
                //codonIndex = (firstCdsBaseInExon + cdsBaseInExon) % CODON_SIZE;

                //txPos is cdsBaseinTranscript
                txPos = String.valueOf((firstCdsBaseInExon + cdsBaseInExon));
            }

            if (transcript.isStrandPlus()) firstCdsBaseInExon += Math.max(0, exon.getEnd() - Math.max(exon.getStart(), cdsStart) + 1);
            else firstCdsBaseInExon += Math.max(0, Math.min(cdsStart, exon.getEnd()) - exon.getStart() + 1);
        }
        change.setTxPos(this.txPos);
        return change;
    }

    /**
     * Create a specific transcript change for a seqChange
     * @param seqChange
     * @param transcript
     * @param changeEffect
     * @return
     */
    TranscriptChange factory(SeqChange seqChange, Transcript transcript, ChangeEffect changeEffect) {
        if (seqChange.isSnp()) return new TranscriptChangeSnp(seqChange, transcript, changeEffect);
        /*if (seqChange.isIns()) return new TranscriptChangeIns(seqChange, transcript, changeEffect);
        if (seqChange.isDel()) return new TranscriptChangeDel(seqChange, transcript, changeEffect);
        if (seqChange.isMnp()) return new TranscriptChangeMnp(seqChange, transcript, changeEffect);
        if (seqChange.isInterval()) return new TranscriptChangeInterval(seqChange, transcript, changeEffect);  */
        throw new RuntimeException("Unimplemented factory for SeqChange: " + seqChange);
    }
}