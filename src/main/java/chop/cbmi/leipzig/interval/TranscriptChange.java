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
    int txPos = -1; //transcript-relative nt position
    String netCdsChange = "";

    public TranscriptChange(SeqChange seqChange, Transcript transcript, ChangeEffect changeEffect) {
        this.seqChange = seqChange;
        this.transcript = transcript;
        this.changeEffect = changeEffect;
    }


    // Get coding start (after 5 prime UTR)
    int cdsStart = transcript.getCdsStart();

    /**
     * Get some details about the effect on this transcript
     * @param seqChange
     * @return
     */
    @Override
    public seqChangeEffect(SeqChange seqChange, ChangeEffect changeEffect) {
        if (!intersects(seqChange)) return ChangeEffect.emptyResults(); // Sanity check

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



    /**
     * Calculate a list of transcript changes at the DNA level
     * and set txPos for use in HGVS-DNA
     * @return
     */
    List<ChangeEffect> transcriptChange() {
        ArrayList<ChangeEffect> changes = new ArrayList<ChangeEffect>();
        if (!transcript.intersects(seqChange)) return changes;

        List<Exon> exons = transcript.sortedStrand();
        List<Intron> introns = transcript.introns();
        List<Utr5prime> utr5primes = transcript.get5primeUtrs();
        List<Utr3prime> utr3primes = transcript.get3primeUtrs();


        //---
        // Hits a UTR region?
        //---
        boolean included = false;
        for (Utr utr : utrs)
            if (utr.intersects(seqChange)) {
                // Calculate the effect
                List<ChangeEffect> chEffList = utr.seqChangeEffect(seqChange, changeEffect.clone());
                if (!chEffList.isEmpty()) changes.addAll(chEffList);
                included |= utr.includes(seqChange); // Is this seqChange fully included in the UTR?
            }
        if (included) return changes; // SeqChange fully included in the UTR? => We are done.




        //3'UTR
        for (Utr3prime utr3prime : utr3primes){
            if(utr3prime.intersects(seqChange)){
                ChangeEffect changeEffectNew = changeEffect.clone();
                changeEffectNew.setMarker(utr3prime);
                changes.add(changeEffectNew);

                //this is just the txPos like any exon
                //except it is off the end of the transcript
                transcript.
                txPos = -(transcript.getFirstCodingExon().distance(seqChange));
                transcript.get
        }
            //5'UTR
            for (Utr5prime utr5prime : utr5primes){
                if(utr5prime.intersects(seqChange)){
                    ChangeEffect changeEffectNew = changeEffect.clone();
                    changeEffectNew.setMarker(utr5prime);
                    changes.add(changeEffectNew);

                    //should be reported as a negative number
                    txPos = -(transcript.getFirstCodingExon().distance(seqChange));
                }

        // Get coding start (after 5 prime UTR)
        int cdsStart = transcript.getCdsStart();

        // We may have to calculate 'netCdsChange', which is the effect on the CDS
        netCdsChange = netCdsChange();

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
                txPos = (firstCdsBaseInExon + cdsBaseInExon);

                // Use appropriate method to calculate codon change
                boolean hasChanged = false; // Was there any change?
                ChangeEffect changeEffectNew = changeEffect.clone();
                changeEffectNew.setMarker(exon); // It is affecting this exon, so we set the marker
                hasChanged = codonChangeSingle(changeEffectNew, exon);

                // Any change? => Add change to list
                if (hasChanged) {
                    codonsAround(seqChange, changeEffectNew, codonNum); // Show codons around change (if required)
                    changes.add(changeEffectNew);
                }

                // Can we return immediately?
                if (returnNow) return changes;
            }

            if (transcript.isStrandPlus()) firstCdsBaseInExon += Math.max(0, exon.getEnd() - Math.max(exon.getStart(), cdsStart) + 1);
            else firstCdsBaseInExon += Math.max(0, Math.min(cdsStart, exon.getEnd()) - exon.getStart() + 1);
        }

        return changes;
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
        if (seqChange.isIns()) return new TranscriptChangeIns(seqChange, transcript, changeEffect);
        if (seqChange.isDel()) return new TranscriptChangeDel(seqChange, transcript, changeEffect);
        if (seqChange.isMnp()) return new TranscriptChangeMnp(seqChange, transcript, changeEffect);
        if (seqChange.isInterval()) return new TranscriptChangeInterval(seqChange, transcript, changeEffect);
        throw new RuntimeException("Unimplemented factory for SeqChange: " + seqChange);
    }
}