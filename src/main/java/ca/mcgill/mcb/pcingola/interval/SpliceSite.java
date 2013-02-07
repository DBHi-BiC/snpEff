package ca.mcgill.mcb.pcingola.interval;

import java.util.List;

import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import chop.cbmi.leipzig.interval.ExonChange;
import chop.cbmi.leipzig.interval.TranscriptChange;

/**
 * Interval for a splice site
 * 
 * Reference: http://en.wikipedia.org/wiki/RNA_splicing
 * 
 * Spliceosomal introns often reside in eukaryotic protein-coding genes. Within the 
 * intron, a 3' splice site, 5' splice site, and branch site are required for splicing. 
 * The 5' splice site or splice donor site includes an almost invariant sequence GU 
 * at the 5' end of the intron, within a larger, less highly conserved consensus region. 
 * The 3' splice site or splice acceptor site terminates the intron with an almost 
 * invariant AG sequence. Upstream (5'-ward) from the AG there is a region high in 
 * pyrimidines (C and U), or polypyrimidine tract. Upstream from the polypyrimidine 
 * tract is the branch point, which includes an adenine nucleotide.
 * 
 * @author pcingola
 *
 */
public abstract class SpliceSite extends Marker {

	private static final long serialVersionUID = 1636197649250882952L;
	public static final int CORE_SPLICE_SITE_SIZE = 2;

	public SpliceSite() {
		super();
	}

	public SpliceSite(Exon parent, int start, int end, int strand, String id) {
		super(parent, start, end, strand, id);
	}

	public SpliceSite(Transcript parent, int start, int end, int strand, String id) {
		super(parent, start, end, strand, id);
	}

	/**
	 * Core splice sites are defined as CORE_SPLICE_SITE_SIZE bases after exon end 
	 * or before exon begins. Usually CORE_SPLICE_SITE_SIZE is 2 bases.
	 * Other spice sites are considered "non-core".
	 * 
	 * @return
	 */
	public abstract boolean intersectsCoreSpliceSite(Marker marker);

	/**
	 * Splice sites are no included in Exons, by definition.
	 */
	@Override
	protected boolean isShowWarningIfParentDoesNotInclude() {
		return false;
	}

	@Override
	public List<ChangeEffect> seqChangeEffect(SeqChange seqChange, ChangeEffect changeEffect) {
		if (!intersects(seqChange)) return ChangeEffect.emptyResults(); // Sanity check


        if(parent instanceof Exon){
            System.out.println("exon"+seqChange.getStart());
            ExonChange ExonChange = new ExonChange(seqChange, (Exon) parent, changeEffect);
            changeEffect = ExonChange.calculate();
        }else if(parent instanceof Transcript){
            System.out.println("tx"+seqChange.getStart());
            Transcript myTxParent = (Transcript) this.parent;
            TranscriptChange transcriptChange = new TranscriptChange(seqChange, myTxParent, changeEffect);
            changeEffect = transcriptChange.calculate();
        }

		changeEffect.set(this, type, "");
		return changeEffect.newList();
	}
}
