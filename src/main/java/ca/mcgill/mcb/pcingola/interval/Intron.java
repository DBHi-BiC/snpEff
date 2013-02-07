package ca.mcgill.mcb.pcingola.interval;

import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;
import chop.cbmi.leipzig.interval.IntronChange;

import java.util.List;

/**
 * Intron
 * 
 * @author pcingola
 *
 */
public class Intron extends Marker {

	private static final long serialVersionUID = -8283322526157264389L;

	int rank; // Exon rank in transcript

	public Intron(Transcript parent, int start, int end, int strand, String id) {
		super(parent, start, end, strand, id);
		type = EffectType.INTRON;
	}

	public int getRank() {
		return rank;
	}

	public void setRank(int rank) {
		this.rank = rank;
	}

    /**
     * Calculate the effect of this seqChange
     * @param seqChange
     * @param changeEffect
     * @return
     */
    @Override
    public List<ChangeEffect> seqChangeEffect(SeqChange seqChange, ChangeEffect changeEffect) {
        if (!intersects(seqChange)) return ChangeEffect.emptyResults(); // Sanity check

        IntronChange intronChange = new IntronChange(seqChange, this, changeEffect);
        changeEffect = intronChange.calculate();
        changeEffect.set(this, type, "");
        return changeEffect.newList();
    }
}
