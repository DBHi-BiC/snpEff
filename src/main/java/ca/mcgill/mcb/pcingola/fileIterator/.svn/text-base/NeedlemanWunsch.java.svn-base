package ca.mcgill.mcb.pcingola.fileIterator;

/**
 * Needleman-Wunsch algorithm for string alignment (short strings, since it's not memory optimized)
 * 
 * @author pcingola
 */
public class NeedlemanWunsch {

	String alignment;
	char a[], b[];
	int score[][];
	int match = 1; // Match score
	int missMatch = -1; // Missmatch score
	int deletion = -2; // Deletion score
	int offset = 0;

	public NeedlemanWunsch(String a, String b) {
		score = new int[a.length() + 1][b.length() + 1];
		this.a = a.toCharArray();
		this.b = b.toCharArray();
	}

	public String align() {
		scoreMatrix();
		calcAlignments();
		return alignment;
	}

	void calcAlignments() {
		int maxLen = Math.max(a.length, b.length);
		char alignmentA[] = new char[maxLen];
		char alignmentB[] = new char[maxLen];

		int i = a.length;
		int j = b.length;
		int h = maxLen - 1;

		while((i > 0) && (j > 0) && (h >= 0)) {
			int s = getScore(i, j);
			int scorediag = getScore(i - 1, j - 1);
			int scoreup = getScore(i, j - 1);
			int scoreleft = getScore(i - 1, j);

			if( s == scoreup + deletion ) {
				alignmentA[h] = '-';
				alignmentB[h] = b[j - 1];
				j--;
			} else if( s == scoreleft + deletion ) {
				alignmentA[h] = a[i - 1];
				alignmentB[h] = '-';
				i--;
			} else if( s == scorediag + simmilarity(i, j) ) {
				alignmentA[h] = ' '; // a[i - 1];
				alignmentB[h] = ' '; // b[j - 1];
				i--;
				j--;
			} else throw new RuntimeException("This should never happen!");

			h--;
		}

		while((i > 0) && (h >= 0)) {
			alignmentA[h] = a[i - 1];
			alignmentB[h] = '-';
			i--;
			h--;
		}

		while((j > 0) && (h >= 0)) {
			alignmentA[h] = '-';
			alignmentB[h] = b[j - 1];
			j--;
			h--;
		}

		// Calculate offset from original position
		for( offset = 0; (offset < maxLen) && (alignmentA[offset] == ' '); offset++ );

		// Create alignment string
		StringBuffer alsb = new StringBuffer();
		char prev = ' ';
		for( i = 0; i < maxLen; i++ ) {
			if( alignmentA[i] == '-' ) {
				if( prev != '-' ) alsb.append('-');
				alsb.append(alignmentB[i]);
				prev = '-';
			} else if( alignmentB[i] == '-' ) {
				if( prev != '+' ) alsb.append('+');
				alsb.append(alignmentA[i]);
				prev = '+';
			}
		}
		alignment = alsb.toString();

	}

	public int getAligmentScore() {
		return getScore(a.length, b.length);
	}

	public String getAlignment() {
		return alignment;
	}

	public int getOffset() {
		return offset;
	}

	int getScore(int i, int j) {
		return score[i][j];
	}

	/**
	 * Calculate score matrix
	 */
	void scoreMatrix() {
		// Initialize
		for( int i = 0; i <= a.length; i++ )
			setScore(i, 0, deletion * i);

		for( int j = 0; j <= b.length; j++ )
			setScore(0, j, deletion * j);

		// Calculate
		for( int i = 1; i <= a.length; i++ )
			for( int j = 1; j <= b.length; j++ ) {
				int match = getScore(i - 1, j - 1) + simmilarity(i, j);
				int del = getScore(i - 1, j) + deletion;
				int ins = getScore(i, j - 1) + deletion;
				int s = Math.max(match, Math.max(del, ins));
				setScore(i, j, s);
			}
	}

	public void setDeletion(int deletion) {
		this.deletion = deletion;
	}

	public void setMatch(int match) {
		this.match = match;
	}

	public void setMissMatch(int missMatch) {
		this.missMatch = missMatch;
	}

	void setScore(int i, int j, int val) {
		score[i][j] = val;
	}

	/**
	 * Simmilarity 'matrix' for bases
	 * @param a
	 * @param b
	 * @return
	 */
	int simmilarity(int i, int j) {
		if( a[i - 1] != b[j - 1] ) return missMatch;
		return match;
	}

}
