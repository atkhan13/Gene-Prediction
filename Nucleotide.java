//Simple class to represent a nucleotide which stores both what the nucleotide is and the state it is
//this is used in the sequencePath algorithm in Viterbi.java to easily determine which regions of a sequence are 
// coding vs non-coding
public class Nucleotide {
	public char nucleotide;
	public char state;
	
	public Nucleotide(char nucleotide, char state)
	{
		this.nucleotide = nucleotide;
		this.state = state;
	}
}
