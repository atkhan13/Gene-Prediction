//Used in viterbiResult function to define one cell in the matrix. This is so that various information about
//each cell can be stored such as the nucleotide, the state of the cell, the maximum probability found and 
//which state before it led to the maximum probability
public class Cell {
	public char nucleotide;
	public char state;
	public double prob;
	public char prevState;
	
	public Cell(char nucleotide, char state)
	{
		this.nucleotide = nucleotide;
		this.state = state;
	}
	public void setProb(double prob) {
		this.prob = prob;
	}
	public void setPrevState(char prevState) {
		this.prevState = prevState;
	}
	

}
