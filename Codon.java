//Used to store the emission probabilities
public class Codon {
	private String codon;
	private double prob;
	
	public Codon(String codon, double prob) {
		this.codon = codon;
		this.prob = prob;
	}
	public String getCodon() {
		return this.codon; 
	}
	public double getProb() {
		return this.prob;
	}

}
