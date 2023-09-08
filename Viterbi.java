import java.io.*;
import java.util.ArrayList;
import java.util.Scanner;
public class Viterbi {
	public static void main(String [] args)
	{
		//Load both the sequences and the file containing the emission and transition probabilities
		File sequences = new File("Vibrio_vulnificus.ASM74310v1.dna.nonchromosomal.fa");
		File info = new File("Details_for_Vitebri.txt");
		try {
			//Call Viterbi algorithm to write to the results.gff3 file which is created and gets updated
			File answer = viterbi(sequences, info);
			if(answer.exists()) {
				System.out.println("Results file updated. View file to see gene prediction using Viterbi algorithm");
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static File viterbi(File sequences, File config) throws Exception
	{
		//This function only sets up the emission and transition matrix from the file inputted
		//It then goes through every sequence from the Fasta file given and executes the Viterbi algorithm
		// by calling the viterbiResult function below. viterbiResult gives the matrix, and then the sequenceStates 
		// function is called to associate each nucleotide within the sequence with a char being either I 
		//(for intergenic) T (for start) G(for genic) or P(for stop). The program then goes through the
		//sequence and notes the nucleotides which contain start genic or stop codons then writes their position
		//to a file which is returned as output
		
		Scanner sc = null;
		try {
			sc = new Scanner(sequences);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		ArrayList<String> mySequences = new ArrayList<String>();
		ArrayList<String> notation = new ArrayList<String>();
		String seq = sc.nextLine();
		notation.add(seq);
		seq = sc.nextLine();
		String thisSeq = "";
		//Code below sets up the emission and transition probability matrix. We convert probabilities to log 
		while(true)
		{
			if(seq.contains(">"))
			{
				notation.add(seq);
				mySequences.add(thisSeq);
				thisSeq = "";
				seq = sc.nextLine();
				continue;
			}
			if(!sc.hasNext()) {
				thisSeq += seq;
				mySequences.add(thisSeq);
				break;
			}
			thisSeq += seq;
			seq = sc.nextLine();
		}
		Scanner sc2 = null;
		try {
			sc2 = new Scanner(config);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		String intergenicLength = sc2.nextLine();
		intergenicLength = sc2.nextLine();
		String genicLength = sc2.nextLine();
		genicLength = sc2.nextLine();
		genicLength = sc2.nextLine();
		double inLen = Double.parseDouble(intergenicLength);
		double genLen = Double.parseDouble(genicLength);
		double[][] transition = new double[4][4];
		transition[0][0] = Math.log10((double)(inLen - 1)/inLen);
		transition[0][1] = Math.log10(1.0/inLen);
		transition[0][2] = Math.log10(0.0);
		transition[0][3] = Math.log10(0.0);
		transition[1][0] = Math.log10(0.0);
		transition[1][1] = Math.log10(0.0);
		transition[1][2] = Math.log10(1);
		transition[1][3] = Math.log10(0.0);
		transition[2][0] = Math.log10(0.0); 
		transition[2][1] = Math.log10(0.0);
		transition[2][2] = Math.log10((double)(genLen-1)/genLen);
		transition[2][3] = Math.log10(1.0/genLen);
		transition[3][0] = Math.log10(1);
		transition[3][1] = Math.log10(0.0);
		transition[3][2] = Math.log10(0.0);
		transition[3][3] = Math.log10(0.0);
		double[] emissionIn = new double[4];
		String curr = sc2.nextLine();
		curr = sc2.nextLine();
		curr = sc2.next();
		int i = 0;
		while(!curr.contains(">"))
		{
			String nuc = sc2.next();
			double nucleotide = Math.log10(Double.parseDouble(nuc));
			emissionIn[i] = nucleotide;
			i++;
			curr = sc2.nextLine();
			curr = sc2.next();
		}
		Codon[] emissionStart = new Codon[64];
		curr = sc2.nextLine();
		curr = sc2.next();
		i = 0;
		while(!curr.contains(">")){
			String codon = curr;
			curr = sc2.next();
			double prob = Math.log10(Double.parseDouble(curr));
			emissionStart[i] = new Codon(codon, prob);
			curr = sc2.next();
			i++;
		}
		Codon[] emissionStop = new Codon[64];
		curr = sc2.nextLine();
		curr = sc2.next();
		i = 0;
		while(!curr.contains(">")) {
			String codon = curr;
			curr = sc2.next();
			double prob = Math.log10(Double.parseDouble(curr));
			emissionStop[i] = new Codon(codon, prob);
			curr = sc2.next();
			i++;
		}
		Codon[] emissionGenic = new Codon[64];
		curr = sc2.nextLine();
		curr = sc2.next();
		i = 0;
		while(!curr.contains(">")) {
			String codon = curr;
			curr = sc2.next();
			double prob = Math.log10(Double.parseDouble(curr));
			emissionGenic[i] = new Codon(codon, prob);
			if(!sc2.hasNext()) {
				break;
			}
			curr = sc2.next();
			i++;
		}
		double [] initial = new double[4];
		initial[0] = Math.log10(1.0);
		initial[1] = 0.0;
		initial[2] = 0.0;
		initial[3] = 0.0;
		
		//code below calls viterbiResult on each sequence and sequenceStates to get the path for each sequence
		//and notes down which regions are genes within the sequence
		try {
			File file = new File( System.getProperty("user.dir")+"/results.gff3");
			boolean result;
			result = file.createNewFile();
			if(result)      // test if successfully created a new file  
			{  
			System.out.println("file created "+file.getCanonicalPath()); //returns the path string  
			}  
			FileWriter myAnswer = new FileWriter(file);
			myAnswer.write("##gff-version 3\n");
			for(int j = 0; j < mySequences.size(); j++) {
				Cell[][] matrix = vitebriResult(mySequences.get(j), transition, emissionIn, emissionStart, emissionStop, emissionGenic, initial);
				Nucleotide [] sequence = sequenceStates(matrix);
				String[] note = notation.get(j).split(" ");
				String[] toWrite = note[0].split(">");
				int start = 0;
				for(int x = 0; x < sequence.length; x++) {
					if(sequence[x].state != 'I') {
						start = x + 1;
						break;
					}
					
				}
				int stop = 0;
				boolean gene = true;
				for(int y = start; y < sequence.length; y++) {
					if(gene) {
						if(sequence[y].state == 'I') {
							stop = y;
							gene = false;
							if(start == 0 && stop == 0) {
								continue;
							}
							myAnswer.write(toWrite[1]+"\tena\tCDS\t"+start+"\t"+stop+"\t.\t+\0\t.\n");
						}
					}
					else {
						if(sequence[y].state != 'I') {
							start = y + 1;
							gene = true;
						}
					}
				}
			}
			myAnswer.close();
			return file;
		} catch (IOException e) {
			e.printStackTrace();
		}
		throw new Exception("File not made");

		
	}
	
	//Program which actually executes the Viterbi algorithm. Note, emission and transition probabilities converted to logs so we use addition 
	//for the probability instead of multiplication
	public static Cell[][] vitebriResult(String sequence, double[][] transition, double[] emIn, Codon[] emStart, Codon[] emStop, Codon [] emGen, double[] init)
	{
		//Initialize the matrix which will be returned
		Cell[][] myMatrix = new Cell[4][sequence.length()];
		double prob = 0.0;
		//First initialization step for the first nucleotide in the sequence
		if(sequence.charAt(0) == 'A'){
			prob = emIn[0] + init[0];
			myMatrix[0][0] = new Cell('A', 'I');
			myMatrix[0][0].setProb(prob);
			myMatrix[1][0] = new Cell('A', 'T');
			myMatrix[1][0].setProb(Math.log(0.0));
			myMatrix[2][0] = new Cell('A', 'G');
			myMatrix[2][0].setProb(Math.log(0.0));
			myMatrix[3][0] = new Cell('A', 'P');
			myMatrix[3][0].setProb(Math.log(0.0));
		}
		else if(sequence.charAt(0) == 'T') {
			prob = emIn[1] + init[0];
			myMatrix[0][0] = new Cell('T', 'I');
			myMatrix[0][0].setProb(prob);
			myMatrix[1][0] = new Cell('T', 'T');
			myMatrix[1][0].setProb(Math.log(0.0));
			myMatrix[2][0] = new Cell('T', 'G');
			myMatrix[2][0].setProb(Math.log(0.0));
			myMatrix[3][0] = new Cell('T', 'P');
			myMatrix[3][0].setProb(Math.log(0.0));
		}
		else if(sequence.charAt(0) == 'C') {
			prob = emIn[2] + init[0];
			myMatrix[0][0] = new Cell('C', 'I');
			myMatrix[0][0].setProb(prob);
			myMatrix[1][0] = new Cell('C', 'T');
			myMatrix[1][0].setProb(Math.log(0.0));
			myMatrix[2][0] = new Cell('C', 'G');
			myMatrix[2][0].setProb(Math.log(0.0));
			myMatrix[3][0] = new Cell('C', 'P');
			myMatrix[3][0].setProb(Math.log(0.0));
		}
		else if(sequence.charAt(0) == 'G'){
			prob = emIn[3] + init[0];
			myMatrix[0][0] = new Cell('G', 'I');
			myMatrix[0][0].setProb(prob);
			myMatrix[1][0] = new Cell('G', 'T');
			myMatrix[1][0].setProb(Math.log(0.0));
			myMatrix[2][0] = new Cell('G', 'G');
			myMatrix[2][0].setProb(Math.log(0.0));
			myMatrix[3][0] = new Cell('G', 'P');
			myMatrix[3][0].setProb(Math.log(0.0));
		}
		//We must find the probabilities for each state for each nucleotide in the sequence
		for(int i = 1; i < sequence.length(); i++)
		{
			//determine the emission probability for the nucleotide at position i given that the state is intergenic
			double em = 0.0;
			if(sequence.charAt(i) == 'A'){
				em = emIn[0];
			}
			else if(sequence.charAt(i) == 'T'){
				em = emIn[1];
			}
			else if(sequence.charAt(i) == 'C') {
				em = emIn[2];
			}
			else if(sequence.charAt(i) == 'G') {
				em = emIn[3];
			}
			//Find the max probability based on if the previous state was intergenic or stop (transition probability
			//from start to intergenic or from genic to intergenic is always 0 so does not need to be calculated since
			//it will never be the max
			double probI = em + transition[0][0] + myMatrix[0][i-1].prob;
			double probP = em + transition[3][0] + myMatrix[3][i-1].prob;
			double probIn = Math.max(probI, probP);
			myMatrix[0][i] = new Cell(sequence.charAt(i), 'I');
			myMatrix[0][i].setProb(probIn);
			//Assign the previous state of the nucleotide in question for the intergenic state as whichever state
			//gave the maximum probability
			if(probI == probIn)
			{
				myMatrix[0][i].setPrevState('I');
			}
			else if(probP == probIn) {
				myMatrix[0][i].setPrevState('P');
			}
			//For start, genic and stop we must look 3 nucleotides back so we set the first three nucleotides to have
			// zero probability for these states.
			if(i < 3)
			{
				myMatrix[1][i] = new Cell(sequence.charAt(i), 'T');
				myMatrix[1][i].setProb(Math.log(0.0));
				myMatrix[2][i] = new Cell(sequence.charAt(i), 'G');
				myMatrix[2][i].setProb(Math.log(0.0));
				myMatrix[3][i] = new Cell(sequence.charAt(i), 'P');
				myMatrix[3][i].setProb(Math.log(0.0));
			}
			else {
				
				//Find and assign each state for nucleotide i the max probability and note down which previous state 
				//led to this max probability which can be used later to determine the path
				String currCodon = new StringBuilder().append(sequence.charAt(i-2)).append(sequence.charAt(i-1)).append(sequence.charAt(i)).toString();
				for(int j = 0; j < emStart.length; j++)
				{
					if(emStart[j].getCodon().equals(currCodon)) {
						em = emStart[j].getProb(); 
						break;
					}
				}
				probI = em + transition[0][1] + myMatrix[0][i-3].prob;
				double probStart = probI;
				myMatrix[1][i] = new Cell(sequence.charAt(i), 'T');
				myMatrix[1][i].setProb(probStart);
				if(probI == probStart)
				{
					myMatrix[1][i].setPrevState('I');
				}
				for(int j = 0; j < emGen.length; j++) {
					if(emGen[j].getCodon().equals(currCodon)) {
						em = emGen[j].getProb();
						break;
					}
				}
				double probT = em + transition[1][2] + myMatrix[1][i-3].prob;
				double probG = em + transition[2][2] + myMatrix[2][i-3].prob;
				double probGen = Math.max(probT, probG);
				myMatrix[2][i] = new Cell(sequence.charAt(i), 'G');
				myMatrix[2][i].setProb(probGen);
				if(probT == probGen) {
					myMatrix[2][i].setPrevState('T');
				}
				else if(probG == probGen) {
					myMatrix[2][i].setPrevState('G');
				}
				for(int j = 0; j < emStop.length; j++)
				{
					if(emStop[j].getCodon().equals(currCodon)) {
						em = emStop[j].getProb();
						break;
					}
				}
				probG = em + transition[2][3] + myMatrix[2][i-3].prob;
				double probStop = probG;
				myMatrix[3][i] = new Cell (sequence.charAt(i), 'P');
				myMatrix[3][i].setProb(probStop);
				if(probG == probStop) {
					myMatrix[3][i].setPrevState('G');
				}
				else if(probP == probStop) {
					myMatrix[1][i].setPrevState('P');
				}
				
			}
	
		}
		return myMatrix;
		
	}
	//This is the traceback step of the algorithm to find the optimal state path
	public static Nucleotide[] sequenceStates(Cell[][] matrix) {
		//Start by finding the maximum probability for the last nucleotide and which state had given this
		//maximum probability
		int lastNucleotide = matrix[0].length-1;
		double maxProb = Math.max(Math.max(matrix[0][lastNucleotide].prob, matrix[1][lastNucleotide].prob), Math.max(matrix[2][lastNucleotide].prob, matrix[3][lastNucleotide].prob));
		Nucleotide[] sequencePath = new Nucleotide[matrix[0].length];
		Cell curr = null;
		int start = lastNucleotide;
		if(maxProb == matrix[0][lastNucleotide].prob) {
			curr = matrix[0][lastNucleotide];
		}
		else if(maxProb == matrix[1][lastNucleotide].prob) {
			curr = matrix[1][lastNucleotide];
		}
		else if(maxProb == matrix[2][lastNucleotide].prob) {
			curr = matrix[2][lastNucleotide];
		}
		else if(maxProb == matrix[3][lastNucleotide].prob) {
			curr = matrix[3][lastNucleotide];
		}
		//Going backwards we assign each nucleotide a state
		while(start >= 0){
			char prevState = curr.prevState;
			char thisState = curr.state;
			if(thisState == 'I') {
				sequencePath[start] = new Nucleotide(curr.nucleotide, curr.state);
			}
			//If the current state we are at is stop, genic, or start, we have to assign this to not only the 
			//nucleotide we are at, but two nucleotides back since in our Viterbi algorithm we looked back three
			//nucleotides as we have to consider codons
			else {
				sequencePath[start] = new Nucleotide(curr.nucleotide, thisState);
				start -=1;
				curr = matrix[0][start];
				sequencePath[start] = new Nucleotide(curr.nucleotide, thisState);
				start-=1;
				curr = matrix[0][start];
				sequencePath[start] = new Nucleotide(curr.nucleotide, thisState);
			}
			if(start == 0) {
				break;
			}
			//We choose the next cell to look at depending on what the previous state of the current
			//nucleotide is (i.e., the previous state that had given the current cell in the matrix its 
			//highest probability)
			start -= 1;
			if(prevState == 'I') {
				curr = matrix[0][start];
			}
			else if(prevState == 'T') {
				curr = matrix[1][start];
			}
			else if(prevState == 'G') {
				curr = matrix[2][start];
			}
			else if(prevState == 'P') {
				curr = matrix[3][start];
			}
			
		}
		return sequencePath;
	}
}
