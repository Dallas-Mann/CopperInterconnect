import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Scanner;

import org.ejml.data.Complex64F;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.EigenDecomposition;

public class Interconnects{
	private int startNode;
	private int size;
	private double length;
	private double riseTime;
	private double RParam;
	private double GParam;
	private double LModifier;
	private DenseMatrix64F LParams;
	private double CModifier;
	private DenseMatrix64F CParams;
	private int n;

	//EntryPoint
	public static void main(String[] args){
		Interconnects interconnects = new Interconnects(args);
	}
	
	public Interconnects(String[] args) {
		this.readParams(args);
		this.getNumSections();
		try{
			this.construct();
			System.out.println("Done!");
		}catch(FileNotFoundException e){
			System.out.println("Unable to construct netlist.");
		}
	}

	private void readParams(String[] args){
		System.out.println("Reading Parameters");
		try{
			String fileName = args[0];
			startNode = Integer.parseInt(args[1]);
			length = convert(args[2]);
			riseTime = convert(args[3]);
			
			Scanner fileReader = new Scanner(new FileReader(new File(fileName)));
			//read in numberOfInterconnects
			size = Integer.parseInt(fileReader.nextLine());
			//read in resistance per meter
			RParam = Double.parseDouble(fileReader.nextLine());
			//read in conductance per meter
			GParam = Double.parseDouble(fileReader.nextLine());
				
			LParams = new DenseMatrix64F(size, size);//new double[size][size];
			CParams = new DenseMatrix64F(size, size);;
			
			fileReader.nextLine();		
			LModifier = convert(fileReader.nextLine());
			//Read in LParams
			for(int row = 0; row < size; row++){
				for(int col = 0; col < size; col++){
					LParams.set(row, col, LModifier * fileReader.nextDouble());
				}
				fileReader.nextLine();
			}
			
			fileReader.nextLine();
			CModifier = convert(fileReader.nextLine());
			//Read in CParams
			for(int row = 0; row < size; row++){
				for(int col = 0; col < size; col++){
					CParams.set(row, col, CModifier * fileReader.nextDouble());
				}
				//fileReader.nextLine();
			}
			fileReader.close();
		}
		catch(Exception e){
			System.out.println("Invalid file format");
			System.exit(-1);
		}
	}
	
	private void getNumSections(){
		EigenDecomposition<DenseMatrix64F> LEig = DecompositionFactory.eig(size, true);
		EigenDecomposition<DenseMatrix64F> CEig = DecompositionFactory.eig(size, true);
		LEig.decompose(LParams);
		CEig.decompose(CParams);
		
		double LEigMax = 0;
		double CEigMax = 0;
		for (int i = 0; i < size; i++) {
			Complex64F l = LEig.getEigenvalue(i);
			Complex64F c = CEig.getEigenvalue(i);
			double LEigMag = magnitude(l);
			double CEigMag = magnitude(c);
			
			if(LEigMag > LEigMax){
				LEigMax = LEigMag;
			}
			
			if(CEigMag > CEigMax){
				CEigMax = CEigMag;
			}
		}
		//TODO make this use max freq instead of riseTime
		n = (int) Math.round(20*Math.sqrt(LEigMax * CEigMax)*length/riseTime);
	}
	
	private double magnitude(Complex64F l){
		return Math.sqrt(Math.pow(l.real, 2) + Math.pow(l.imaginary, 2));
	}

	private static double convert(String token){
		double baseNum;
		int indexOfModifier = token.length();
		for(int i = 0; i < token.length(); i++){
			char temp = token.charAt(i);
			if(temp != 'E' && Character.isAlphabetic(temp)){
				indexOfModifier = i;
				break;
			}
		}
		String modifier = token.substring(indexOfModifier);
		String value = token.substring(0, indexOfModifier);
		baseNum = Double.parseDouble(value);
		//there should be a trailing modifier after the number
		/*
		F	E-15	femto
		P	E-12	pico
		N	E-9		nano
		U	E-6		micro
		M	E-3		milli
		K	E+3		kilo
		MEG E+6 	mega
		G 	E+9 	giga
		T 	E+12 	tera
		 */
		switch(modifier){
			case "f":
				return baseNum *= Math.pow(10, -15);
			case "p":
				return baseNum *= Math.pow(10, -12);
			case "n":
				return baseNum *= Math.pow(10, -9);
			case "u":
				return baseNum *= Math.pow(10, -6);
			case "m":
				return baseNum *= Math.pow(10, -3);
			case "k":
				return baseNum *= Math.pow(10, 3);
			case "meg":
				return baseNum *= Math.pow(10, 6);
			case "g":
				return baseNum *= Math.pow(10, 9);
			case "t":
				return baseNum *= Math.pow(10, 12);
			default:
				try{
					if(token.chars().allMatch(Character::isDigit)){
						baseNum = Double.parseDouble(token);
					}
				}
				catch(Exception e){
					System.out.println("Invalid Modifier");
				}
				return baseNum;	
		}
	}
	
	private void construct() throws FileNotFoundException{
		double RPS = RParam*length/n;
		double RGPS = 1.0/(GParam*length/n);
		DenseMatrix64F LPS = new DenseMatrix64F(n,n);
		DenseMatrix64F CPS = new DenseMatrix64F(n,n);
		//create per section params
		for(int row = 0; row < size; row++){
			for(int col = 0; col < size; col++){
				LPS.set(row, col, length * LParams.get(row, col) / n);
				CPS.set(row, col, length * CParams.get(row, col) / n);
			}
		}
		System.out.println("Constructing: " + size + " interconnects of "+ n + " sections");
		
		//interconnect #, section #, component #, 2 nodes
		int[][][][] sectionComponents;
		int currentNode = startNode;
		int numberCoupled = n*(n-1)/2;
		int[][][] coupledComponents = new int[numberCoupled][2][2];
		int[] sectionStartNodes = new int[n];
		int[][][] sectionInductanceNodes = new int[size][n][2];
		
		//TODO implement conductance to ground being zero or non-zero
		/*
		if(RGPS != 0 && Double.isFinite(RGPS)){
			sectionComponents = new int[size][n][4+n-1][2];
		}
		else{
			sectionComponents = new int[size][n][3+n-1][2];
		}
		*/
		
		sectionComponents = new int[size][n][3+n-1][2];
		for(int cI = 0; cI < size; cI++){
			for(int cS = 0; cS < n; cS++){
				//set sectionStartNodes so coupled components know which nodes to connect to
				if(cS == 0){
					sectionStartNodes[cI] = currentNode;
				}
				//conductance to ground
				sectionComponents[cI][cS][0][0] = currentNode;
				sectionComponents[cI][cS][0][1] = 0;
				
				//capacitance to ground
				sectionComponents[cI][cS][1][0] = currentNode;
				sectionComponents[cI][cS][1][1] = 0;
				
				//resistance per section
				sectionComponents[cI][cS][2][0] = currentNode;
				sectionComponents[cI][cS][2][1] = ++currentNode;
				
				//inductance per section
				sectionComponents[cI][cS][3][0] = currentNode;
				sectionComponents[cI][cS][3][1] = ++currentNode;
				sectionInductanceNodes[cI][cS][0] = sectionComponents[cI][cS][2][0];
				sectionInductanceNodes[cI][cS][1] = sectionComponents[cI][cS][2][1];
				
				//vcvs per section
				for(int vs = 0; vs < size-1; vs++){
					sectionComponents[cI][cS][vs+4][0] = currentNode;
					sectionComponents[cI][cS][vs+4][1] = ++currentNode;
				}
			}
			//increment currentNode to separate interconnects
			currentNode++;
		}
		
		//n(n-1)/2 coupled connections
		int indexOfCoupled = 0;
		for(int current = 0; current < size-1; current++){
			for(int next = current+1; next < size; next++){
				//coupled resistance
				coupledComponents[indexOfCoupled][0][0] = sectionStartNodes[current];
				coupledComponents[indexOfCoupled][0][1] = sectionStartNodes[next];
				
				//coupled capacitance
				coupledComponents[indexOfCoupled][1][0] = sectionStartNodes[current];
				coupledComponents[indexOfCoupled][1][1] = sectionStartNodes[next];
				
				indexOfCoupled++;
			}
		}
		
		PrintWriter writer = new PrintWriter(new File("interconnect"));
		//nodes all set up, start creating netlist with values.
		int resCount = 1;
		int capCount = 1;
		int indCount = 1;
		int vcvsCount = 1;
		
		for(int cI = 0; cI < size; cI++){
			for(int cS = 0; cS < n; cS++){
				//TODO conductance to ground
				//writer.println("R"+(resCount++)+" "+sectionComponents[cI][cS][0][0]+" "+sectionComponents[cI][cS][0][1]+" "+RGPS);

				//capacitance to ground
				writer.println("C"+(capCount++)+" "+sectionComponents[cI][cS][1][0]+" "+sectionComponents[cI][cS][1][1]+" "+CPS.get(cI, cI));
				
				//resistance per section
				writer.println("R"+(resCount++)+" "+sectionComponents[cI][cS][2][0]+" "+sectionComponents[cI][cS][2][1]+" "+RPS);
				
				//inductance per section
				writer.println("L"+(indCount++)+" "+sectionComponents[cI][cS][3][0]+" "+sectionComponents[cI][cS][3][1]+" "+LPS.get(cI, cI));
				
				//VCVS per section, black magic happens here.
				//I feel like Linus Torvalds. "you aren't expected to understand this"
				int indexVCVS = 0;
				for(int cS2 = 0; cS2 < size; cS2++){
					if(cI != cS2){
						//inductance per section
						double inductance = LPS.get(cI,cS2)/LPS.get(cS2, cS2);
						writer.println("E"+(vcvsCount++)+" "+sectionComponents[cI][cS2][3][0]+" "+sectionComponents[cI][cS2][3][1]+" "+
						sectionComponents[cI][cS][4+indexVCVS][0]+" "+sectionComponents[cI][cS][4+indexVCVS][1]+" "+inductance);
						indexVCVS++;
					}
				}
			}
			//separate each interconnect with a blank line
			writer.println();
		}
		
		//n(n-1)/2 coupled connections
		for(int current = 0; current < size-1; current++){
			for(int next = current+1; next < size; next++){
				//TODO coupled conductance
				//writer.println("R"+(resCount++)+" "+sectionStartNodes[current]+" "+sectionStartNodes[next]+" "+RGPS);

				//coupled capacitance
				writer.println("C"+(capCount++)+" "+sectionStartNodes[current]+" "+sectionStartNodes[next]+" "+CPS.get(current, next));
			}
		}
		writer.close();
	}
}
