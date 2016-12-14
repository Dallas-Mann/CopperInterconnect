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
		System.out.println("Constructing: " + n + " Sections");
		PrintWriter writer = new PrintWriter(new File("interconnect"));
		
		int[] beginningNodes = new int[size];
		
		int currentNode = startNode;
		
		int resCount = 1;
		int capCount = 1;
		int indCount = 1;
		int vcvsCount = 1;
		for(int cI = 0; cI < size; cI++){
			beginningNodes[cI] = currentNode;
			for(int currentSec = 0; currentSec < n; currentSec++){
				//TODO make if elseif else to catch case of both being zero, no time right now
				if(RGPS != 0 && Double.isFinite(RGPS)){
					writer.println("R"+(resCount++)+" "+currentNode+" 0 "+RGPS);
				}
				if(CPS.get(cI, cI) != 0){
					writer.println("C"+(capCount++)+" "+currentNode+" 0 "+CPS.get(cI, cI));
				}
				if(RPS != 0){
					writer.println("R"+(resCount++)+" "+(currentNode++)+" "+currentNode+" "+RPS);
				}
				if(LPS.get(cI, cI) != 0){
					writer.println("L"+(indCount++)+" "+(currentNode++)+" "+currentNode+" "+LPS.get(cI, cI));
				}
				if(size > 1){
					//TODO finish VCVS sources
					//for loop to insert vcvs from coupling
					for(int numVCVS = 1; numVCVS < size; numVCVS++){
						writer.println("E"+(vcvsCount++)+" "+(currentNode++)+" "+currentNode+" "+LPS.get(cI, cI));
					}
				}
			}
		}
		
		if(size > 1){
			//TODO finish coupling of capacitance/conductance
			//couple the interconnects
			
		}
		writer.close();
	}
}
