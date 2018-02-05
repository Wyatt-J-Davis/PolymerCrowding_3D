package org.opensourcephysics.sip.CPM;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.text.DecimalFormat;
import java.util.Calendar;
import java.util.Scanner;

public class DataFile {
	private static Calendar calendar;
	private final static String baseDir = "data";
	private String initialConfigurations;
	private String name;
	private String data;
	private BufferedWriter bw;
	private Path filePath;
	public enum FileIdentifier { 
		SIZE,
		FRACTION,
		SIZE_AND_FRACTION
	};
	
	public DataFile(String name, String configurations, FileIdentifier type, boolean multiple){
		this.name = name;
		if(DataFile.calendar == null){
			calendar = Calendar.getInstance();
		}
		Scanner scan = new Scanner(configurations);
		String phiN = "";
		String runNumber = "";
		String q = "";

		while(scan.hasNext()){
			String token = scan.next();
			if(token.equals("Fraction:")){
				phiN = scan.next();
				double phi_N = Double.parseDouble(phiN);
				DecimalFormat threeDecimal = new DecimalFormat("#0.###");
				DecimalFormat largeDecimal = new DecimalFormat("0.##E0");
				DecimalFormat fourDecimal = new DecimalFormat("#0.####");

				if(phi_N < 0.0001){
					phiN = largeDecimal.format(phi_N);
				} else if(phi_N < 0.001){
					phiN = fourDecimal.format(phi_N);
				} else {
					phiN = threeDecimal.format(phi_N);
				}

			}
			
			if(token.equals("Ratio:")){
				q = scan.next();
			}
			
			if(token.equals("Number:")){
				runNumber = scan.next();
			}
		}
		scan.close();
		
		String prefix;
		switch(type){
			case SIZE: 
				prefix = "phiN=" + phiN; break;
			case FRACTION:
				prefix = "q=" + q; break;
			case SIZE_AND_FRACTION:
				prefix = "phiN=" + phiN + ",q=" + q; break;
			default:
				prefix = "phiN=" + phiN; break;
		}

		String creationTime = calendar.get(Calendar.HOUR_OF_DAY) + "." + 
							  calendar.get(Calendar.MINUTE);
		String creationDate = calendar.get(Calendar.DAY_OF_MONTH) + "-" + 
				  (calendar.get(Calendar.MONTH) + 1) + "-" + // calendar.MONTH starts counting at 0
				  calendar.get(Calendar.YEAR) + " " + creationTime;
		Path dir;
		if(multiple){
			dir = Paths.get(DataFile.baseDir + "/" +  prefix + "_" + creationDate + "/" + runNumber);
		} else {
			dir = Paths.get(DataFile.baseDir + "/");
		}
		
		// create directory to store datafiles, categorized by date if it doesn't exist
		try {
			Files.createDirectories(dir);
		} catch (IOException e) {
			e.printStackTrace();
		}

		
		if(multiple){
			filePath = Paths.get(dir + "/" + prefix + "_" + name + ".dat");
		}else{
			filePath = Paths.get(dir + "/" + prefix + "_" + name + "." + creationDate +".dat");
		}
		
		try {
			bw = Files.newBufferedWriter(filePath, 
					Charset.forName("US-ASCII"), filePath.toFile().exists()? StandardOpenOption.APPEND : StandardOpenOption.CREATE);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		record(configurations);
	}
	
	/**
	 * Adds the data to the buffer.
	 * @param data Data to be added
	 */
	public void record(String data){
		try {
			bw.write(data);
			bw.newLine();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (NullPointerException e){
			System.out.println("record method called without calling constructor!");
		}
	}
	
	/**
	 * Writes the data to the file.
	 */
	public void write(){
		try{
			bw.flush();
		} catch (IOException e){
			e.printStackTrace();
			System.out.println("Error writing to file.");
		}
	}
	
	/**
	 * Closes the stream for writing.
	 */
	public void close(){
		try{
			bw.close();
		} catch(IOException e){
			e.printStackTrace();
			System.out.println("Problem closing the buffer stream.");
		}
	}
}
