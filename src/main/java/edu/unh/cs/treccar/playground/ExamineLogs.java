package edu.unh.cs.treccar.playground;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

public class ExamineLogs {
	static final String LOG_DIR = "/home/sumanta/Documents/new_research/unh/test200-v1.4results/detective_results";

	public void printSectionsMapped(String logPath){
		BufferedReader br;
		try{
			br = new BufferedReader(new FileReader(logPath));
			String line;
			while((line = br.readLine()) != null){
				if(line.contains("paras mapped")){
					System.out.println(line.split(" ")[0]+" "+line.split(" ")[1]);
				}
			}
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		ExamineLogs el = new ExamineLogs();
		System.out.println("Topic model sections mapped");
		el.printSectionsMapped(ExamineLogs.LOG_DIR+"/log1");
		System.out.println("\n\nKMeans sections mapped");
		el.printSectionsMapped(ExamineLogs.LOG_DIR+"/log2");
	}

}
