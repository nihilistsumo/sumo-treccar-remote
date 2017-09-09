package edu.unh.cs.treccar.playground;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;

public class ExamineLogs {
	static final String LOG_DIR = "/home/sumanta/Documents/new_research/unh/test200-v1.4results/detective_results2";

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
	public void measureSpreadContMat(String logPath){
		BufferedReader br;
		HashMap<String, double[]> pageSpreadMap = new HashMap<String, double[]>();
		try{
			br = new BufferedReader(new FileReader(logPath));
			String line, pageid;
			int row, col;
			int[][] contmat, invContmat;
			double rand = 0, adrand = 0;
			while((line = br.readLine()) != null){
				if(line.startsWith("RAND ="))
					rand = new Double(line.split(" = ")[1]);
				else if(line.contains("Contingency table for")){
					pageid = line.split(" ")[3];
					System.out.println(pageid);
					ArrayList<String> contmatStrings = new ArrayList<String>();
					while(!((line = br.readLine()).startsWith("n:"))){
						System.out.println(line);
						contmatStrings.add(line);
					}
					adrand = new Double(line.split(" = ")[1]);
					row = contmatStrings.size();
					col = contmatStrings.get(0).split(" ").length;
					contmat = new int[row][col];
					invContmat = new int[col][row];
					for(int i=0; i<row; i++){
						String[] controw = contmatStrings.get(i).split(" ");
						for(int j=0; j<col; j++)
							contmat[i][j] = Integer.parseInt(controw[j]);
					}
					for(int i=0; i<row; i++){
						for(int j=0; j<col; j++)
							invContmat[j][i] = contmat[i][j];
					}
					double[] rowSpreadScores = new double[row];
					double[] colSpreadScores = new double[col];
					for(int i=0; i<row; i++)
						rowSpreadScores[i] = calculateSpread(contmat[i]);
					for(int j=0; j<col; j++)
						colSpreadScores[j] = calculateSpread(invContmat[j]);
					System.out.println("Row spread: ");
					for(int i=0; i<row; i++)
						System.out.print(rowSpreadScores[i]+" ");
					System.out.println();
					System.out.println("Column spread: ");
					for(int j=0; j<col; j++)
						System.out.print(colSpreadScores[j]+" ");
					System.out.println();
					double spread = 0, rspread = 0, cspread = 0;
					for(int i=0; i<row; i++)
						rspread+= rowSpreadScores[i];
					for(int j=0; j<col; j++)
						cspread+= colSpreadScores[j];
					rspread = rspread/row;
					cspread = cspread/col;
					spread = (rspread+cspread)/2;
					System.out.println("Spread score: "+spread);
					double[] scores = new double[5];
					scores[0] = rspread;
					scores[1] = cspread;
					scores[2] = spread;
					scores[3] = rand;
					scores[4] = adrand;
					pageSpreadMap.put(pageid, scores);
				}
			}
		}
		catch(Exception e){
			e.printStackTrace();
		}
		for(String pid:pageSpreadMap.keySet()){
			double[] scr = pageSpreadMap.get(pid);
			System.out.println(pid+" "+scr[0]+" "+scr[1]+" "+scr[2]+" "+scr[3]+" "+scr[4]);
		}
	}
	private double calculateSpread(int[] cmValArray){
		double max = 0, sum = 0;
		int n = cmValArray.length;
		for(int i=0; i<n; i++){
			if(cmValArray[i]>max)
				max = cmValArray[i];
			sum+= cmValArray[i];
		}
		if(sum==0)
			return 0;
		return (n*max-sum)/(n*sum-max);
	}
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		ExamineLogs el = new ExamineLogs();
		el.measureSpreadContMat(ExamineLogs.LOG_DIR+"/ldakl/ldakl_logs");
		/*
		System.out.println("Topic model sections mapped");
		el.printSectionsMapped(ExamineLogs.LOG_DIR+"/log1");
		System.out.println("\n\nKMeans sections mapped");
		el.printSectionsMapped(ExamineLogs.LOG_DIR+"/log2");
		*/
	}

}
