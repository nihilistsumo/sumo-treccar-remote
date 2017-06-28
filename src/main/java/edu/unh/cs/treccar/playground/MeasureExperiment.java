package edu.unh.cs.treccar.playground;

import java.util.ArrayList;
import java.util.HashMap;

public class MeasureExperiment {
	private final HashMap<String, ResultForPage> resultFromExperiment;
	private final HashMap<String, ArrayList<String>> gtMap;
	public HashMap<String, ResultForPage> getResultFromExperiment() {
		return resultFromExperiment;
	}
	public MeasureExperiment(HashMap<String, ResultForPage> result, HashMap<String, ArrayList<String>> gtMap){
		this.resultFromExperiment = result;
		this.gtMap = gtMap;
	}
	public HashMap<String, ArrayList<String>> getGtMap() {
		return gtMap;
	}
	public HashMap<String, Double> calculateRANDPerPage(){
		return calculateRANDPerPage(this.getResultFromExperiment(), this.getGtMap());
	}
	public HashMap<String, Double> calculateRANDPerPage(HashMap<String, ResultForPage> results, 
			HashMap<String, ArrayList<String>> gt){
		HashMap<String, Double> pageRANDMap = new HashMap<String, Double>();
		for(String pageID:results.keySet()){
			ResultForPage r = results.get(pageID);
			ArrayList<ArrayList<String>> clusters = r.getParaClusters();
			ArrayList<ArrayList<String>> currGTClusters = new ArrayList<ArrayList<String>>();
			for(String gtQuery:gt.keySet()){
				if(gtQuery.startsWith(pageID))
					currGTClusters.add(gt.get(gtQuery));
			}
			
			String resultString = "";
			
			
			int[][] contingencyMatrix = new int[currGTClusters.size()][clusters.size()];
			//double randIndex = 0.0;
			ArrayList<String> correctParas = new ArrayList<String>();
			ArrayList<String> candParas = new ArrayList<String>();
			for(int i=0; i<currGTClusters.size(); i++){
				for(int j=0; j<clusters.size(); j++){
					int matchCount = 0;
					correctParas = currGTClusters.get(i);
					candParas = clusters.get(j);
					if(correctParas == null){
						System.out.println("We have null in correctParas!");
					} else if(candParas != null){
						for(String candPara : candParas){
							if(correctParas.contains(candPara)){
								matchCount++;
							}
						}
					}
					contingencyMatrix[i][j] = matchCount;
				}
			}
			printContingencyMatrix(contingencyMatrix);
			Double rand = computeRand(contingencyMatrix);
			if(rand.isNaN()){
				System.out.println("Adjusted Rand index could not be computed!");
				pageRANDMap.put(pageID, null);
			} else{
				pageRANDMap.put(pageID, rand);
			}
		}
		return pageRANDMap;
	}
	private double computeRand(int[][] contMat){
		double score = 0.0;
		int sumnij=0, sumni=0, sumnj=0, nC2=0, nrow=0, ncol=0, n=0;		
		ncol = contMat[0].length;
		nrow = contMat.length;
		int[] njvals = new int[ncol];
		//nC2 = this.nC2(ncol+nrow);
		for(int r=0; r<nrow; r++){
			for(int c=0; c<ncol; c++)
				n+=contMat[r][c];
		}
		nC2 = this.nC2(n);
		for(int i=0; i<nrow; i++){
			int ni=0;
			for(int j=0; j<ncol; j++){
				sumnij+=this.nC2(contMat[i][j]);
				ni+=contMat[i][j];
				njvals[j]+=contMat[i][j];
			}
			sumni+=this.nC2(ni);
		}
		for(int j=0; j<njvals.length; j++){
			sumnj+=this.nC2(njvals[j]);
		}
		
		/* ################### 
		 * This code is for simple Rand index
		
		int a=0, b=0, c=0, d=0;
		a = sumnij;
		b = sumni - sumnij;
		c = sumnj - sumnij;
		d = nC2 - (a+b+c);
		score = ((double)(a+d))/(a+b+c+d);
		
		#################### */
		
		double denom = ((double)(sumni+sumnj))/2-((double)sumni*sumnj/nC2);
		double nom = (sumnij-((double)(sumni*sumnj))/nC2);
		System.out.println("n: "+n+", sumnij: "+sumnij+", sumni: "+sumni+", sumnj: "+sumnj+", nC2: "+nC2+", nom: "+nom+", denom: "+denom);
		score = nom/denom;
		return score;
	}
	private int nC2(int n){
		if(n<2) return 0;
		else if(n==2) return 1;
		else{
			return n*(n-1)/2;
		}
	}
	private void printContingencyMatrix(int[][] contingency){
		int colNum = contingency[0].length;
		int rowNum = contingency.length;
		for(int i=0; i<rowNum; i++){
			for(int j=0; j<colNum; j++){
				System.out.print(contingency[i][j]+" ");
			}
			System.out.println();
		}
	}
}
