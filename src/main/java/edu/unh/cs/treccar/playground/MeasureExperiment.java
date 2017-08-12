package edu.unh.cs.treccar.playground;

import java.util.ArrayList;
import java.util.Collections;
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
			System.out.println();
			ResultForPage r = results.get(pageID);
			ArrayList<ArrayList<String>> clusters = r.getParaClusters();
			ArrayList<ArrayList<String>> currGTClusters = new ArrayList<ArrayList<String>>();
			for(String gtQuery:gt.keySet()){
				if(pageID.equals("Whole%20Corpus")){
					currGTClusters = new ArrayList<ArrayList<String>>(gt.values());
					break;
				}
				if(gtQuery.startsWith(pageID))
					currGTClusters.add(gt.get(gtQuery));
			}
			
			int[][] contingencyMatrix = new int[currGTClusters.size()][clusters.size()];
			//double randIndex = 0.0;
			ArrayList<String> correctParas = new ArrayList<String>();
			ArrayList<String> candParas = new ArrayList<String>();
			ArrayList<String> allParas = new ArrayList<String>();
			for(ArrayList<String> cl:clusters)
				allParas.addAll(cl);
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
			clusteringInfo(clusters, currGTClusters, allParas);
			System.out.println("Contingency table for "+pageID);
			printContingencyMatrix(contingencyMatrix);
			Double rand = computeRand(contingencyMatrix);
			if(rand.isNaN()){
				System.out.println("Adjusted Rand index could not be computed!");
				pageRANDMap.put(pageID, -99.0);
			} else{
				pageRANDMap.put(pageID, rand);
			}
		}
		return pageRANDMap;
	}
	public HashMap<String, Integer> paraLabelMap(ArrayList<ArrayList<String>> clusters){
		HashMap<String, Integer> map = new HashMap<String, Integer>();
		for(int label=0; label<clusters.size(); label++){
			for(String p:clusters.get(label))
				map.put(p, label);
		}
		return map;
	}
	public void clusteringInfo(ArrayList<ArrayList<String>> clusters, 
			ArrayList<ArrayList<String>> gt, ArrayList<String> paras){
		HashMap<String, Integer> clmap = paraLabelMap(clusters);
		HashMap<String, Integer> gtmap = paraLabelMap(gt);
		int tp=0, tn=0, fp=0, fn=0;
		for(int i=0; i<paras.size()-1; i++){
			for(int j=i+1; j<paras.size(); j++){
				if(clmap.get(paras.get(i))==clmap.get(paras.get(j)) && gtmap.get(paras.get(i))==gtmap.get(paras.get(j)))
					tp++;
				else if(clmap.get(paras.get(i))!=clmap.get(paras.get(j)) && gtmap.get(paras.get(i))!=gtmap.get(paras.get(j)))
					tn++;
				else if(clmap.get(paras.get(i))==clmap.get(paras.get(j)) && gtmap.get(paras.get(i))!=gtmap.get(paras.get(j)))
					fp++;
				else if(clmap.get(paras.get(i))!=clmap.get(paras.get(j)) && gtmap.get(paras.get(i))==gtmap.get(paras.get(j)))
					fn++;
			}
		}
		System.out.println("No of paras: "+paras.size()+", True +ve: "+tp+", True -ve: "+tn+", False +ve: "+fp+", False -ve: "+fn);
		System.out.println("RAND = "+((double)(tp+tn))/(tp+tn+fp+fn));
	}
	public HashMap<String, Double> calculatePurityPerPage(){
		return this.calculatePurityPerPage(this.getResultFromExperiment(), this.getGtMap());
	}
	public HashMap<String, Double> calculatePurityPerPage(HashMap<String, ResultForPage> results, 
			HashMap<String, ArrayList<String>> gt){
		HashMap<String, Double> pagePurityMap = new HashMap<String, Double>();
		for(String pageID:results.keySet()){
			ResultForPage r = results.get(pageID);
			ArrayList<ArrayList<String>> clusters = r.getParaClusters();
			ArrayList<ArrayList<String>> currGTClusters = new ArrayList<ArrayList<String>>();
			for(String gtQuery:gt.keySet()){
				if(pageID.equals("Whole%20Corpus")){
					currGTClusters = new ArrayList<ArrayList<String>>(gt.values());
					break;
				}
				if(gtQuery.startsWith(pageID))
					currGTClusters.add(gt.get(gtQuery));
			}
			Double purity = computePurity(clusters, currGTClusters);
			if(purity.isNaN()){
				System.out.println("Purity could not be computed");
				pagePurityMap.put(pageID, purity);
			} else
				pagePurityMap.put(pageID, purity);
		}
		return pagePurityMap;
	}
	private double computeRand(int[][] contMat){
		double score = 0.0;
		int sumnij=0, sumni=0, sumnj=0, nC2=0, nrow=0, ncol=0, n=0;
		int a=0, b=0, c=0, d=0;
		ncol = contMat[0].length;
		nrow = contMat.length;
		int[] nivals = new int[nrow];
		int[] njvals = new int[ncol];
		//nC2 = this.nC2(ncol+nrow);
		for(int r=0; r<nrow; r++){
			for(int co=0; co<ncol; co++)
				n+=contMat[r][co];
		}
		nC2 = this.nC2(n);
		for(int i=0; i<nrow; i++){
			int ni=0;
			for(int j=0; j<ncol; j++){
				sumnij+=this.nC2(contMat[i][j]);
				ni+=contMat[i][j];
				njvals[j]+=contMat[i][j];
			}
			nivals[i]=ni;
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
		score = nom/denom;
		System.out.println("n: "+n+", nom: "+nom+", denom: "+denom+", Adjusted RAND = "+score);
		return score;
	}
	public double computePurity(ArrayList<ArrayList<String>> candidateClusters, ArrayList<ArrayList<String>> gtClusters){
		String result = "";
		int n=0, sum=0, tempMaxCount, maxCount;
		double score=0.0;
		for(ArrayList<String> currCandCluster:candidateClusters){
			maxCount=0;
			n+=currCandCluster.size();
			for(ArrayList<String> currGTCluster:gtClusters){
				tempMaxCount = 0;
				for(String candPara:currCandCluster){
					if(currGTCluster.contains(candPara))
						tempMaxCount++;
				}
				if(tempMaxCount>maxCount)
					maxCount = tempMaxCount;
			}
			sum+=maxCount;
		}
		score = ((double)sum)/n;
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
