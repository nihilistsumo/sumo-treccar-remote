package edu.unh.cs.treccar.playground;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;

public class DifficultySortingCorpus {
	public String[] argsForSingleRun;
	public static int DIFFICULTY_PERC = 10;
	public String k = "0";
	public String numIter = "500";
	public String model = "99";
	public String tw = "0";
	public String alphaSum = "1.3";
	public String beta = "0.3";

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		DifficultySortingCorpus diff = new DifficultySortingCorpus();
		diff.argsForSingleRun = 
				new String[]{diff.k,diff.numIter,diff.model,diff.tw,diff.alphaSum,diff.beta,
        		"/home/sumanta/Documents/new_research/unh/test200-v1.4/all.test200.cbor.paragraphs",
        		"/home/sumanta/Documents/new_research/unh/test200-v1.4/all.test200.cbor.outlines",
        		"/home/sumanta/Documents/new_research/unh/test200-v1.4/all.test200.cbor.toplevel.qrels",
        		"/home/sumanta/Documents/new_research/unh/test200-v1.4/all.test200.cbor.article.qrels",
        		"/home/sumanta/Documents/new_research/unh/test200-v1.4results/custom_lda_and_km_results"};
		SingleRun sr = new SingleRun(diff.argsForSingleRun);
		sr.runExperiment();
		HashMap<String, Double> pagerand = sr.pageRAND;
		LinkedHashMap<String, Double> sortedPageRAND = new LinkedHashMap<String, Double>();
		ArrayList<Double> randValues = new ArrayList<Double>(pagerand.values());
		Collections.sort(randValues);
		Iterator<Double> valItr = randValues.iterator();
		Double rand;
		while(valItr.hasNext()){
			rand = valItr.next();
			for(String pageid:pagerand.keySet()){
				if(Math.abs(pagerand.get(pageid)-rand)<0.00001){
					sortedPageRAND.put(pageid, rand);
					break;
				}
			}
		}
		System.out.println("\n\n\n");
		System.out.println("k="+diff.k+",iter="+diff.numIter+",model="+diff.model+",tw="+diff.tw+",alphaSum="+diff.alphaSum+",beta="+diff.beta);
		System.out.println("##############################");
		System.out.println("Difficulty sorting of pages\n\n");
		if(DifficultySortingCorpus.DIFFICULTY_PERC == 0){
			System.out.println("All sorted pages");
			System.out.println("----------------\n");
			for(String pid:sortedPageRAND.keySet())
				System.out.println(pid+" "+sortedPageRAND.get(pid));
		}
		else{
			int cutoff = DifficultySortingCorpus.DIFFICULTY_PERC*sortedPageRAND.size()/100;
			int count = 0;
			boolean printGap = true;
			for(String pageid:sortedPageRAND.keySet()){
				if(count<=cutoff)
					System.out.println(pageid+" "+sortedPageRAND.get(pageid));
				else if(count>=(sortedPageRAND.size()-cutoff))
					System.out.println(pageid+" "+sortedPageRAND.get(pageid));
				else{
					if(printGap){
						System.out.println("\n\n");
						printGap = false;
					}
				}
				count++;
			}
		}

	}

}
