package edu.unh.cs.treccar.playground;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.DoubleSummaryStatistics;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.regex.Pattern;

import cc.mallet.cluster.Clustering;
import cc.mallet.pipe.CharSequence2TokenSequence;
import cc.mallet.pipe.FeatureSequence2FeatureVector;
import cc.mallet.pipe.Input2CharSequence;
import cc.mallet.pipe.Pipe;
import cc.mallet.pipe.SerialPipes;
import cc.mallet.pipe.TokenSequence2FeatureSequence;
import cc.mallet.pipe.TokenSequenceLowercase;
import cc.mallet.pipe.TokenSequenceRemoveStopwords;
import cc.mallet.topics.TopicAssignment;
import cc.mallet.topics.TopicInferencer;
import cc.mallet.types.Instance;
import cc.mallet.types.InstanceList;
import cc.mallet.types.Metric;
import cc.mallet.types.NormalizedDotProductMetric;
import cc.mallet.types.SparseVector;
import co.nstant.in.cbor.CborException;
import edu.unh.cs.treccar.Data;
import edu.unh.cs.treccar.playground.cluster.CustomKMeans;
import edu.unh.cs.treccar.playground.topics.CustomLDA;
import edu.unh.cs.treccar.playground.topics.UnigramTopicInferencer;
import edu.unh.cs.treccar.playground.topics.UnigramTopicModel;
import edu.unh.cs.treccar.read_data.DeserializeData;
import edu.unh.cs.treccar.read_data.DeserializeData.RuntimeCborException;

public class SingleRun {
	int k, numIter, model, tw;
	double alphaSum, betaSum;
	String paraPath, outlinePath, gtPath, articleQrelsPath, outputPath;
	private HashMap<String, ArrayList<String>> gtSecParaMap = null;
	private HashMap<String, ArrayList<String>> articleParaMap = null;
	private ArrayList<Data.Paragraph> paraList = null;
	private ArrayList<Data.Page> pageList = null;
	public HashMap<String, Double> pageRAND;
	public HashMap<String, Double> pagePurity;
	
	
	public SingleRun(String[] args){
		if(args.length!=11){
			System.out.println("Wrong no. of arguments passed!");
		} else{
			this.k = new Integer(args[0]);
			this.numIter = new Integer(args[1]);
			this.model = new Integer(args[2]);
			this.tw = new Integer(args[3]);
			this.alphaSum = new Double(args[4]);
			this.betaSum = new Double(args[5]);
			this.paraPath = args[6];
			this.outlinePath = args[7];
			this.gtPath = args[8];
			this.articleQrelsPath = args[9];
			this.outputPath = args[10];
			this.paraList = getParaListFromPath(this.paraPath);
			this.pageList = getPageListFromPath(this.outlinePath);
			this.gtSecParaMap = getGTMapFromPath(this.gtPath);
			this.articleParaMap = getArticleMapFromPath(this.articleQrelsPath);
			System.out.println("SingleRun initialized with k="+this.k+
					",numIter="+this.numIter+
					",model="+this.model+
					",tw="+this.tw+
					",alphaSum="+this.alphaSum+
					",betaSum="+this.betaSum);
		}
	}
	public void runExperiment(){
		HashMap<String, ResultForPage> resultPerPageID = new HashMap<String, ResultForPage>();
		int[] luckyPageids = {5, 25, 35};
		int pageid = -1;
		for(Data.Page page:this.pageList){
			/*
			pageid++;
			if(pageid!=luckyPageids[0] && pageid!=luckyPageids[1] && pageid!=luckyPageids[2])
				continue;
			*/
			String pageID = page.getPageId();
			System.out.println("PAGE ID: "+pageID);
			ArrayList<String> paraIDs = this.articleParaMap.get(pageID);
			ArrayList<Data.Paragraph> paraObjects = getParasFromIDs(paraIDs);
			InstanceList paraIList = convertParasToIList(paraObjects);
			ArrayList<String> queryIDs = getqueryIDs(page);
			InstanceList qIList = convertQueriesToIList(queryIDs);
			ResultForPage resultPage = getAssignment(paraIList, qIList, page);
			resultPerPageID.put(pageID, resultPage);
			
			//-- Error Checking and book keeping codes --//
			System.out.println("Size of paraIList = "+paraIList.size());
			for(String q:resultPage.getQueryParaAssignment().keySet()){
				ArrayList<String> paras = resultPage.getQueryParaAssignment().get(q);
				System.out.println(q+" ("+paras.size()+" paras mapped) --> "+paras.toString());
			}
			System.out.println("----------------------------------------------");
			
			ArrayList<String> checkParas = new ArrayList<String>();
			for(int topic=0; topic<resultPage.getParaClusters().size(); topic++){
				ArrayList<String> paras = resultPage.getParaClusters().get(topic);
				System.out.println("Cluster"+topic+" size = "+paras.size()+": "+paras.toString());
				checkParas.addAll(paras);
			}
			Set<String> checkParaSet = new HashSet<String>(checkParas);
			if(checkParaSet.size()<checkParas.size() && this.model!=98){
				try {
					throw new Exception("Duplicate para assigned in "+page.getPageId());
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
					System.exit(1);
				}
			}
			//-- ############################################# --//
		}
		if(RunExperiment.SAVE_RESULT){
			try{
				String runid;
				if(this.model==3)
					runid = "run"+this.k+this.numIter+this.model+this.betaSum+RunExperiment.SMOOTHED_UMM;
				else
					runid = "run"+this.k+this.numIter+this.model+this.tw+this.alphaSum+this.betaSum;
				FileWriter fw = new FileWriter(this.outputPath+"/"+RunExperiment.TRECEVAL_ASSIGN_FILENAME, true);
				for(String p:resultPerPageID.keySet()){
					HashMap<String, ArrayList<String>> currQParaAssign = resultPerPageID.get(p).getQueryParaAssignment();
					for(String q:currQParaAssign.keySet()){
						ArrayList<String> candParaIDs = currQParaAssign.get(q);
						for(String para:candParaIDs)
							fw.write(q+" 0 "+para+" 0 1 "+runid+"\n");
					}
				}
				fw.close();
			} catch(IOException e){
				e.printStackTrace();
			}
		}
		// Measure performance from resultPerPageID and store result
		MeasureExperiment me = new MeasureExperiment(resultPerPageID, this.gtSecParaMap);
		pageRAND = me.calculateRANDPerPage();
		System.out.println("RAND scores\n-----------\n");
		for(String p:pageRAND.keySet()){
			System.out.println(p+" "+pageRAND.get(p));
		}
		pagePurity = me.calculatePurityPerPage();
		System.out.println("Purity scores\n-------------\n");
		for(String p:pagePurity.keySet()){
			System.out.println(p+" "+pagePurity.get(p));
		}
		double sumRAND = 0,sumSqrDev = 0, meanRAND, stdDevRAND, stderrRAND;
		for(Double randVal:pageRAND.values()){
			// We may have some null rand values ( = -99) for which rand index could not be computed
			// (e.g. nom=0.0, denom=0.0), we simply dont add them in our final calc
			if(randVal>-1.01)
				sumRAND+=randVal;
		}
		meanRAND = sumRAND/pageRAND.size();
		for(Double randVal:pageRAND.values())
			if(randVal>-1.01)
				sumSqrDev+=Math.pow(randVal - meanRAND, 2);
		stdDevRAND = Math.sqrt((sumSqrDev/pageRAND.size()));
		stderrRAND = stdDevRAND/Math.sqrt(pageRAND.size());
		
		double sumPurity = 0,sumSqrDevP = 0, meanPurity, stdDevP, stderrPurity;
		for(Double pVal:pagePurity.values()){
			if(pVal!=null)
				sumPurity+=pVal;
		}
		meanPurity = sumPurity/pagePurity.size();
		for(Double pVal:pagePurity.values())
			if(pVal!=null)
				sumSqrDevP+=Math.pow(pVal - meanPurity, 2);
		stdDevP = Math.sqrt((sumSqrDevP/pagePurity.size()));
		stderrPurity = stdDevP/Math.sqrt(pagePurity.size());
		System.out.println("Mean RAND = "+meanRAND+", mean Purity = "+meanPurity);
		if(RunExperiment.SAVE_RESULT){
			try {
				FileWriter fw = new FileWriter(this.outputPath+"/"+RunExperiment.CLUSTERING_MEASURE_FILENAME, true);
				fw.write(this.k+" "+this.numIter+" "+this.model+" "+this.tw+" "+this.alphaSum+" "+this.betaSum+" "+meanRAND+" "+stderrRAND+" "+meanPurity+" "+stderrPurity+"\n");
				fw.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	public void runExperimentWholeCorpus(){
		HashMap<String, ResultForPage> resultForCorpus = new HashMap<String, ResultForPage>();
		List<Data.PageSkeleton> dummySkeleton = pageList.get(0).getSkeleton();
		Data.Page corpus = new Data.Page("Whole Corpus", "Whole%20Corpus", dummySkeleton);
		HashMap<String, ArrayList<String>> corpusGT = this.articleParaMap;
		ArrayList<String> queryids = new ArrayList<String>();
		for(Data.Page p:this.pageList)
			queryids.add(p.getPageId());
		InstanceList paraIList = convertParasToIList(this.paraList);
		InstanceList qIList = convertQueriesToIList(queryids);
		resultForCorpus.put(corpus.getPageId(), getAssignment(paraIList, qIList, corpus));
		HashMap<String, ArrayList<String>> currQParaAssign = resultForCorpus.get(corpus.getPageId()).getQueryParaAssignment();
		
		// Measure performance from resultPerPageID and store result
		MeasureExperiment me = new MeasureExperiment(resultForCorpus, corpusGT);
		double randCorpus = (Double)me.calculateRANDPerPage().values().toArray()[0];
		double purityCorpus = (Double)me.calculatePurityPerPage().values().toArray()[0];
		System.out.println("Corpus RAND score = "+randCorpus+", Purity score = "+purityCorpus);
		if(RunExperiment.SAVE_RESULT){
			try {
				FileWriter fw = new FileWriter(this.outputPath+"/"+RunExperiment.CLUSTERING_MEASURE_FILENAME, true);
				fw.write(this.k+" "+this.numIter+" "+this.model+" "+this.tw+" "+this.alphaSum+" "+this.betaSum+" "+randCorpus+" 0 "+purityCorpus+" 0\n");
				fw.close();
				
				String runid = "runWholeCorpus"+this.k+this.numIter+this.model+this.tw+this.alphaSum+this.betaSum;
				FileWriter fw2 = new FileWriter(this.outputPath+"/"+RunExperiment.TRECEVAL_ASSIGN_FILENAME, true);	
				for(String q:currQParaAssign.keySet()){
					ArrayList<String> candParaIDs = currQParaAssign.get(q);
					for(String para:candParaIDs)
						fw2.write(q.split("/")[1]+" 0 "+para+" 0 1 "+runid+"\n");
				}
				fw2.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	private ResultForPage getAssignment(InstanceList paraIList, InstanceList queryIList, Data.Page currPage){
		//HashMap<String, ArrayList<String>> assignment = new HashMap<String, ArrayList<String>>();
		ResultForPage result = new ResultForPage();
		switch(this.model){
		case 1:
			//topic
			int numTopics;
			if(this.k==0){
				/*
				if(queryIList.size()>paraIList.size()){
					numTopics = paraIList.size();
					System.out.println("K was "+queryIList.size()+", but we've reduced it to be paralist size "+numTopics);
				} else{
					numTopics = queryIList.size();
					System.out.println("K "+numTopics);
				}
				*/
				numTopics = queryIList.size();
				System.out.println("K "+numTopics);
			} else{
				/*
				if(this.k>paraIList.size()){
					numTopics = paraIList.size();
					System.out.println("K was "+this.k+", but we've reduced it to be paralist size "+numTopics);
				} else{
					numTopics = this.k;
					System.out.println("K "+numTopics);
				}
				*/
				if(this.k>=paraIList.size()){
					numTopics = paraIList.size();
					System.out.println("K "+numTopics);
				} else{
					numTopics = this.k;
					System.out.println("K "+numTopics);
				}
			}
			CustomLDA lda = new CustomLDA(numTopics, this.alphaSum, this.betaSum);
			try {
				lda.addInstances(paraIList);
				lda.sample(this.numIter);
				//lda.addInstancesUnigram(paraIList);
				//lda.sampleUMM(this.numIter);
				//System.out.println(lda.topWords(10));
				result = assignUsingLDA(paraIList, queryIList, lda, currPage);
				
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			break;
		case 2:
			//cluster
			int numClusters;
			if(this.k==0){
				if(queryIList.size()<paraIList.size())
					numClusters = queryIList.size();
				else
					numClusters = paraIList.size();
				System.out.println("K "+numClusters);
			} else{
				if(this.k>=paraIList.size()){
					numClusters = paraIList.size();
					System.out.println("K "+numClusters);
				} else{
					numClusters = this.k;
					System.out.println("K "+numClusters);
				}
			}
			CustomKMeans kmeans = new CustomKMeans(paraIList.getPipe(), numClusters, new NormalizedDotProductMetric(), CustomKMeans.EMPTY_DROP);
			//ArrayList<ArrayList<String>> clusters = formatClusterData(kmeans.cluster(paraIList));
			result = assignUsingKMeans(paraIList, queryIList, kmeans.cluster(paraIList), kmeans, currPage);
			break;
		case 3:
			//unigram topic
			int numTopicsUMM;
			if(this.k==0){
				/*
				if(queryIList.size()>paraIList.size()){
					numTopics = paraIList.size();
					System.out.println("K was "+queryIList.size()+", but we've reduced it to be paralist size "+numTopics);
				} else{
					numTopics = queryIList.size();
					System.out.println("K "+numTopics);
				}
				*/
				numTopicsUMM = queryIList.size();
				System.out.println("K "+numTopicsUMM);
			} else{
				/*
				if(this.k>paraIList.size()){
					numTopics = paraIList.size();
					System.out.println("K was "+this.k+", but we've reduced it to be paralist size "+numTopics);
				} else{
					numTopics = this.k;
					System.out.println("K "+numTopics);
				}
				*/
				if(this.k>=paraIList.size()){
					numTopicsUMM = paraIList.size();
					System.out.println("K "+numTopicsUMM);
				} else{
					numTopicsUMM = this.k;
					System.out.println("K "+numTopicsUMM);
				}
			}
			UnigramTopicModel ummlda = new UnigramTopicModel(numTopicsUMM, this.alphaSum, this.betaSum);
			try {
				ummlda.addInstances(paraIList);
				ummlda.sample(this.numIter);
				result = assignUsingUMM(paraIList, queryIList, ummlda, currPage);
				
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			break;
		case 98:
			//all correct
			result = assignUsingAllCorrect(paraIList, queryIList, currPage);
			break;
		case 99:
			//random
			result = assignUsingRandom(paraIList, queryIList, currPage);
			break;
		default:
			System.out.println("Wrong model option!");
		}
		return result;
	}
	private HashMap<String, ArrayList<String>> convInsAssignToIDAssign(HashMap<Instance, ArrayList<Instance>> insAssign, Data.Page page){
		HashMap<String, ArrayList<String>> idAssign = new HashMap<String, ArrayList<String>>();
		String pageID = page.getPageId();
		for(Instance qIns:insAssign.keySet()){
			String queryID = qIns.getName().toString();
			ArrayList<Instance> paraInsList = insAssign.get(qIns);
			ArrayList<String> paraIDs = new ArrayList<String>();
			if(queryID.equals(pageID))
				queryID = pageID;
			else
				queryID = pageID+"/"+queryID;
			for(Instance paraIns:paraInsList)
				paraIDs.add(paraIns.getName().toString());
			idAssign.put(queryID, paraIDs);
		}
		return idAssign;
	}
	private ResultForPage assignUsingAllCorrect(InstanceList paraIList, InstanceList queryIList, Data.Page page){
		HashMap<String, ArrayList<String>> assign = new HashMap<String, ArrayList<String>>();
		ArrayList<ArrayList<String>> paraClusters = new ArrayList<ArrayList<String>>();
		ResultForPage r = new ResultForPage();
		ArrayList<String> queryIDs = new ArrayList<String>();
		ArrayList<String> paraIDs = new ArrayList<String>();
		for(Instance qIns:queryIList){
			String pageid = page.getPageId();
			String q = qIns.getName().toString();
			if(q.equals(pageid))
				q = pageid;
			else
				q = pageid+"/"+q;
			assign.put(q, new ArrayList<String>());
			queryIDs.add(q);
		}
		for(String q:queryIDs){
			paraIDs = this.gtSecParaMap.get(q);
			if(paraIDs==null){
				System.out.println("No entry for query "+q+" in gt!");
				continue;
			}
			assign.get(q).addAll(paraIDs);
		}
		for(String qry:assign.keySet())
			paraClusters.add(assign.get(qry));
		r.setQueryParaAssignment(assign);
		r.setParaClusters(paraClusters);
		return r;
	}
	private ResultForPage assignUsingRandom(InstanceList paraIList, InstanceList queryIList, Data.Page page){
		Random rand = new Random();
		HashMap<String, ArrayList<String>> assign = new HashMap<String, ArrayList<String>>();
		ArrayList<ArrayList<String>> paraClusters = new ArrayList<ArrayList<String>>();
		ResultForPage r = new ResultForPage();
		ArrayList<String> queryIDs = new ArrayList<String>();
		ArrayList<String> paraIDs = new ArrayList<String>();
		for(Instance qIns:queryIList){
			String pageid = page.getPageId();
			String q = qIns.getName().toString();
			if(q.equals(pageid))
				q = pageid;
			else
				q = pageid+"/"+q;
			assign.put(q, new ArrayList<String>());
			queryIDs.add(q);
		}
		for(Instance pIns:paraIList)
			paraIDs.add(pIns.getName().toString());
		int qindex=0, pindex;
		// To guarantee that every sec/que got some para //
		while(!paraIDs.isEmpty() && qindex<queryIDs.size()){
			if(paraIDs.size()>1)
				pindex = rand.nextInt(paraIDs.size()-1);
			else
				pindex = 0;
			assign.get(queryIDs.get(qindex)).add(paraIDs.get(pindex));
			paraIDs.remove(pindex);
			qindex++;
		}
		// ---------------------------------------------- //
		while(!paraIDs.isEmpty()){
			qindex = rand.nextInt(queryIDs.size()-1);
			if(paraIDs.size()>1)
				pindex = rand.nextInt(paraIDs.size()-1);
			else
				pindex = 0;
			
			assign.get(queryIDs.get(qindex)).add(paraIDs.get(pindex));
			paraIDs.remove(pindex);
		}
		for(String qry:assign.keySet())
			paraClusters.add(assign.get(qry));
		r.setQueryParaAssignment(assign);
		r.setParaClusters(paraClusters);
		return r;
	}
	private ResultForPage assignUsingLDA(InstanceList paraIList, InstanceList queryIList, CustomLDA lda, Data.Page page) throws Exception{
		HashMap<Instance, ArrayList<Instance>> assign = new HashMap<Instance, ArrayList<Instance>>();
		for(Instance qIns:queryIList)
			assign.put(qIns, new ArrayList<Instance>());
		
		TopicInferencer inf = lda.getInferencer();
		int numIterForInf = 30;
		int thinningForInf = 1;
		int burninForInf = 5;
		double[][] queryTopicProbMatrix = new double[queryIList.size()][lda.getNumTopics()];
		double[][] paraTopicProbMatrix = new double[paraIList.size()][lda.getNumTopics()];
		for(int i=0; i<queryIList.size(); i++){
			Instance queryIns = queryIList.get(i);
			//System.out.println("Page ID: "+page.getPageId());
			double[] currQueryTopicProbDist = inf.getSampledDistribution(queryIns, numIterForInf, thinningForInf, burninForInf);
			queryTopicProbMatrix[i] = currQueryTopicProbDist;
		}
		if(paraIList.size()!=lda.getData().size())
			throw new Exception("paralist size and lda topic assignment size dont match!");
		for(int j=0; j<paraIList.size(); j++){
			// Following if block ensures that we are picking the correct topic dist for current para instance
			if(paraIList.get(j).getName()!=lda.getData().get(j).instance.getName())
				throw new Exception("paraIList indices are not following the same order as lda instances");
			paraTopicProbMatrix[j] = lda.getTopicProbabilities(j);
		}
		boolean[] isParaAssigned = new boolean[paraIList.size()];
		for(int pos=0; pos<isParaAssigned.length; pos++)
			isParaAssigned[pos] = false;
		for(int m=0; m<queryIList.size(); m++){
			Instance queryIns = queryIList.get(m);
			int bestParaInsIndex = 0;
			double minKLDiv = 99999.0;
			for(int n=0; n<paraIList.size(); n++){
				if(!isParaAssigned[n]){
					double currKLDiv = getKLdiv(queryTopicProbMatrix[m], paraTopicProbMatrix[n]);
					if(currKLDiv<minKLDiv){
						minKLDiv = currKLDiv;
						bestParaInsIndex = n;
					}
				}
			}
			assign.get(queryIns).add(paraIList.get(bestParaInsIndex));
			isParaAssigned[bestParaInsIndex] = true;
		}
		for(int p=0; p<paraIList.size(); p++){
			if(!isParaAssigned[p]){
				Instance bestQueryIns = queryIList.get(0);
				double minKLDiv = 99999.0;
				for(int q=0; q<queryIList.size(); q++){
					double currKLDiv = getKLdiv(queryTopicProbMatrix[q], paraTopicProbMatrix[p]);
					if(currKLDiv<minKLDiv){
						minKLDiv = currKLDiv;
						bestQueryIns = queryIList.get(q);
					}
				}
				assign.get(bestQueryIns).add(paraIList.get(p));
			}
		}
		HashMap<String, ArrayList<String>> assignment = convInsAssignToIDAssign(assign, page);
		ResultForPage r = new ResultForPage();
		r.setQueryParaAssignment(assignment);
		r.setParaClusters(getParasClustersFromMatrix(paraTopicProbMatrix, paraIList, lda.getNumTopics()));
		return r;
	}
	private ResultForPage assignUsingUMM(InstanceList paraIList, InstanceList queryIList, UnigramTopicModel umm, Data.Page page) throws Exception{
		HashMap<Instance, ArrayList<Instance>> assign = new HashMap<Instance, ArrayList<Instance>>();
		for(Instance qIns:queryIList)
			assign.put(qIns, new ArrayList<Instance>());
		
		UnigramTopicInferencer inf = umm.getInferencer();
		int numIterForInf = 30;
		int thinningForInf = 1;
		int burninForInf = 5;
		int[] queryTopics = new int[queryIList.size()];
		int[] paraTopics = new int[paraIList.size()];
		for(int i=0; i<queryIList.size(); i++){
			Instance queryIns = queryIList.get(i);
			//System.out.println("Page ID: "+page.getPageId());
			int currQueryTopic = inf.inferInstanceTopic(queryIns, numIterForInf, thinningForInf, burninForInf);
			queryTopics[i] = currQueryTopic;
		}
		if(paraIList.size()!=umm.getData().size())
			throw new Exception("paralist size and lda topic assignment size dont match!");
		for(int j=0; j<paraIList.size(); j++){
			// Following if block ensures that we are picking the correct topic dist for current para instance
			if(paraIList.get(j).getName()!=umm.getData().get(j).instance.getName())
				throw new Exception("paraIList indices are not following the same order as ummlda instances");
			paraTopics[j] = umm.getTopicOfInstance(j);
		}
		boolean[] isParaAssigned = new boolean[paraIList.size()];
		for(int pos=0; pos<isParaAssigned.length; pos++)
			isParaAssigned[pos] = false;
		/*for(int m=0; m<queryIList.size(); m++){
			Instance queryIns = queryIList.get(m);
			int bestParaInsIndex = 0;
			double minKLDiv = 99999.0;
			for(int n=0; n<paraIList.size(); n++){
				if(!isParaAssigned[n]){
					double currKLDiv = getKLdiv(queryTopicProbMatrix[m], paraTopicProbMatrix[n]);
					if(currKLDiv<minKLDiv){
						minKLDiv = currKLDiv;
						bestParaInsIndex = n;
					}
				}
			}
			assign.get(queryIns).add(paraIList.get(bestParaInsIndex));
			isParaAssigned[bestParaInsIndex] = true;
		}
		for(int p=0; p<paraIList.size(); p++){
			if(!isParaAssigned[p]){
				Instance bestQueryIns = queryIList.get(0);
				double minKLDiv = 99999.0;
				for(int q=0; q<queryIList.size(); q++){
					double currKLDiv = getKLdiv(queryTopicProbMatrix[q], paraTopicProbMatrix[p]);
					if(currKLDiv<minKLDiv){
						minKLDiv = currKLDiv;
						bestQueryIns = queryIList.get(q);
					}
				}
				assign.get(bestQueryIns).add(paraIList.get(p));
			}
		}*/
		for(int q=0; q<queryTopics.length; q++){
			Instance currQIns = queryIList.get(q);
			for(int p=0; p<paraTopics.length; p++){
				if(!isParaAssigned[p]){
					if(queryTopics[q]==paraTopics[p]){
						assign.get(currQIns).add(paraIList.get(p));
						isParaAssigned[p] = true;
					}
				}
			}
		}
		HashMap<String, ArrayList<String>> assignment = convInsAssignToIDAssign(assign, page);
		ResultForPage r = new ResultForPage();
		r.setQueryParaAssignment(assignment);
		r.setParaClusters(getParasClustersFromParatopics(paraTopics, paraIList, umm.getNumTopics()));
		return r;
	}
	private ResultForPage assignUsingKMeans(InstanceList paraIList, InstanceList queryIList, 
			Clustering clusters, CustomKMeans kmeans, Data.Page page){
		ResultForPage r = new ResultForPage();
		HashMap<Instance, ArrayList<Instance>> insAssign = new HashMap<Instance, ArrayList<Instance>>();
		double[][] queryClusterDistMat = getQueryClusterDistances(queryIList, clusters, kmeans, new NormalizedDotProductMetric());
		boolean[] isClusterAssigned = new boolean[clusters.getNumClusters()];
		for(int i=0; i<isClusterAssigned.length; i++)
			isClusterAssigned[i] = false;
		while(hasFalse(isClusterAssigned)){
			int qIndex=0, cIndex=0;
			double minVal=9999999.0, val;
			for(int i=0; i<queryClusterDistMat.length; i++){
				for(int j=0; j<queryClusterDistMat[0].length; j++){
					val = queryClusterDistMat[i][j];
					if(val>0 && val<minVal){
						minVal = val;
						qIndex = i;
						cIndex = j;
					}
				}
			}
			insAssign.put(queryIList.get(qIndex), clusters.getCluster(clusters.getLabel(cIndex)));
			for(int i=0; i<queryClusterDistMat.length; i++)
				queryClusterDistMat[i][cIndex] = -1;
			isClusterAssigned[cIndex] = true;
		}
		HashMap<String, ArrayList<String>> assignment = convInsAssignToIDAssign(insAssign, page);
		r.setQueryParaAssignment(assignment);
		r.setParaClusters(formatClusterData(clusters));
		return r;
	}
	private boolean hasFalse(boolean[] boolArray){
		for(int i=0; i<boolArray.length; i++){
			if(!boolArray[i])
				return true;
		}
		return false;
	}
	private double[][] getQueryClusterDistances(InstanceList querylist, Clustering clusters, CustomKMeans kmeans, Metric metric){
		double[][] distanceMatrix = new double[querylist.size()][clusters.getNumClusters()];
		ArrayList<SparseVector> clusterMeans = kmeans.getClusterMeans();
		Instance query;
		for(int q=0; q<distanceMatrix.length; q++){
			for(int c=0; c<distanceMatrix[0].length; c++)
				distanceMatrix[q][c] = 
				metric.distance(clusterMeans.get(c), (SparseVector) querylist.get(q).getData());
		}
		return distanceMatrix;
	}
	private ArrayList<ArrayList<String>> formatClusterData(Clustering rawClusterData){
		ArrayList<ArrayList<String>> finalClusterData = new ArrayList<ArrayList<String>>();
		for(InstanceList iList : rawClusterData.getClusters()){
			ArrayList<String> paraIDList = new ArrayList<String>();
			for(Instance i : iList){
				paraIDList.add(i.getName().toString());
			}
			finalClusterData.add(paraIDList);
		}
		return finalClusterData;
	}
	private ArrayList<ArrayList<String>> getParasClustersFromMatrix(double[][] paraTopicMatrix, InstanceList pIList, int numTopics) throws Exception{
		ArrayList<ArrayList<String>> clusters = new ArrayList<ArrayList<String>>();
		for(int topicPos=0; topicPos<numTopics; topicPos++)
			clusters.add(new ArrayList<String>());
		for(int i=0; i<paraTopicMatrix.length; i++){
			Instance p = pIList.get(i);
			double[] currParaTopicDist = paraTopicMatrix[i];
			if(currParaTopicDist.length!=numTopics)
				throw new Exception("Current topic dist array size does not match k!");
			double maxProb = 0;
			int maxTopicIndex = 0;
			for(int j=0; j<currParaTopicDist.length; j++){
				if(currParaTopicDist[j]>maxProb){
					maxTopicIndex = j;
					maxProb = currParaTopicDist[j];
				}
			}
			clusters.get(maxTopicIndex).add(p.getName().toString());
		}
		return clusters;
	}
	private ArrayList<ArrayList<String>> getParasClustersFromParatopics(int[] paraTopics, InstanceList pIList, int numTopics){
		ArrayList<ArrayList<String>> clusters = new ArrayList<ArrayList<String>>();
		for(int topicPos=0; topicPos<numTopics; topicPos++)
			clusters.add(new ArrayList<String>());
		for(int paraIndex=0; paraIndex<paraTopics.length; paraIndex++){
			Instance p = pIList.get(paraIndex);
			clusters.get(paraTopics[paraIndex]).add(p.getName().toString());
		}
		return clusters;
	}
	private double getKLdiv(double[] p, double[] q){
		double result = 0;
		for(int i=0; i<p.length; i++){
			if(q[i]<0.0000001 || p[i]<0.0000001){
				continue;
			}
			result+=p[i]*Math.log(p[i]/q[i]);
		}
		return result;
	}
	private InstanceList convertQueriesToIList(ArrayList<String> queries){
		InstanceList qIList = new InstanceList(SingleRun.buildPipe(this.model));
		for(String q:queries){
			Instance secIns = new Instance(q.replaceAll("%20|/", " "), null, q, q);
			qIList.addThruPipe(secIns);
		}
		return qIList;
	}
	private InstanceList convertParasToIList(ArrayList<Data.Paragraph> paraObjs){
		InstanceList iListPara = new InstanceList(SingleRun.buildPipe(this.model));
		for(Data.Paragraph paraObj:paraObjs){
			Instance paraIns = new Instance(paraObj.getTextOnly(), null, paraObj.getParaId(), paraObj.getTextOnly());
			iListPara.addThruPipe(paraIns);
		}
		return iListPara;
	}
	private ArrayList<String> getqueryIDs(Data.Page page){
		ArrayList<String> queries = new ArrayList<String>();
		queries.add(page.getPageId());
		for(Data.Section sec:page.getChildSections())
			queries.add(sec.getHeadingId());
		return queries;
	}
	private ArrayList<Data.Paragraph> getParasFromIDs(ArrayList<String> paraIDs){
		ArrayList<Data.Paragraph> paras = new ArrayList<Data.Paragraph>();
		for(String pID:paraIDs){
			for(Data.Paragraph p:this.paraList){
				if(pID.equals(p.getParaId()))
					paras.add(p);
			}
		}
		return paras;
	}
	private HashMap<String, ArrayList<String>> getArticleMapFromPath(String path){
		HashMap<String, ArrayList<String>> articleMap = new HashMap<String, ArrayList<String>>();
		BufferedReader br;
		try{
			br = new BufferedReader(new FileReader(path));
			String line;
			String[] lineData = new String[4];
			while((line = br.readLine()) != null){
				lineData = line.split(" ");
				if(articleMap.containsKey(lineData[0])){
					articleMap.get(lineData[0]).add(lineData[2]);
				} else{
					ArrayList<String> paraList = new ArrayList<String>();
					paraList.add(lineData[2]);
					articleMap.put(lineData[0], paraList);
				}	
			}
			br.close();
		} catch (Exception e){
			e.printStackTrace();
		}
		return articleMap;
	}
	private ArrayList<Data.Paragraph> getParaListFromPath(String path){
		ArrayList<Data.Paragraph> paraList = new ArrayList<Data.Paragraph>();
		try {
			FileInputStream fis = new FileInputStream(new File(path));
			for(Data.Paragraph para:DeserializeData.iterableParagraphs(fis))
				paraList.add(para);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (CborException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return paraList;
	}
	private ArrayList<Data.Page> getPageListFromPath(String path){
		ArrayList<Data.Page> pageList = new ArrayList<Data.Page>();
		try {
			FileInputStream fis = new FileInputStream(new File(path));
			for(Data.Page page: DeserializeData.iterableAnnotations(fis))
				pageList.add(page);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (RuntimeCborException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return pageList;
	}
	private HashMap<String, ArrayList<String>> getGTMapFromPath(String path){
		HashMap<String, ArrayList<String>> gtMap = new HashMap<String, ArrayList<String>>();
		BufferedReader br;
		try{
			br = new BufferedReader(new FileReader(path));
			String line;
			String[] lineData = new String[4];
			while((line = br.readLine()) != null){
				lineData = line.split(" ");
				if(gtMap.containsKey(lineData[0])){
					gtMap.get(lineData[0]).add(lineData[2]);
				} else{
					ArrayList<String> paraList = new ArrayList<String>();
					paraList.add(lineData[2]);
					gtMap.put(lineData[0], paraList);
				}	
			}
			br.close();
		} catch (Exception e){
			e.printStackTrace();
		}
		return gtMap;
	}
	public static Pipe buildPipe(int model){
		if(model==1 || model==3)
			return buildPipeForLDA();
		else
			return buildPipeDefault();
			
	}
	public static Pipe buildPipeForLDA(){
		ArrayList pipeList = new ArrayList();

        // Read data from File objects
        pipeList.add(new Input2CharSequence("UTF-8"));

        // Regular expression for what constitutes a token.
        //  This pattern includes Unicode letters, Unicode numbers, 
        //   and the underscore character. Alternatives:
        //    "\\S+"   (anything not whitespace)
        //    "\\w+"    ( A-Z, a-z, 0-9, _ )
        //    "[\\p{L}\\p{N}_]+|[\\p{P}]+"   (a group of only letters and numbers OR
        //                                    a group of only punctuation marks)
        Pattern tokenPattern = Pattern.compile("[\\p{L}\\p{N}_]+");

        // Tokenize raw strings
        pipeList.add(new CharSequence2TokenSequence(tokenPattern));

        // Normalize all tokens to all lowercase
        pipeList.add(new TokenSequenceLowercase());

        // Remove stopwords from a standard English stoplist.
        //  options: [case sensitive] [mark deletions]
        pipeList.add(new TokenSequenceRemoveStopwords(false, false));

        // Rather than storing tokens as strings, convert 
        //  them to integers by looking them up in an alphabet.
        pipeList.add(new TokenSequence2FeatureSequence());

        // Do the same thing for the "target" field: 
        //  convert a class label string to a Label object,
        //  which has an index in a Label alphabet.
        //pipeList.add(new Target2Label());

        // Now convert the sequence of features to a sparse vector,
        //  mapping feature IDs to counts.
        //pipeList.add(new FeatureSequence2FeatureVector());
        
        // Print out the features and the label
        //pipeList.add(new PrintInputAndTarget());
        
        return (new SerialPipes(pipeList));
	}
	public static Pipe buildPipeDefault(){
		ArrayList pipeList = new ArrayList();

        // Read data from File objects
        pipeList.add(new Input2CharSequence("UTF-8"));

        // Regular expression for what constitutes a token.
        //  This pattern includes Unicode letters, Unicode numbers, 
        //   and the underscore character. Alternatives:
        //    "\\S+"   (anything not whitespace)
        //    "\\w+"    ( A-Z, a-z, 0-9, _ )
        //    "[\\p{L}\\p{N}_]+|[\\p{P}]+"   (a group of only letters and numbers OR
        //                                    a group of only punctuation marks)
        Pattern tokenPattern = Pattern.compile("[\\p{L}\\p{N}_]+");

        // Tokenize raw strings
        pipeList.add(new CharSequence2TokenSequence(tokenPattern));

        // Normalize all tokens to all lowercase
        pipeList.add(new TokenSequenceLowercase());

        // Remove stopwords from a standard English stoplist.
        //  options: [case sensitive] [mark deletions]
        pipeList.add(new TokenSequenceRemoveStopwords(false, false));

        // Rather than storing tokens as strings, convert 
        //  them to integers by looking them up in an alphabet.
        pipeList.add(new TokenSequence2FeatureSequence());

        // Do the same thing for the "target" field: 
        //  convert a class label string to a Label object,
        //  which has an index in a Label alphabet.
        //pipeList.add(new Target2Label());

        // Now convert the sequence of features to a sparse vector,
        //  mapping feature IDs to counts.
        pipeList.add(new FeatureSequence2FeatureVector());
        
        // Print out the features and the label
        //pipeList.add(new PrintInputAndTarget());
        
        return (new SerialPipes(pipeList));
	}
}
