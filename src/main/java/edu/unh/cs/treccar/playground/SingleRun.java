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
import java.util.Iterator;
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
import cc.mallet.types.FeatureSequence;
import cc.mallet.types.FeatureVector;
import cc.mallet.types.Instance;
import cc.mallet.types.InstanceList;
import cc.mallet.types.Metric;
import cc.mallet.types.NormalizedDotProductMetric;
import cc.mallet.types.SparseVector;
import cc.mallet.util.VectorStats;
import co.nstant.in.cbor.CborException;
import edu.unh.cs.treccar.Data;
import edu.unh.cs.treccar.playground.cluster.CustomKMeans;
import edu.unh.cs.treccar.playground.topics.CustomLDA;
import edu.unh.cs.treccar.playground.topics.UnigramTopicInferencer;
import edu.unh.cs.treccar.playground.topics.UnigramTopicModel;
import edu.unh.cs.treccar.read_data.DeserializeData;
import edu.unh.cs.treccar.read_data.DeserializeData.RuntimeCborException;

public class SingleRun {
	int k, numIter, model, tw, countNan=0, removed=0;
	ArrayList<String> nanIns = new ArrayList<String>();
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
	private void removeEmptyInstances(InstanceList iList){
		Iterator<Instance> iListIter = iList.iterator();
		if(iList.getDataClass()==FeatureSequence.class){
			FeatureSequence fet;
			
			while(iListIter.hasNext()){
				fet = (FeatureSequence) iListIter.next().getData();
				if(fet.getLength()==0){
					iListIter.remove();
					removed++;
				}
			}
		} else{
			SparseVector sv;
			while(iListIter.hasNext()){
				sv = (SparseVector) iListIter.next().getData();
				if(sv.numLocations()==0){
					iListIter.remove();
					removed++;
				}
			}
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
			if(pageID.startsWith("Moment%20magnitude"))
				System.out.println();
			System.out.println("PAGE ID: "+pageID);
			ArrayList<String> paraIDs = this.articleParaMap.get(pageID);
			ArrayList<Data.Paragraph> paraObjects = getParasFromIDs(paraIDs);
			InstanceList paraIList = convertParasToIList(paraObjects);
			ArrayList<String> queryIDs = getqueryIDs(page);
			InstanceList qIList = convertQueriesToIList(queryIDs);
			removeEmptyInstances(qIList);
			removeEmptyInstances(paraIList);
			ResultForPage resultPage = modelAndAssign(paraIList, qIList, page); // clustering and assignments done here
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
		//System.out.println("NAN Instances: "+this.countNan);
		//for(String nan:this.nanIns)
			//System.out.println(nan);
		System.out.println("No. of removed instances: "+removed);
		if(RunExperiment.SAVE_RESULT){
			try {
				FileWriter fw = new FileWriter(this.outputPath+"/"+RunExperiment.CLUSTERING_MEASURE_FILENAME, true);
				fw.write(this.k+" "+this.numIter+" "+this.model+" "+this.tw+" "+this.alphaSum+" "+this.betaSum+" "+meanRAND+" "+stderrRAND+" "+meanPurity+" "+stderrPurity+"\n");
				fw.close();
				
				String runid;
				if(this.model==3)
					runid = "run"+this.k+this.numIter+this.model+this.betaSum+RunExperiment.SMOOTHED_UMM;
				else
					runid = "run"+this.k+this.numIter+this.model+this.tw+this.alphaSum+this.betaSum;
				FileWriter fw2 = new FileWriter(this.outputPath+"/"+RunExperiment.TRECEVAL_ASSIGN_FILENAME, true);
				for(String p:resultPerPageID.keySet()){
					if(p.equals("Miko"))
						System.out.println("Debug here");
					HashMap<String, ArrayList<String>> currQParaAssign = resultPerPageID.get(p).getQueryParaAssignment();
					HashMap<ArrayList<String>, Double> currQParaRank = resultPerPageID.get(p).getQueryParaRank();
					if(!RunExperiment.ASSIGN_RANK_MODE){
						for(String q:currQParaAssign.keySet()){
							ArrayList<String> candParaIDs = currQParaAssign.get(q);
							for(String para:candParaIDs)
								fw2.write(q+" 0 "+para+" 0 1 "+runid+"\n");
							// (query-id ignore document-id ignore rank-score runid)
						}
					} else{
						for(ArrayList<String> qp:currQParaRank.keySet())
							fw2.write(qp.get(0)+" 0 "+qp.get(1)+" 0 "+currQParaRank.get(qp)+" "+runid+"\n");
					}
				}
				fw2.close();
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
		removeEmptyInstances(qIList);
		resultForCorpus.put(corpus.getPageId(), modelAndAssign(paraIList, qIList, corpus));
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
	private ResultForPage modelAndAssign(InstanceList paraIList, InstanceList queryIList, Data.Page currPage){
		//HashMap<String, ArrayList<String>> assignment = new HashMap<String, ArrayList<String>>();
		ResultForPage result = new ResultForPage();
		switch(this.model){
		case 1:
			//topic
			int numTopics;
			if(this.k==0){
				
				numTopics = queryIList.size();
				System.out.println("K "+numTopics);
			} else{
				
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
				
				numTopicsUMM = queryIList.size();
				System.out.println("K "+numTopicsUMM);
			} else{
				
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
	private HashMap<ArrayList<String>, Double> convRanksToIDRanks(HashMap<ArrayList<Instance>, Double> ranks, Data.Page page){
		HashMap<ArrayList<String>, Double> idRanks = new HashMap<ArrayList<String>, Double>();
		String pageID = page.getPageId();
		String q,p;
		ArrayList<String> rankKey;
		for(ArrayList<Instance> sp:ranks.keySet()){
			q = sp.get(0).getName().toString();
			p = sp.get(1).getName().toString();
			if(!q.equals(pageID))
				q = pageID+"/"+q;
			rankKey = new ArrayList<String>();
			rankKey.add(0, q);
			rankKey.add(1, p);
			idRanks.put(rankKey, ranks.get(sp));
		}
		return idRanks;
	}
	private ResultForPage assignUsingAllCorrect(InstanceList paraIList, InstanceList queryIList, Data.Page page){
		HashMap<String, ArrayList<String>> assign = new HashMap<String, ArrayList<String>>();
		HashMap<ArrayList<Instance>, Double> ranks = new HashMap<ArrayList<Instance>, Double>();
		ArrayList<ArrayList<String>> paraClusters = new ArrayList<ArrayList<String>>();
		ResultForPage r = new ResultForPage();
		ArrayList<String> queryIDs = new ArrayList<String>();
		ArrayList<String> paraIDs = new ArrayList<String>();
		String pageid = page.getPageId();
		for(Instance qIns:queryIList){
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
		ArrayList<Instance> rankKey;
		String qid, pid;
		for(Instance q:queryIList){
			qid = q.getName().toString();
			if(qid.equals(pageid))
				qid = pageid;
			else
				qid = pageid+"/"+qid;
			for(Instance p:paraIList){
				pid = p.getName().toString();
				rankKey = new ArrayList<Instance>();
				rankKey.add(0, q);
				rankKey.add(1, p);
				if(assign.get(qid).contains(pid))
					ranks.put(rankKey, 1.0);
				else
					ranks.put(rankKey, 0.0);
			}
		}
		for(String qry:assign.keySet())
			paraClusters.add(assign.get(qry));
		r.setQueryParaAssignment(assign);
		r.setQueryParaRank(convRanksToIDRanks(ranks, page));
		r.setParaClusters(paraClusters);
		return r;
	}
	private ResultForPage assignUsingRandom(InstanceList paraIList, InstanceList queryIList, Data.Page page){
		Random rand = new Random();
		HashMap<String, ArrayList<String>> assign = new HashMap<String, ArrayList<String>>();
		HashMap<ArrayList<Instance>, Double> ranks = new HashMap<ArrayList<Instance>, Double>();
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
		ArrayList<Instance> rankKey;
		for(Instance q:queryIList){
			for(Instance p:paraIList){
				rankKey = new ArrayList<Instance>();
				rankKey.add(0, q);
				rankKey.add(1, p);
				ranks.put(rankKey, rand.nextDouble());
			}
		}
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
		r.setQueryParaRank(convRanksToIDRanks(ranks, page));
		r.setParaClusters(paraClusters);
		return r;
	}
	private ResultForPage assignUsingLDA(InstanceList paraIList, InstanceList queryIList, CustomLDA lda, Data.Page page) throws Exception{
		HashMap<Instance, ArrayList<Instance>> assign = new HashMap<Instance, ArrayList<Instance>>();
		HashMap<ArrayList<Instance>, Double> ranks = new HashMap<ArrayList<Instance>, Double>();
		FeatureSequence fet;
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
			queryTopicProbMatrix[i] = inf.getSampledDistribution(queryIns, numIterForInf, thinningForInf, burninForInf);
		}
		if(paraIList.size()!=lda.getData().size())
			throw new Exception("paralist size and lda topic assignment size dont match!");
		for(int j=0; j<paraIList.size(); j++){
			// Following if block ensures that we are picking the correct topic dist for current para instance
			if(paraIList.get(j).getName()!=lda.getData().get(j).instance.getName())
				throw new Exception("paraIList indices are not following the same order as lda instances");
			paraTopicProbMatrix[j] = lda.getTopicProbabilities(j);
		}
		matrixAssignment(assign, ranks, queryIList, paraIList, queryTopicProbMatrix, paraTopicProbMatrix);
		HashMap<String, ArrayList<String>> assignment = convInsAssignToIDAssign(assign, page);
		HashMap<ArrayList<String>, Double> idRanks = convRanksToIDRanks(ranks, page);
		ResultForPage r = new ResultForPage();
		r.setQueryParaAssignment(assignment);
		r.setQueryParaRank(idRanks);
		r.setParaClusters(getParasClustersFromMatrix(paraTopicProbMatrix, paraIList, lda.getNumTopics()));
		return r;
	}
	private ResultForPage assignUsingUMM(InstanceList paraIList, InstanceList queryIList, UnigramTopicModel umm, Data.Page page) throws Exception{
		HashMap<Instance, ArrayList<Instance>> assign = new HashMap<Instance, ArrayList<Instance>>();
		HashMap<ArrayList<Instance>, Double> ranks = new HashMap<ArrayList<Instance>, Double>();
		for(Instance qIns:queryIList)
			assign.put(qIns, new ArrayList<Instance>());
		UnigramTopicInferencer inf = umm.getInferencer();
		int numIterForInf = 30;
		int thinningForInf = 1;
		int burninForInf = 5;
		int[] queryTopics = new int[queryIList.size()];
		double[][] queryTopicScores = new double[queryIList.size()][umm.getNumTopics()];
		int[] paraTopics = new int[paraIList.size()];
		double[][] paraTopicScores = new double[paraIList.size()][umm.getNumTopics()];
		for(int i=0; i<queryIList.size(); i++){
			Instance queryIns = queryIList.get(i);
			/*
			fet = (FeatureSequence) queryIns.getData();
			if(fet.size()==0){
				this.countNan++;
				this.nanIns.add(queryIns.getName().toString());
			}
			*/
			//System.out.println("Page ID: "+page.getPageId());
			queryTopics[i] = inf.inferInstanceTopic(queryIns, numIterForInf, thinningForInf, burninForInf);
			queryTopicScores[i] = inf.inferInstanceTopicScores(queryIns, numIterForInf, thinningForInf, burninForInf);
		}
		if(paraIList.size()!=umm.getData().size())
			throw new Exception("paralist size and lda topic assignment size dont match!");
		for(int j=0; j<paraIList.size(); j++){
			// Following if block ensures that we are picking the correct topic dist for current para instance
			if(paraIList.get(j).getName()!=umm.getData().get(j).instance.getName())
				throw new Exception("paraIList indices are not following the same order as ummlda instances");
			paraTopics[j] = umm.getTopicOfInstance(j);
			paraTopicScores[j] = umm.getTopicScoresOfInstance(j);
		}
		matrixAssignment(assign, ranks, queryIList, paraIList, queryTopicScores, paraTopicScores);
		
		HashMap<String, ArrayList<String>> assignment = convInsAssignToIDAssign(assign, page);
		HashMap<ArrayList<String>, Double> idRanks = convRanksToIDRanks(ranks, page);
		ResultForPage r = new ResultForPage();
		r.setQueryParaAssignment(assignment);
		r.setQueryParaRank(idRanks);
		r.setParaClusters(getParasClustersFromParatopics(paraTopics, paraIList, umm.getNumTopics()));
		return r;
	}
	private ResultForPage assignUsingKMeans(InstanceList paraIList, InstanceList queryIList, 
			Clustering clusters, CustomKMeans kmeans, Data.Page page){
		boolean debug = false;
		ResultForPage r = new ResultForPage();
		HashMap<Instance, ArrayList<Instance>> insAssign = new HashMap<Instance, ArrayList<Instance>>();
		HashMap<ArrayList<Instance>, Double> ranks = new HashMap<ArrayList<Instance>, Double>();
		for(Instance qIns:queryIList)
			insAssign.put(qIns, new ArrayList<Instance>());
		double[][] queryClusterDistMat = transposeMatrix(getClusterQueryDistances(queryIList, clusters, kmeans, new NormalizedDotProductMetric()));
		double[][] paraClusterDistMat = transposeMatrix(getClusterQueryDistances(paraIList, clusters, kmeans, new NormalizedDotProductMetric()));
	
		//matrixAssignment(insAssign, ranks, queryIList, paraIList, queryClusterDistMat, paraClusterDistMat);
		matrixAssignmentKM(insAssign, ranks, queryIList, paraIList, queryClusterDistMat, paraClusterDistMat);
		
		HashMap<String, ArrayList<String>> assignment = convInsAssignToIDAssign(insAssign, page);
		HashMap<ArrayList<String>, Double> idRanks = convRanksToIDRanks(ranks, page);
		r.setQueryParaAssignment(assignment);
		r.setQueryParaRank(idRanks);
		r.setParaClusters(formatClusterData(clusters));
		return r;
	}
	public HashMap<Integer, SparseVector> getLabelledClusterMeans(Clustering clusters, CustomKMeans kmeans){
		  HashMap<Integer, SparseVector> labelledMeans = new HashMap<Integer, SparseVector>();
		  ArrayList<SparseVector> clusterMeans = kmeans.getClusterMeans();
		  Metric metric = kmeans.getMetric();
		  InstanceList currCluster;
		  SparseVector mean;
		  for(int c=0; c<clusters.getNumClusters(); c++){
			  // c is label here
			  currCluster = clusters.getCluster(c);
			  mean = VectorStats.mean(currCluster);
			  for(SparseVector cMean:clusterMeans){
				  if(metric.distance(mean, cMean)<0.0000001)
					  labelledMeans.put(c, cMean);
			  }
		  }
		  return labelledMeans;
	  }
	private double[][] getClusterQueryDistances(InstanceList inslist, Clustering clusters, CustomKMeans kmeans, Metric metric){
		SparseVector vec;
		Instance i;
		for(int ins=0; ins<inslist.size(); ins++){
			i = inslist.get(ins);
			vec = (SparseVector) i.getData();
			if(vec.numLocations()==0){
				this.countNan++;
				this.nanIns.add(i.getName().toString());
			}
		}
		double[][] distanceMatrix = new double [clusters.getNumClusters()][inslist.size()];
		ArrayList<SparseVector> clusterMeans = kmeans.getClusterMeans();
		double dist;
		for(int c=0; c<distanceMatrix.length; c++){
			for(int ins=0; ins<distanceMatrix[0].length; ins++){
				dist = metric.distance(clusterMeans.get(c), (SparseVector) inslist.get(ins).getData());
				if(dist<0.0000001)
					System.out.println("ins no. "+ins+":"+
				inslist.get(ins).getName().toString()+" -- "+"cluster "+c+" dist is "+dist);
				distanceMatrix[c][ins] = dist;
			}
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
	private double getKSD(double[] p, double[] q){
		double result = 0, sump = 0, sumq = 0, sumpTemp = 0, sumqTemp = 0, diff;
		for(int i=0; i<p.length; i++){
			sump+=p[i];
			sumq+=q[i];
		}
		result = 0;
		for(int i=0; i<p.length; i++){
			sumpTemp+=p[i];
			sumqTemp+=q[i];
			diff = Math.abs(sumpTemp/sump - sumqTemp/sumq);
			if(diff>result)
				result = diff;
		}
		return result;
	}
	private double getBD(double[] p, double[] q){
		double result = 0;
		for(int i=0; i<p.length; i++){
			if(p[i]<0 || q[i]<0)
				continue;
			result+=Math.sqrt(p[i]*q[i]);
		}
		result = - Math.log(result);
		return result;
	}
	private double getChiD(double[] p, double[] q){
		double result = 0;
		for(int i=0; i<p.length; i++)
			result+= ((p[i]-q[i])*(p[i]-q[i]))/(p[i]+q[i]);
		return result;
	}
	private InstanceList convertQueriesToIList(ArrayList<String> queries){
		InstanceList qIList = new InstanceList(SingleRun.buildPipe(this.model));
		for(String q:queries){
			Instance secIns = new Instance(q.replaceAll("%20|/", " "), null, q, q);
			qIList.addThruPipe(secIns);
		}
		if(this.tw == 1){ //TD
			if(qIList.getDataClass()==FeatureVector.class){
				SparseVector vec;
				for(Instance p:qIList){
					double d = 0;
					vec = (SparseVector) p.getData();
					for(double freq:vec.getValues())
						d+=freq;
					for(int index:vec.getIndices())
						vec.setValue(index, vec.value(index)/d);
				}
			}
			else{
				System.out.println("TF variants are only implemented for SparseVectors. Using default tw value 0");
			}
		} else if(this.tw == 2){ //TDS
			final double beta = 10;
			double v = 0;
			double[] d = new double[qIList.size()];
			if(qIList.getDataClass()==FeatureVector.class){
				SparseVector vec;
				for(int i=0; i<qIList.size(); i++){
					vec = (SparseVector) qIList.get(i).getData();
					for(double freq:vec.getValues())
						d[i]+=freq;
					v+=d[i];
				}
				for(int i=0; i<qIList.size(); i++){
					vec = (SparseVector) qIList.get(i).getData();
					for(int index:vec.getIndices())
						vec.setValue(index, (vec.value(index)+beta)/(v*d[i]*beta));
				}
			}
			else{
				System.out.println("TF variants are only implemented for SparseVectors. Using default tw value 0");
			}
		}
		return qIList;
	}
	private InstanceList convertParasToIList(ArrayList<Data.Paragraph> paraObjs){
		InstanceList iListPara = new InstanceList(SingleRun.buildPipe(this.model));
		for(Data.Paragraph paraObj:paraObjs){
			Instance paraIns = new Instance(paraObj.getTextOnly(), null, paraObj.getParaId(), paraObj.getTextOnly());
			iListPara.addThruPipe(paraIns);
		}
		if(this.tw == 1){ //TD
			if(iListPara.getDataClass()==FeatureVector.class){
				SparseVector vec;
				for(Instance p:iListPara){
					double d = 0;
					vec = (SparseVector) p.getData();
					for(double freq:vec.getValues())
						d+=freq;
					for(int index:vec.getIndices())
						vec.setValue(index, vec.value(index)/d);
				}
			}
			else{
				System.out.println("TF variants are only implemented for SparseVectors. Using default tw value 0");
			}
		} else if(this.tw == 2){ //TDS
			final double beta = 10;
			double v = 0;
			double[] d = new double[iListPara.size()];
			if(iListPara.getDataClass()==FeatureVector.class){
				SparseVector vec;
				for(int i=0; i<iListPara.size(); i++){
					vec = (SparseVector) iListPara.get(i).getData();
					for(double freq:vec.getValues())
						d[i]+=freq;
					v+=d[i];
				}
				for(int i=0; i<iListPara.size(); i++){
					vec = (SparseVector) iListPara.get(i).getData();
					for(int index:vec.getIndices())
						vec.setValue(index, (vec.value(index)+beta)/(v*d[i]*beta));
				}
			}
			else{
				System.out.println("TF variants are only implemented for SparseVectors. Using default tw value 0");
			}
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
	public int[] getMinValIndex(double[][] kldivMat, double maxVal){
		int[] minindex = new int[2];
		double currVal, minVal=maxVal+2;
		for(int i=0; i<kldivMat.length; i++){
			for(int j=0; j<kldivMat[0].length; j++){
				currVal = kldivMat[i][j];
				if(currVal<minVal){
					minVal = currVal;
					minindex[0] = i;
					minindex[1] = j;
				}
			}
		}
		return minindex;
	}
	public double getMaxVal(double[][] kldivMat){
		double currVal, maxVal=0;
		for(int i=0; i<kldivMat.length; i++){
			for(int j=0; j<kldivMat[0].length; j++){
				currVal = kldivMat[i][j];
				if(currVal>maxVal)
					maxVal = currVal;
			}
		}
		return maxVal;
	}
	public boolean checkIfDone(double[][] kldivMat, double maxVal){
		for(int i=0; i<kldivMat.length; i++){
			for(int j=0; j<kldivMat[0].length; j++){
				if(kldivMat[i][j]<=maxVal)
					return false;
			}
		}
		return true;
	}
	public boolean checkIfDone(double[][] distMat){
		for(int i=0; i<distMat.length; i++){
			for(int j=0; j<distMat[0].length; j++){
				if(!Double.isNaN(distMat[i][j]) && distMat[i][j]>0)
					return false;
			}
		}
		return true;
	}
	private double[][] transposeMatrix(double[][] matrix){
		double[][] transposed = new double[matrix[0].length][matrix.length];
		for(int r=0; r<transposed.length; r++){
			for(int c=0; c<transposed[0].length; c++)
				transposed[r][c] = matrix[c][r];
		}
		return transposed;
	}
	private double[][] normalizeMatrix(double[][] mat){ // normalize each row of mat to 1
		int rowCnt = mat.length;
		int colCnt = mat[0].length;
		double rowSum;
		double[][] norm = new double[rowCnt][colCnt];
		for(int i=0; i<rowCnt; i++){
			rowSum = 0;
			for(int j=0; j<colCnt; j++)
				rowSum+=mat[i][j];
			for(int j=0; j<colCnt; j++)
				norm[i][j] = mat[i][j]/rowSum;
		}
		return norm;
	}
	private double[][] reciprocMatrix(double[][] mat){ // recip[x][y] = (sum of mat[x])/mat[x][y]
		int rowCnt = mat.length;
		int colCnt = mat[0].length;
		double rowSum;
		boolean gotNaN;
		double[][] recip = new double[rowCnt][colCnt];
		for(int i=0; i<rowCnt; i++){
			rowSum = 0;
			gotNaN = false;
			for(int j=0; j<colCnt; j++)
				rowSum+=mat[i][j];
			for(int j=0; j<colCnt; j++){
				if(mat[i][j]<0.00001){
					gotNaN = true;
					break;
				}
				recip[i][j] = rowSum/mat[i][j];
			}
			if(gotNaN){
				for(int j=0; j<colCnt; j++){
					if(mat[i][j]<0.00001)
						recip[i][j] = 0.99;
					else
						recip[i][j] = 0.01;
				}
			}
		}
		return recip;
	}
	public void matrixAssignmentKM(HashMap<Instance, ArrayList<Instance>> assign, HashMap<ArrayList<Instance>, Double> ranks,
			InstanceList queryIList, InstanceList paraIList, double[][] queryDistMat, double[][] paraDistMat){
		double[][] queryMat = this.normalizeMatrix(queryDistMat);
		double[][] paraMat = this.normalizeMatrix(paraDistMat);
		
		// We should be able to treat normalized distance dist. as prob. dist.
		// because KLDiv is suppose to find similarity (information lost) between
		// two dist. Right??
		
		double[][] paraQueryMat = new double[paraIList.size()][queryIList.size()];
		double matVal;
		for(int n=0; n<paraIList.size(); n++){
			for(int m=0; m<queryIList.size(); m++){
				switch(RunExperiment.ASSIGN_METHOD){
				case 1:
					matVal = getKLdiv(paraMat[n], queryMat[m]); // This will give perfect KLDiv results
					//klVal = getKLdiv(queryMat[m], paraMat[n]); // This will give some -ve KLDive results
					paraQueryMat[n][m] = matVal;
					break;
				case 2:
					matVal = getKSD(paraMat[n], queryMat[m]);
					paraQueryMat[n][m] = matVal;
					break;
				case 3:
					matVal = getBD(paraMat[n], queryMat[m]);
					paraQueryMat[n][m] = matVal;
					break;
				case 4:
					matVal = getChiD(paraMat[n], queryMat[m]);
					paraQueryMat[n][m] = matVal;
					break;
				default:
					System.out.println("Wrong assign method");
				}
			}
		}
		System.out.println("\nparaQuery matrix\n----------------\n");
		for(int n=0; n<paraIList.size(); n++){
			for(int m=0; m<queryIList.size(); m++){
				System.out.printf("%16f ",paraQueryMat[n][m]);
			}
			System.out.println();
		}
		matrixAssignment(assign, ranks, queryIList, paraIList, paraQueryMat);
	}
	public void matrixAssignment(HashMap<Instance, ArrayList<Instance>> assign, HashMap<ArrayList<Instance>, Double> ranks,
			InstanceList queryIList, InstanceList paraIList, double[][] queryMat, double[][] paraMat){
		double[][] paraQueryMat = new double[paraIList.size()][queryIList.size()];
		double matVal;
		for(int n=0; n<paraIList.size(); n++){
			for(int m=0; m<queryIList.size(); m++){
				switch(RunExperiment.ASSIGN_METHOD){
				case 1:
					matVal = getKLdiv(paraMat[n], queryMat[m]); // This will give perfect KLDiv results
					//klVal = getKLdiv(queryMat[m], paraMat[n]); // This will give some -ve KLDive results
					paraQueryMat[n][m] = matVal;
					break;
				case 2:
					matVal = getKSD(paraMat[n], queryMat[m]);
					paraQueryMat[n][m] = matVal;
					break;
				case 3:
					matVal = getBD(paraMat[n], queryMat[m]);
					paraQueryMat[n][m] = matVal;
					break;
				case 4:
					matVal = getChiD(paraMat[n], queryMat[m]);
					paraQueryMat[n][m] = matVal;
					break;
				default:
					System.out.println("Wrong assign method");
				}
			}
		}
		System.out.println("\nparaQuery matrix\n----------------\n");
		for(int n=0; n<paraIList.size(); n++){
			
			for(int m=0; m<queryIList.size(); m++){
				System.out.printf("%16f ",paraQueryMat[n][m]);
			}
			System.out.println();
		}
		matrixAssignment(assign, ranks, queryIList, paraIList, paraQueryMat);
	}
	public void matrixAssignment(HashMap<Instance, ArrayList<Instance>> assign, HashMap<ArrayList<Instance>, Double> ranks,
			InstanceList queryIList, InstanceList paraIList, double[][] paraQueryMatrix){
		ArrayList<Instance> ranksKey;
		Instance secIns, pIns;
		double matrixVal;
		for(int n=0; n<paraIList.size(); n++){
			for(int m=0; m<queryIList.size(); m++){
				secIns = queryIList.get(m);
				pIns = paraIList.get(n);
				ranksKey = new ArrayList<Instance>();
				ranksKey.add(0, secIns);
				ranksKey.add(1, pIns);
				matrixVal = paraQueryMatrix[n][m];
				ranks.put(ranksKey, 1/matrixVal);
				/*
				if(matrixVal>0.0001){
					ranks.put(ranksKey, 1/matrixVal);
				} else{
					System.out.println("Debug here");
				}
				*/
			}
		}
		boolean[] isParaAssigned = new boolean[paraIList.size()];
		for(int pos=0; pos<isParaAssigned.length; pos++)
			isParaAssigned[pos] = false;
		int[] bestQueryForPara = new int[paraIList.size()];
		double[] queryVals = new double[queryIList.size()];
		for(int n=0; n<paraIList.size(); n++){
			int bestQueryInd=0;
			double minVal = 99999.0;
			queryVals = paraQueryMatrix[n];
			for(int m=0; m<queryIList.size(); m++){
				if(queryVals[m]<minVal){
					bestQueryInd = m;
					minVal = queryVals[m];
				}
			}
			bestQueryForPara[n] = bestQueryInd;
		}
		boolean isDone = false;
		double maxKLVal = getMaxVal(paraQueryMatrix);
		int[] minInd = new int[2];
		while(!isDone){
			minInd = getMinValIndex(paraQueryMatrix, maxKLVal);
			assign.get(queryIList.get(minInd[1])).add(paraIList.get(minInd[0]));
			isParaAssigned[minInd[0]] = true;
			for(int r=0; r<paraQueryMatrix.length; r++){
				for(int c=0; c<paraQueryMatrix[0].length; c++){
					if(r==minInd[0] || c==minInd[1])
						paraQueryMatrix[r][c] = maxKLVal+1;
				}
			}	
			isDone = checkIfDone(paraQueryMatrix, maxKLVal);
		}
		for(int n=0; n<isParaAssigned.length; n++){
			if(!isParaAssigned[n])
				assign.get(queryIList.get(bestQueryForPara[n])).add(paraIList.get(n));
		}
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
