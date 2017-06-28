package edu.unh.cs.treccar.playground;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Formatter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Locale;
import java.util.TreeSet;
import java.util.regex.Pattern;

import cc.mallet.pipe.CharSequence2TokenSequence;
import cc.mallet.pipe.FeatureSequence2FeatureVector;
import cc.mallet.pipe.Input2CharSequence;
import cc.mallet.pipe.Pipe;
import cc.mallet.pipe.PrintInputAndTarget;
import cc.mallet.pipe.SerialPipes;
import cc.mallet.pipe.TokenSequence2FeatureSequence;
import cc.mallet.pipe.TokenSequenceLowercase;
import cc.mallet.pipe.TokenSequenceRemoveStopwords;
import cc.mallet.topics.ParallelTopicModel;
import cc.mallet.topics.TopicInferencer;
import cc.mallet.types.Alphabet;
import cc.mallet.types.FeatureSequence;
import cc.mallet.types.IDSorter;
import cc.mallet.types.Instance;
import cc.mallet.types.InstanceList;
import cc.mallet.types.LabelSequence;
import edu.unh.cs.treccar.Data;
import edu.unh.cs.treccar.playground.AssignParagraphs.SectionPathID;

public class TopicModel {
	public int topicNum = MeasurePerformanceMain.TOPICNUMBER;
	private ArrayList<AssignParagraphs.SectionPathID> sectionList;
	private Iterable<Data.Paragraph> paraList;
	private HashMap<AssignParagraphs.SectionPathID, ArrayList<Data.Paragraph>> candidate = new HashMap<AssignParagraphs.SectionPathID, ArrayList<Data.Paragraph>>();
	
	public TopicModel(ArrayList<AssignParagraphs.SectionPathID> sectionList, Iterable<Data.Paragraph> paraList){
		this.sectionList = sectionList;
		this.paraList = paraList;
	}
	public HashMap<AssignParagraphs.SectionPathID, ArrayList<Data.Paragraph>> modelTopics(){
		HashMap<AssignParagraphs.SectionPathID, double[]> secToTopicMap = new HashMap<AssignParagraphs.SectionPathID, double[]>();
		HashMap<Data.Paragraph, double[]> paraToTopicMap = new HashMap<Data.Paragraph, double[]>();
		InstanceList iList = new InstanceList(TopicModel.buildPipe());
		ArrayList<Instance> rawInstanceList = new ArrayList<Instance>();
		ArrayList<Data.Paragraph> listOfParagraphs = new ArrayList<Data.Paragraph>();
		for(Data.Paragraph para : this.paraList){
			Instance instance = new Instance(para.getTextOnly(), "", para.getParaId(), para.getTextOnly());
			rawInstanceList.add(instance);
			listOfParagraphs.add(para);
		}
		iList.addThruPipe(rawInstanceList.iterator());
		int numTopics = this.topicNum;
		ParallelTopicModel model = new ParallelTopicModel(numTopics, 1.0, 0.01);
		model.addInstances(iList);
		model.setNumThreads(2);
		try {
			model.estimate();
			
	        //Infer topics of section headings
	        InstanceList sectionInstances = new InstanceList(iList.getPipe());
	        TopicInferencer inferencer = model.getInferencer();
	        for(AssignParagraphs.SectionPathID sec : this.sectionList){
	        	sectionInstances.addThruPipe(new Instance(sec.getSectionPathID(), "", sec.getSectionPathID(), ""));
	        }
	        for(int i=0; i<sectionInstances.size(); i++){
	        	//get topic dist of current instance
	        	double[] topicDist = inferencer.getSampledDistribution(sectionInstances.get(i), 10, 1, 5);
	        	//get sectionpathid object from sectionList using instance name
	        	//map it with the highest topic index
	        	for(AssignParagraphs.SectionPathID secID : this.sectionList){
	        		String candSecId = secID.getSectionPathID();
	        		String instanceSecId = sectionInstances.get(i).getName().toString();
	        		if(candSecId.equals(instanceSecId)){
	        			secToTopicMap.put(secID, topicDist);
	        			break;
	        		}
	        	}
	        }
	        for(int i=0; i<iList.size(); i++){
	        	double[] topicDistPara = model.getTopicProbabilities(i);
	        	/*
	        	int indexHighestPara = 0;
	        	double maxProbPara = 0.0;
	        	for(int j=1; j<topicDistPara.length; j++){
	        		if(topicDistPara[j]>maxProbPara){
	        			indexHighestPara = j;
	        			maxProbPara = topicDistPara[j];
	        		}
	        	}
	        	*/
	        	for(Data.Paragraph p : listOfParagraphs){
	        		String candPId = p.getParaId();
	        		String instancePId = iList.get(i).getName().toString();
	        		if(candPId.equals(instancePId)){
	        			paraToTopicMap.put(p, topicDistPara);
	        			break;
	        		}
	        	}
	        }
	        prepareCandidate(secToTopicMap, paraToTopicMap);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return this.candidate;
	}
	private void prepareCandidate(HashMap<AssignParagraphs.SectionPathID, double[]> secTopic, HashMap<Data.Paragraph, double[]> paraTopic){
		HashMap<AssignParagraphs.SectionPathID, ArrayList<Data.Paragraph>> cand = new HashMap<AssignParagraphs.SectionPathID, ArrayList<Data.Paragraph>>();
		for(Data.Paragraph para : paraTopic.keySet()){
			double[] paraTopicDist = paraTopic.get(para);
			AssignParagraphs.SectionPathID bestSecID = (AssignParagraphs.SectionPathID)secTopic.keySet().toArray()[0];
			double minKLdiv = 99999.0;
			for(AssignParagraphs.SectionPathID sec : secTopic.keySet()){
				double currentKLdiv = getKLdiv(paraTopic.get(para), secTopic.get(sec));
				if(minKLdiv>currentKLdiv){
					minKLdiv=currentKLdiv;
					bestSecID = sec;
				}
			}
			//map bestSecID with para
			if(cand.containsKey(bestSecID)){
				cand.get(bestSecID).add(para);
			} else{
				ArrayList<Data.Paragraph> currentParaList = new ArrayList<Data.Paragraph>();
				currentParaList.add(para);
				cand.put(bestSecID, currentParaList);
			}
		}
		this.candidate = cand;
	}
	public HashMap<AssignParagraphs.SectionPathID, ArrayList<Data.Paragraph>> getCandidate() {
		return this.candidate;
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
	public static Pipe buildPipe(){
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
}
