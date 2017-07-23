package edu.unh.cs.treccar.playground.topics;

import java.util.*;
import java.util.logging.*;
import java.util.zip.*;
import java.io.*;
import java.text.NumberFormat;

import cc.mallet.topics.*;
import cc.mallet.types.*;
import cc.mallet.util.*;
import edu.unh.cs.treccar.playground.RunExperiment;

public class UnigramTopicModel implements Serializable {

	private static Logger logger = MalletLogger.getLogger(UnigramTopicModel.class.getName());
	
	// the training instances and their topic assignments
	protected ArrayList<TopicAssignment> data;  

	// the alphabet for the input data
	protected Alphabet alphabet; 

	// the alphabet for the topics
	protected LabelAlphabet topicAlphabet; 
	
	// The number of topics requested
	protected int numTopics;
	
	// These values are used to encode type/topic counts as
	//  count/topic pairs in a single int.
	public int topicMask;
	public int topicBits;
	
	// The size of the vocabulary
	protected int numTypes;

	// Prior parameters
	//protected double alpha;	 // Dirichlet(alpha,alpha,...) is the distribution over topics
	//protected double alphaSum;
	protected double beta;   // Prior on per-topic multinomial distribution over words
	protected double betaSum;
	public static final double DEFAULT_BETASUM = 260; // if beta<0 then beta = BETA_NOM/V
	protected double lambda = 0.0;
	
	// An array to put the topic counts for the current document. 
	// Initialized locally below.  Defined here to avoid
	// garbage collection overhead.
	protected int[] oneDocTopicCounts; // indexed by <document index, topic index>

	// Statistics needed for sampling.
	protected int[][] typeTopicCounts; // indexed by <feature index, topic index>
	protected int[] tokensPerTopic; // indexed by <topic index>
	protected int[] docTopicSequence; // indexed by <doc index> for unigram
	protected double[][] docTopicScores; // indexed by <doc index, topic index>

	public int showTopicsInterval = 50;
	public int wordsPerTopic = 10;
	
	protected Randoms random;
	protected NumberFormat formatter;
	protected boolean printLogLikelihood = false;
	protected boolean printMessages = false;
	
	public UnigramTopicModel (int numberOfTopics) {
		this (numberOfTopics, numberOfTopics, DEFAULT_BETASUM);
	}
	
	public UnigramTopicModel (int numberOfTopics, double alphaSum, double beta) {
		this (numberOfTopics, alphaSum, beta, new Randoms());
	}
	
	private static LabelAlphabet newLabelAlphabet (int numTopics) {
		LabelAlphabet ret = new LabelAlphabet();
		for (int i = 0; i < numTopics; i++)
			ret.lookupIndex("topic"+i);
		return ret;
	}
	
	public UnigramTopicModel (int numberOfTopics, double alphaSum, double beta, Randoms random) {
		this (newLabelAlphabet (numberOfTopics), alphaSum, beta, random);
	}
	
	public UnigramTopicModel (LabelAlphabet topicAlphabet, double alphaSum, double betaSum, Randoms random)
	{
		this.data = new ArrayList<TopicAssignment>();
		this.topicAlphabet = topicAlphabet;
		this.numTopics = topicAlphabet.size();

		//this.alphaSum = alphaSum;
		//this.alpha = alphaSum / numTopics;
		this.betaSum = betaSum;
		this.random = random;
		
		if (Integer.bitCount(numTopics) == 1) {
			// exact power of 2
			topicMask = numTopics - 1;
			topicBits = Integer.bitCount(topicMask);
		}
		else {
			// otherwise add an extra bit
			topicMask = Integer.highestOneBit(numTopics) * 2 - 1;
			topicBits = Integer.bitCount(topicMask);
		}
		
		oneDocTopicCounts = new int[numTopics];
		tokensPerTopic = new int[numTopics];
		
		formatter = NumberFormat.getInstance();
		formatter.setMaximumFractionDigits(5);

		logger.info("Unigram Topic Model: " + numTopics + " topics");
	}
	
	public Alphabet getAlphabet() { return alphabet; }
	public LabelAlphabet getTopicAlphabet() { return topicAlphabet; }
	public int getNumTopics() { return numTopics; }
	public ArrayList<TopicAssignment> getData() { return data; }
	
	public void setTopicDisplay(int interval, int n) {
		this.showTopicsInterval = interval;
		this.wordsPerTopic = n;
	}

	public void setRandomSeed(int seed) {
		random = new Randoms(seed);
	}
	
	public int[][] getTypeTopicCounts() { return typeTopicCounts; }
	public int[] getTopicTotals() { return tokensPerTopic; }

	public UnigramTopicInferencer getInferencer(){
		UnigramTopicInferencer inf = new UnigramTopicInferencer(this.typeTopicCounts, this.tokensPerTopic, 
				this.alphabet, this.betaSum);
		return inf;
	}
	
	public int getTopicOfInstance(int instanceIndex){
		// equivalent to getTopicProbabilities() of CustomLDA
		return this.docTopicSequence[instanceIndex];
	}
	public double[] getTopicScoresOfInstance(int instanceIndex){
		return this.docTopicScores[instanceIndex];
	}
	
	// addInstances method for unigram style //
	public void addInstances (InstanceList training) {

		alphabet = training.getDataAlphabet();
		numTypes = alphabet.size();
		
		beta = betaSum / numTypes;

		typeTopicCounts = new int[numTypes][numTopics];
		docTopicSequence = new int[training.size()];
		docTopicScores = new double[training.size()][numTopics];
		initAll();

		System.out.println("UMM betaSum="+this.betaSum+", beta="+this.beta);
		int doc = 0;

		for (Instance instance : training) {
			FeatureSequence tokens = (FeatureSequence) instance.getData();
			LabelSequence topicSequence =
					new LabelSequence(topicAlphabet, new int[ tokens.size() ]);

			int[] topics = topicSequence.getFeatures();
			int topic = random.nextInt(numTopics);
			docTopicSequence[doc] = topic;
			docTopicScores[doc][topic] = 1.0;
			for (int position = 0; position < tokens.size(); position++) {
				topics[position] = topic;
				tokensPerTopic[topic]++;

				int type = tokens.getIndexAtPosition(position);
				typeTopicCounts[type][topic]+=topicMask+1;
			}

			TopicAssignment t = new TopicAssignment (instance, topicSequence);
			data.add (t);
			doc++;
		}
	}
	
	private void initAll(){
		for(int topic=0; topic<numTopics; topic++){
			for(int doc=0; doc<docTopicScores.length; doc++)
				docTopicScores[doc][topic]=0.0;
			for(int type=0; type<numTypes; type++)
				typeTopicCounts[type][topic] = topic;
		}
	}
	
	public void sample (int iterations) throws IOException {
		int lambdaStep = iterations/11;
		int lambdaCount = 0;
		for (int iteration = 1; iteration <= iterations; iteration++) {
			long iterationStart = System.currentTimeMillis();
			lambdaCount++;
			if(lambdaCount==lambdaStep){
				lambdaCount = 0;
				if(this.lambda<0.99)
					this.lambda+=0.1;
				else if(this.lambda>0.99 && this.lambda<1.0)
					this.lambda=1.0;
				//System.out.println("iter="+iteration+",lambda="+this.lambda);
			}
			if(printMessages){
				System.out.print("Iteration "+iteration+", topicDocSeq array=[");
				for(int doci=0; doci<docTopicSequence.length; doci++)
					System.out.print(docTopicSequence[doci]+" ");
				System.out.print("]\n");
			}
			// Loop over every document in the corpus
			for (int doc = 0; doc < data.size(); doc++) {
				FeatureSequence tokenSequence =
					(FeatureSequence) data.get(doc).instance.getData();
				LabelSequence topicSequence =
					(LabelSequence) data.get(doc).topicSequence;

				sampleTopicsForOneDoc (tokenSequence, topicSequence, doc, this.lambda);
				if(printMessages){
					System.out.print("Para "+doc+" topic seq: [");
					for(int i=0; i<data.get(doc).topicSequence.getLength(); i++){
						System.out.print(data.get(doc).topicSequence.getIndexAtPosition(i)+" ");
					}
					System.out.print("]\n");
				}
				
			}
		
            long elapsedMillis = System.currentTimeMillis() - iterationStart;
			logger.fine(iteration + "\t" + elapsedMillis + "ms\t");

			/*
			// Occasionally print more information
			if (showTopicsInterval != 0 && iteration % showTopicsInterval == 0) {
				logger.info("<" + iteration + "> Log Likelihood: " + modelLogLikelihood() + "\n" +
							topWords (wordsPerTopic));
			}
			*/

		}
	}
	
	protected void sampleTopicsForOneDoc (FeatureSequence tokenSequence, LabelSequence topicSequence,
			  int docIndex, double lambda) {
		int oldDocTopic = docTopicSequence[docIndex];
		
		int newTopic;
		int docLength = tokenSequence.getLength(); // T
		
		for(int position=0; position<docLength; position++){
			int currType = tokenSequence.getIndexAtPosition(position);
			typeTopicCounts[currType][oldDocTopic]-=topicMask+1;
		}
		// now typeTopicCounts is n_{k,-i}^(t)
		tokensPerTopic[oldDocTopic]-=docLength;
		// now tokensPerTopic is tc_{k,-i}
		docTopicSequence[docIndex] = -1;
		
		double[] topicScores = new double[numTopics];
		double[] smoothedScores = new double[numTopics];
		double[] logpkList = new double[numTopics];
		double topicScoreSum = 0, normalizingC, currentMin, currentMax, minLog, maxLog;
		
		//-- p'(k) distribution is calculated in log space with normalization --//
		
		if(printMessages)
			System.out.println("beta="+beta+", betaSum="+betaSum+", V="+numTypes);
		for(int k=0; k<numTopics; k++){
			double sumOfLog = 0, logOfpk;
			for(int pos=0; pos<docLength; pos++)
				sumOfLog+=Math.log((typeTopicCounts[tokenSequence.getIndexAtPosition(pos)][k]>>topicBits)+beta);
			logpkList[k] = sumOfLog-docLength*Math.log(tokensPerTopic[k]+betaSum);
			
		}
		minLog = logpkList[0];
		maxLog = logpkList[0];
		for(int k=0; k<numTopics; k++){
			currentMin = logpkList[k];
			currentMax = logpkList[k];
			if(currentMin<minLog)
				minLog = currentMin;
			if(currentMax>maxLog)
				maxLog = currentMax;
		}
		normalizingC = (minLog+maxLog)/2;
		if(printMessages){
			System.out.println("minLog="+minLog+", maxLog="+maxLog);
			System.out.println("Normalizing C="+normalizingC);
		}
		for(int k=0; k<numTopics; k++){	
			topicScores[k] = Math.exp(logpkList[k]-normalizingC);
			topicScoreSum+=topicScores[k];
		}
		for(int k=0; k<numTopics; k++){
			topicScores[k] = topicScores[k]/topicScoreSum;
			if(printMessages)
				System.out.println("log p'("+k+")="+logpkList[k]+", topicScore("+k+")="+topicScores[k]);
		}
		
		
		//------------------------------------------------------------------//
		
		//-- p'(k) distribution is calculated directly without any normalization --//
		
		/*
		for(int k=0; k<numTopics; k++){
			double prod = 1.0;
			for(int t=0; t<docLength; t++)
				prod*=((typeTopicCounts[t][k]>>topicBits)+beta)/(tokensPerTopic[k]+betaSum);
			topicScores[k] = prod;
			topicScoreSum+=topicScores[k];
		}
		*/
		
		//-------------------------------------------------------------------------//
		if(RunExperiment.SMOOTHED_UMM){
			double sumSmooth = numTopics-(numTopics-1)*lambda;
			for(int k=0; k<numTopics; k++){
				smoothedScores[k] = (lambda*topicScores[k]+(1-lambda))/sumSmooth;
				if(printMessages)
					System.out.println("smoothedScore("+k+")="+smoothedScores[k]);
			}
			docTopicScores[docIndex] = smoothedScores;
		} else{
			docTopicScores[docIndex] = topicScores;
		}
		//double sample = random.nextUniform()*topicScoreSum;
		double sample = random.nextUniform();
		if(printMessages)
			System.out.println("Sample="+sample);
		newTopic = -1;
		while (sample > 0.0) {
			newTopic++;
			if(RunExperiment.SMOOTHED_UMM)
				sample -= smoothedScores[newTopic];
			else
				sample -= topicScores[newTopic];
		}
		if(newTopic==-1)
			throw new IllegalStateException("New topic not sampled");
		if(printMessages)
			System.out.println("New topic="+newTopic);
		int[] oneDocTopics = topicSequence.getFeatures();
		for(int position=0; position<docLength; position++){
			oneDocTopics[position] = newTopic;
			int currType = tokenSequence.getIndexAtPosition(position);
			typeTopicCounts[currType][newTopic]+=topicMask+1;
		}
		// now typeTopicCounts is n_{k}^(t)
		tokensPerTopic[newTopic]+=docLength;
		// now tokensPerTopic is tc_{k}
		docTopicSequence[docIndex] = newTopic;
	}
	private static final long serialVersionUID = 1;
}
