package edu.unh.cs.treccar.playground.topics;

import cc.mallet.types.Alphabet;
import cc.mallet.types.FeatureSequence;
import cc.mallet.types.Instance;
import cc.mallet.util.Randoms;

public class UnigramTopicInferencer {
	protected int numTopics; 

	// These values are used to encode type/topic counts as
	//  count/topic pairs in a single int.
	protected int topicMask;
	protected int topicBits;
	
	protected int numTypes;
	
	//protected double[] alpha;
	protected double beta;
	protected double betaSum;
	
	protected int[][] typeTopicCounts;
	protected int[] tokensPerTopic;
	
	Alphabet alphabet;
	
	protected Randoms random = null;
	
	//double smoothingOnlyMass = 0.0;
	//double[] cachedCoefficients;
	
	public UnigramTopicInferencer (int[][] typeTopicCounts, int[] tokensPerTopic, Alphabet alphabet, double betaSum) {

		this.tokensPerTopic = tokensPerTopic;
		this.typeTopicCounts = typeTopicCounts;

		this.alphabet = alphabet;

		numTopics = tokensPerTopic.length;
		numTypes = typeTopicCounts.length;

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
		this.betaSum = betaSum;
		this.beta = this.betaSum/this.alphabet.size();

		/*
		cachedCoefficients = new double[numTopics];

		for (int topic=0; topic < numTopics; topic++) {
			smoothingOnlyMass += alpha[topic] * beta / (tokensPerTopic[topic] + betaSum);
			cachedCoefficients[topic] =  alpha[topic] / (tokensPerTopic[topic] + betaSum);
		}
		*/

		random = new Randoms();
}
	
	// equivalent to getSampledDistribution() of TopicInferencer
	public int inferInstanceTopic(Instance instance, int numIterations, int thinning, int burnIn) {
		double[] topicScores = new double[numTopics];
		double[] smoothedScores = new double[numTopics];
		double[] logpkList = new double[numTopics];
		int inferredTopic = -1;
		boolean printMessages = true;
		
		FeatureSequence tokens = (FeatureSequence) instance.getData();
		int docLength = tokens.size();
		
		double topicScoreSum = 0, normalizingC, currentMin, currentMax, minLog, maxLog;
		for(int k=0; k<numTopics; k++){
			double sumOfLog = 0;
			for(int t=0; t<docLength; t++)
				sumOfLog+=Math.log((typeTopicCounts[tokens.getIndexAtPosition(t)][k]>>topicBits)+beta);
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
		/*
		double sumSmooth = numTopics-(numTopics-1)*UnigramTopicModel.LAMBDA;
		for(int k=0; k<numTopics; k++){
			smoothedScores[k] = (UnigramTopicModel.LAMBDA*topicScores[k]+(1-UnigramTopicModel.LAMBDA))/sumSmooth;
			if(printMessages)
				System.out.println("smoothedScore("+k+")="+smoothedScores[k]);
		}
		*/
		double sample = random.nextUniform();
		if(printMessages)
			System.out.println("Sample="+sample);
		while (sample > 0.0) {
			inferredTopic++;
			//sample -= smoothedScores[inferredTopic];
			sample -= topicScores[inferredTopic];
		}
		if(inferredTopic==-1)
			throw new IllegalStateException("New topic not sampled");
		if(printMessages)
			System.out.println("New topic="+inferredTopic);
		return inferredTopic;
	}

}
