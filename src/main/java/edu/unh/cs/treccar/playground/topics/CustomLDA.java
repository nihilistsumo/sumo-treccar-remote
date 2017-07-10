/* Copyright (C) 2005 Univ. of Massachusetts Amherst, Computer Science Dept.
   This file is part of "MALLET" (MAchine Learning for LanguagE Toolkit).
   http://www.cs.umass.edu/~mccallum/mallet
   This software is provided under the terms of the Common Public License,
   version 1.0, as published by http://www.opensource.org.	For further
   information, see the file `LICENSE' included with this distribution. */

package edu.unh.cs.treccar.playground.topics;

import java.util.*;
import java.util.logging.*;
import java.util.zip.*;
import java.io.*;
import java.text.NumberFormat;

import cc.mallet.topics.*;
import cc.mallet.types.*;
import cc.mallet.util.*;


/**
 * A simple implementation of Latent Dirichlet Allocation using Gibbs sampling.
 * This code is slower than the regular Mallet LDA implementation, but provides a 
 *  better starting place for understanding how sampling works and for 
 *  building new topic models.
 * 
 * @author David Mimno, Andrew McCallum
 */

public class CustomLDA implements Serializable {

	private static Logger logger = MalletLogger.getLogger(CustomLDA.class.getName());
	
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
	protected double alpha;	 // Dirichlet(alpha,alpha,...) is the distribution over topics
	protected double alphaSum;
	protected double beta;   // Prior on per-topic multinomial distribution over words
	protected double betaSum;
	public static final double DEFAULT_BETA = 0.01;
	public static final double BETA_SUM = 260; // if beta<0 then beta = BETA_NOM/V
	public static final double LAMBDA = 0.9; // for smoothing in UMM
	
	// An array to put the topic counts for the current document. 
	// Initialized locally below.  Defined here to avoid
	// garbage collection overhead.
	protected int[] oneDocTopicCounts; // indexed by <document index, topic index>

	// Statistics needed for sampling.
	protected int[][] typeTopicCounts; // indexed by <feature index, topic index>
	protected int[] tokensPerTopic; // indexed by <topic index>
	protected int[] topicDocSequence; // indexed by <doc index> for unigram

	public int showTopicsInterval = 50;
	public int wordsPerTopic = 10;
	
	protected Randoms random;
	protected NumberFormat formatter;
	protected boolean printLogLikelihood = false;
	protected boolean printMessages = false;
	
	public CustomLDA (int numberOfTopics) {
		this (numberOfTopics, numberOfTopics, DEFAULT_BETA);
	}
	
	public CustomLDA (int numberOfTopics, double alphaSum, double beta) {
		this (numberOfTopics, alphaSum, beta, new Randoms());
	}
	
	private static LabelAlphabet newLabelAlphabet (int numTopics) {
		LabelAlphabet ret = new LabelAlphabet();
		for (int i = 0; i < numTopics; i++)
			ret.lookupIndex("topic"+i);
		return ret;
	}
	
	public CustomLDA (int numberOfTopics, double alphaSum, double beta, Randoms random) {
		this (newLabelAlphabet (numberOfTopics), alphaSum, beta, random);
	}
	
	public CustomLDA (LabelAlphabet topicAlphabet, double alphaSum, double beta, Randoms random)
	{
		this.data = new ArrayList<TopicAssignment>();
		this.topicAlphabet = topicAlphabet;
		this.numTopics = topicAlphabet.size();

		this.alphaSum = alphaSum;
		this.alpha = alphaSum / numTopics;
		this.beta = beta;
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

		logger.info("Simple LDA: " + numTopics + " topics");
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

	public void addInstances (InstanceList training) {

		alphabet = training.getDataAlphabet();
		numTypes = alphabet.size();
		
		if(beta<0)
			beta = CustomLDA.BETA_SUM/numTypes;
		betaSum = beta * numTypes;
		
		typeTopicCounts = new int[numTypes][numTopics];
		initTypeTopicCounts();

		int doc = 0;

		for (Instance instance : training) {
			doc++;

			FeatureSequence tokens = (FeatureSequence) instance.getData();
			LabelSequence topicSequence =
				new LabelSequence(topicAlphabet, new int[ tokens.size() ]);
			
			int[] topics = topicSequence.getFeatures();
			for (int position = 0; position < tokens.size(); position++) {

				int topic = random.nextInt(numTopics);
				topics[position] = topic;
				tokensPerTopic[topic]++;
				
				int type = tokens.getIndexAtPosition(position);
				//typeTopicCounts[type][topic]++;
				typeTopicCounts[type][topic] += topicMask+1;
			}

			TopicAssignment t = new TopicAssignment (instance, topicSequence);
			data.add (t);
		}

	}
	
	public void sample (int iterations) throws IOException {

		for (int iteration = 1; iteration <= iterations; iteration++) {

			long iterationStart = System.currentTimeMillis();

			// Loop over every document in the corpus
			for (int doc = 0; doc < data.size(); doc++) {
				FeatureSequence tokenSequence =
					(FeatureSequence) data.get(doc).instance.getData();
				LabelSequence topicSequence =
					(LabelSequence) data.get(doc).topicSequence;

				sampleTopicsForOneDoc (tokenSequence, topicSequence);
			}
		
            long elapsedMillis = System.currentTimeMillis() - iterationStart;
			logger.fine(iteration + "\t" + elapsedMillis + "ms\t");

			// Occasionally print more information
			if (showTopicsInterval != 0 && iteration % showTopicsInterval == 0) {
				logger.info("<" + iteration + "> Log Likelihood: " + modelLogLikelihood() + "\n" +
							topWords (wordsPerTopic));
			}

		}
	}
	
	protected void sampleTopicsForOneDoc (FeatureSequence tokenSequence,
										  FeatureSequence topicSequence) {

		int[] oneDocTopics = topicSequence.getFeatures();

		int[] currentTypeTopicCounts;
		int type, oldTopic, newTopic;
		double topicWeightsSum;
		int docLength = tokenSequence.getLength();

		int[] localTopicCounts = new int[numTopics];

		//		populate topic counts
		for (int position = 0; position < docLength; position++) {
			localTopicCounts[oneDocTopics[position]]++;
		}

		double score, sum;
		double[] topicTermScores = new double[numTopics];

		//	Iterate over the positions (words) in the document 
		for (int position = 0; position < docLength; position++) {
			type = tokenSequence.getIndexAtPosition(position);
			oldTopic = oneDocTopics[position];

			// Grab the relevant row from our two-dimensional array
			currentTypeTopicCounts = typeTopicCounts[type];

			//	Remove this token from all counts. 
			localTopicCounts[oldTopic]--;
			tokensPerTopic[oldTopic]--;
			assert(tokensPerTopic[oldTopic] >= 0) : "old Topic " + oldTopic + " below 0";
			currentTypeTopicCounts[oldTopic]-=topicMask+1;

			// Now calculate and add up the scores for each topic for this word
			sum = 0.0;
			
			// Here's where the math happens! Note that overall performance is 
			//  dominated by what you do in this loop.
			for (int topic = 0; topic < numTopics; topic++) {
				score =
					(alpha + localTopicCounts[topic]) *
					((beta + (currentTypeTopicCounts[topic] >> topicBits)) /
					 (betaSum + tokensPerTopic[topic]));
				sum += score;
				topicTermScores[topic] = score;
			}
			
			// Choose a random point between 0 and the sum of all topic scores
			double sample = random.nextUniform() * sum;

			// Figure out which topic contains that point
			newTopic = -1;
			while (sample > 0.0) {
				newTopic++;
				sample -= topicTermScores[newTopic];
			}

			// Make sure we actually sampled a topic
			if (newTopic == -1) {
				throw new IllegalStateException ("SimpleLDA: New topic not sampled.");
			}

			// Put that new topic into the counts
			oneDocTopics[position] = newTopic;
			localTopicCounts[newTopic]++;
			tokensPerTopic[newTopic]++;
			currentTypeTopicCounts[newTopic]+=topicMask+1;
		}
	}
	
	public double modelLogLikelihood() {
		double logLikelihood = 0.0;

		// The likelihood of the model is a combination of a 
		// Dirichlet-multinomial for the words in each topic
		// and a Dirichlet-multinomial for the topics in each
		// document.

		// The likelihood function of a dirichlet multinomial is
		//	 Gamma( sum_i alpha_i )	 prod_i Gamma( alpha_i + N_i )
		//	prod_i Gamma( alpha_i )	  Gamma( sum_i (alpha_i + N_i) )

		// So the log likelihood is 
		//	logGamma ( sum_i alpha_i ) - logGamma ( sum_i (alpha_i + N_i) ) + 
		//	 sum_i [ logGamma( alpha_i + N_i) - logGamma( alpha_i ) ]

		// Do the documents first

		int[] topicCounts = new int[numTopics];
		double[] topicLogGammas = new double[numTopics];
		int[] docTopics;

		for (int topic=0; topic < numTopics; topic++) {
			topicLogGammas[ topic ] = Dirichlet.logGamma( alpha );
		}
	
		for (int doc=0; doc < data.size(); doc++) {
			LabelSequence topicSequence = (LabelSequence) data.get(doc).topicSequence;

			docTopics = topicSequence.getFeatures();

			for (int token=0; token < docTopics.length; token++) {
				topicCounts[ docTopics[token] ]++;
			}

			for (int topic=0; topic < numTopics; topic++) {
				if (topicCounts[topic] > 0) {
					logLikelihood += (Dirichlet.logGamma(alpha + topicCounts[topic]) -
									  topicLogGammas[ topic ]);
				}
			}

			// subtract the (count + parameter) sum term
			logLikelihood -= Dirichlet.logGamma(alphaSum + docTopics.length);

			Arrays.fill(topicCounts, 0);
		}
	
		// add the parameter sum term
		logLikelihood += data.size() * Dirichlet.logGamma(alphaSum);

		// And the topics

		double logGammaBeta = Dirichlet.logGamma(beta);

		for (int type=0; type < numTypes; type++) {
			// reuse this array as a pointer

			topicCounts = typeTopicCounts[type];

			for (int topic = 0; topic < numTopics; topic++) {
				if (topicCounts[topic] == 0) { continue; }
				
				logLikelihood += Dirichlet.logGamma(beta + topicCounts[topic]) -
					logGammaBeta;

				if (Double.isNaN(logLikelihood)) {
					System.out.println(topicCounts[topic]);
					System.exit(1);
				}
			}
		}
	
		for (int topic=0; topic < numTopics; topic++) {
			logLikelihood -= 
				Dirichlet.logGamma( (beta * numTypes) +
											tokensPerTopic[ topic ] );
			if (Double.isNaN(logLikelihood)) {
				System.out.println("after topic " + topic + " " + tokensPerTopic[ topic ]);
				System.exit(1);
			}

		}
	
		logLikelihood += 
			numTopics * Dirichlet.logGamma(beta * numTypes);

		if (Double.isNaN(logLikelihood)) {
			System.out.println("at the end");
			System.exit(1);
		}


		return logLikelihood;
	}

	// 
	// Methods for displaying and saving results
	//

	public String topWords (int numWords) {

		StringBuilder output = new StringBuilder();

		IDSorter[] sortedWords = new IDSorter[numTypes];

		for (int topic = 0; topic < numTopics; topic++) {
			for (int type = 0; type < numTypes; type++) {
				sortedWords[type] = new IDSorter(type, typeTopicCounts[type][topic]);
			}

			Arrays.sort(sortedWords);
			
			output.append(topic + "\t" + tokensPerTopic[topic] + "\t");
			for (int i=0; i < numWords; i++) {
				output.append(alphabet.lookupObject(sortedWords[i].getID()) + " ");
			}
			output.append("\n");
		}

		return output.toString();
	}

	/**
	 *  @param file        The filename to print to
	 *  @param threshold   Only print topics with proportion greater than this number
	 *  @param max         Print no more than this many topics
	 */
	public void printDocumentTopics (File file, double threshold, int max) throws IOException {
		PrintWriter out = new PrintWriter(file);

		out.print ("#doc source topic proportion ...\n");
		int docLen;
		int[] topicCounts = new int[ numTopics ];

		IDSorter[] sortedTopics = new IDSorter[ numTopics ];
		for (int topic = 0; topic < numTopics; topic++) {
			// Initialize the sorters with dummy values
			sortedTopics[topic] = new IDSorter(topic, topic);
		}

		if (max < 0 || max > numTopics) {
			max = numTopics;
		}

		for (int doc = 0; doc < data.size(); doc++) {
			LabelSequence topicSequence = (LabelSequence) data.get(doc).topicSequence;
			int[] currentDocTopics = topicSequence.getFeatures();

			out.print (doc); out.print (' ');

			if (data.get(doc).instance.getSource() != null) {
				out.print (data.get(doc).instance.getSource()); 
			}
			else {
				out.print ("null-source");
			}

			out.print (' ');
			docLen = currentDocTopics.length;

			// Count up the tokens
			for (int token=0; token < docLen; token++) {
				topicCounts[ currentDocTopics[token] ]++;
			}

			// And normalize
			for (int topic = 0; topic < numTopics; topic++) {
				sortedTopics[topic].set(topic, (float) topicCounts[topic] / docLen);
			}
			
			Arrays.sort(sortedTopics);

			for (int i = 0; i < max; i++) {
				if (sortedTopics[i].getWeight() < threshold) { break; }
				
				out.print (sortedTopics[i].getID() + " " + 
						  sortedTopics[i].getWeight() + " ");
			}
			out.print (" \n");

			Arrays.fill(topicCounts, 0);
		}
		
	}
	
	public void printState (File f) throws IOException {
		PrintStream out =
			new PrintStream(new GZIPOutputStream(new BufferedOutputStream(new FileOutputStream(f))));
		printState(out);
		out.close();
	}
	
	public void printState (PrintStream out) {

		out.println ("#doc source pos typeindex type topic");

		for (int doc = 0; doc < data.size(); doc++) {
			FeatureSequence tokenSequence =	(FeatureSequence) data.get(doc).instance.getData();
			LabelSequence topicSequence =	(LabelSequence) data.get(doc).topicSequence;

			String source = "NA";
			if (data.get(doc).instance.getSource() != null) {
				source = data.get(doc).instance.getSource().toString();
			}

			for (int position = 0; position < topicSequence.getLength(); position++) {
				int type = tokenSequence.getIndexAtPosition(position);
				int topic = topicSequence.getIndexAtPosition(position);
				out.print(doc); out.print(' ');
				out.print(source); out.print(' '); 
				out.print(position); out.print(' ');
				out.print(type); out.print(' ');
				out.print(alphabet.lookupObject(type)); out.print(' ');
				out.print(topic); out.println();
			}
		}
	}
	
	
	// Serialization
	
	private static final long serialVersionUID = 1;
	private static final int CURRENT_SERIAL_VERSION = 0;
	private static final int NULL_INTEGER = -1;
	
	public void write (File f) {
		try {
			ObjectOutputStream oos = new ObjectOutputStream (new FileOutputStream(f));
			oos.writeObject(this);
			oos.close();
		}
		catch (IOException e) {
			System.err.println("Exception writing file " + f + ": " + e);
		}
	}
	
	private void writeObject (ObjectOutputStream out) throws IOException {
		out.writeInt (CURRENT_SERIAL_VERSION);

		// Instance lists
		out.writeObject (data);
		out.writeObject (alphabet);
		out.writeObject (topicAlphabet);

		out.writeInt (numTopics);
		out.writeObject (alpha);
		out.writeDouble (beta);
		out.writeDouble (betaSum);

		out.writeInt(showTopicsInterval);
		out.writeInt(wordsPerTopic);

		out.writeObject(random);
		out.writeObject(formatter);
		out.writeBoolean(printLogLikelihood);

		out.writeObject (typeTopicCounts);

		for (int ti = 0; ti < numTopics; ti++) {
			out.writeInt (tokensPerTopic[ti]);
		}
	}
	
	private void readObject (ObjectInputStream in) throws IOException, ClassNotFoundException {
		int featuresLength;
		int version = in.readInt ();

		data = (ArrayList<TopicAssignment>) in.readObject ();
		alphabet = (Alphabet) in.readObject();
		topicAlphabet = (LabelAlphabet) in.readObject();

		numTopics = in.readInt();
		alpha = in.readDouble();
		alphaSum = alpha * numTopics;
		beta = in.readDouble();
		betaSum = in.readDouble();

		showTopicsInterval = in.readInt();
		wordsPerTopic = in.readInt();

		random = (Randoms) in.readObject();
		formatter = (NumberFormat) in.readObject();
		printLogLikelihood = in.readBoolean();
		
		int numDocs = data.size();
		this.numTypes = alphabet.size();

		typeTopicCounts = (int[][]) in.readObject();
		tokensPerTopic = new int[numTopics];
		for (int ti = 0; ti < numTopics; ti++) {
			tokensPerTopic[ti] = in.readInt();
		}
	}

	public static void main (String[] args) throws IOException {

		InstanceList training = InstanceList.load (new File(args[0]));

		int numTopics = args.length > 1 ? Integer.parseInt(args[1]) : 200;

		CustomLDA lda = new CustomLDA (numTopics, 50.0, 0.01);
		lda.addInstances(training);
		lda.sample(1000);
	}
	
	// Custom methods: Author- Sumanta
	public TopicInferencer getInferencer(){
		double[] alphaVals = new double[this.numTopics];
		for(int i=0; i<this.numTopics; i++)
			alphaVals[i] = this.alpha;
		TopicInferencer inf = new TopicInferencer(this.typeTopicCounts, this.tokensPerTopic, this.alphabet, alphaVals, this.beta, this.betaSum);
		return inf;
	}
	public double[] getTopicProbabilities(int instanceIndex){
		double[] topicDist = new double[this.numTopics];
		TopicAssignment currTopicAssign = this.data.get(instanceIndex);
		LabelSequence topicSeqForCurr = currTopicAssign.topicSequence;
		for(int tokenPos=0; tokenPos<topicSeqForCurr.size(); tokenPos++){
			int currTopic = topicSeqForCurr.getIndexAtPosition(tokenPos);
			topicDist[currTopic]++;
		}
		
		// Add the smoothing parameters and normalize
		// Also here we are using same alpha for each topic
		double sum = 0.0;
		for (int topic = 0; topic < numTopics; topic++) {
			topicDist[topic] += this.alpha;
			sum += topicDist[topic];
		}

		// And normalize
		for (int topic = 0; topic < numTopics; topic++) {
			topicDist[topic] /= sum;
		}
		return topicDist;
	}
	// addInstances method for unigram style //
	public void addInstancesUnigram (InstanceList training) {

		alphabet = training.getDataAlphabet();
		numTypes = alphabet.size();
		
		if(beta<0)
			beta = CustomLDA.BETA_SUM/numTypes;
		betaSum = beta * numTypes;

		typeTopicCounts = new int[numTypes][numTopics];
		initTypeTopicCounts();
		topicDocSequence = new int[training.size()];

		int doc = 0;

		for (Instance instance : training) {
			FeatureSequence tokens = (FeatureSequence) instance.getData();
			LabelSequence topicSequence =
					new LabelSequence(topicAlphabet, new int[ tokens.size() ]);

			int[] topics = topicSequence.getFeatures();
			int topic = random.nextInt(numTopics);
			topicDocSequence[doc] = topic;
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
	public void sampleUMM (int iterations) throws IOException {

		for (int iteration = 1; iteration <= iterations; iteration++) {

			long iterationStart = System.currentTimeMillis();
			if(printMessages){
				System.out.print("Iteration "+iteration+", topicDocSeq array=[");
				for(int doci=0; doci<topicDocSequence.length; doci++)
					System.out.print(topicDocSequence[doci]+" ");
				System.out.print("]\n");
			}
			// Loop over every document in the corpus
			for (int doc = 0; doc < data.size(); doc++) {
				FeatureSequence tokenSequence =
					(FeatureSequence) data.get(doc).instance.getData();
				LabelSequence topicSequence =
					(LabelSequence) data.get(doc).topicSequence;

				sampleTopicsForOneDocUnigram (tokenSequence, topicSequence, doc);
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

			// Occasionally print more information
			if (showTopicsInterval != 0 && iteration % showTopicsInterval == 0) {
				logger.info("<" + iteration + "> Log Likelihood: " + modelLogLikelihood() + "\n" +
							topWords (wordsPerTopic));
			}

		}
	}
	private void initTypeTopicCounts(){
		for(int type=0; type<numTypes; type++){
			for(int topic=0; topic<numTopics; topic++)
				typeTopicCounts[type][topic] = topic;
		}
	}
	protected void sampleTopicsForOneDocUnigram (FeatureSequence tokenSequence, LabelSequence topicSequence,
			  int docIndex) {
		int oldDocTopic = topicDocSequence[docIndex];
		
		int newTopic;
		int docLength = tokenSequence.getLength(); // T
		
		for(int position=0; position<docLength; position++){
			int currType = tokenSequence.getIndexAtPosition(position);
			typeTopicCounts[currType][oldDocTopic]-=topicMask+1;
		}
		// now typeTopicCounts is n_{k,-i}^(t)
		tokensPerTopic[oldDocTopic]-=docLength;
		// now tokensPerTopic is tc_{k,-i}
		topicDocSequence[docIndex] = -1;
		
		double[] topicScores = new double[numTopics];
		double[] smoothedScores = new double[numTopics];
		double[] logpkList = new double[numTopics];
		double topicScoreSum = 0, normalizingC, currentMin, currentMax, minLog, maxLog;
		
		//-- p'(k) distribution is calculated in log space with normalization --//
		
		if(printMessages)
			System.out.println("beta="+beta+", betaSum="+betaSum+", V="+numTypes);
		for(int k=0; k<numTopics; k++){
			double sumOfLog = 0, logOfpk;
			for(int t=0; t<docLength; t++)
				sumOfLog+=Math.log((typeTopicCounts[t][k]>>topicBits)+beta);
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
		double sumSmooth = numTopics-(numTopics-1)*CustomLDA.LAMBDA;
		for(int k=0; k<numTopics; k++){
			smoothedScores[k] = (CustomLDA.LAMBDA*topicScores[k]+(1-CustomLDA.LAMBDA))/sumSmooth;
			if(printMessages)
				System.out.println("smoothedScore("+k+")="+smoothedScores[k]);
		}
		//double sample = random.nextUniform()*topicScoreSum;
		double sample = random.nextUniform();
		if(printMessages)
			System.out.println("Sample="+sample);
		newTopic = -1;
		while (sample > 0.0) {
			newTopic++;
			sample -= smoothedScores[newTopic];
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
		topicDocSequence[docIndex] = newTopic;
	}
}
