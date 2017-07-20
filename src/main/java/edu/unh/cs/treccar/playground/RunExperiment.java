package edu.unh.cs.treccar.playground;

import java.util.ArrayList;

import com.trolltech.qt.gui.*;

public class RunExperiment {
	// k0beta_umm
	public static final boolean SAVE_RESULT = false;
	public static final String CLUSTERING_MEASURE_FILENAME = "umm_iter";
	public static final String TRECEVAL_ASSIGN_FILENAME = "umm_iter_trec";
	public static final boolean RUN_BY_PAGE = true;
	//public static final double LAMBDA_UMM = 0.99;
	public static final boolean SMOOTHED_UMM = true;
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		/*
		QApplication.initialize(args);

        QPushButton hello = new QPushButton("Hello World!");
        hello.resize(120, 40);
        hello.setWindowTitle("Hello World");
        hello.show();

        QApplication.execStatic();
        */
		/*
        SingleRun sr = new SingleRun(new String[]{"0","20","1","0","1.0","0.3",
        		"/home/sumanta/Documents/new_research/unh/test200-v1.4/all.test200.cbor.paragraphs",
        		"/home/sumanta/Documents/new_research/unh/test200-v1.4/all.test200.cbor.outlines",
        		"/home/sumanta/Documents/new_research/unh/test200-v1.4/all.test200.cbor.toplevel.qrels",
        		"/home/sumanta/Documents/new_research/unh/test200-v1.4/all.test200.cbor.article.qrels",
        		""});
        sr.runExperiment();
        */
        String model = "1";
        // 1 for LDA, 2 for kmeans, 3 for unigram tm
        String tw = "0"; // currently ignored
        String startK="0", startIter="300", startAlpha="1.0", startBeta="260"; 
        // treat alpha and beta as alphaSum and betaSum
        String stopK=startK, stopIter=startIter, stopAlpha=startAlpha, stopBeta=startBeta;
        String stepK="1", stepIter="1", stepAlpha="1", stepBeta="1";
        boolean[] isVar = {false, false, false, false};
        if(isVar[0]){
        	stopK = "20";
        	stepK = "1";
        }
        if(isVar[1]){
        	stopIter = "30";
        	stepIter = "1";
        }
        if(isVar[2]){
        	stopAlpha = "5.0";
        	stepAlpha = "0.05";
        }
        if(isVar[3]){
        	stopBeta = "150";
        	stepBeta = "25";
        }
        for(int a=Integer.parseInt(startK); !(a>Integer.parseInt(stopK)); a+=Integer.parseInt(stepK)){
        	for(int b=Integer.parseInt(startIter); !(b>Integer.parseInt(stopIter)); b+=Integer.parseInt(stepIter)){
        		for(double c=Double.parseDouble(startAlpha); !(c>Double.parseDouble(stopAlpha)); c+=Double.parseDouble(stepAlpha)){
        			for(double d=Double.parseDouble(startBeta); !(d>Double.parseDouble(stopBeta)); d+=Double.parseDouble(stepBeta)){
        				System.out.println("K="+a+",Iter="+b+",AlphaSum="+c+",Beta="+d);
        				SingleRun sr = new SingleRun(new String[]{Integer.toString(a),Integer.toString(b),model,tw,Double.toString(c),Double.toString(d),
        		        		"/home/sumanta/Documents/new_research/unh/test200-v1.4/all.test200.cbor.paragraphs",
        		        		"/home/sumanta/Documents/new_research/unh/test200-v1.4/all.test200.cbor.outlines",
        		        		"/home/sumanta/Documents/new_research/unh/test200-v1.4/all.test200.cbor.toplevel.qrels",
        		        		"/home/sumanta/Documents/new_research/unh/test200-v1.4/all.test200.cbor.article.qrels",
        		        		"/home/sumanta/Documents/new_research/unh/test200-v1.4results/results4paper"});
        		        if(RunExperiment.RUN_BY_PAGE)
        		        	sr.runExperiment();
        		        else
        		        	sr.runExperimentWholeCorpus();
        			}
        		}
        	}
        }
        
	}

}
