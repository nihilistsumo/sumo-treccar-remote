package edu.unh.cs.treccar.playground;

import java.util.ArrayList;

import com.trolltech.qt.gui.*;

public class RunExperiment {
	// k0beta_umm
	public static final boolean SAVE_RESULT = true;
	public static final String CLUSTERING_MEASURE_FILENAME = "chi_km_tds";
	public static final String TRECEVAL_ASSIGN_FILENAME = "chi_km_tds_trec";
	public static final boolean RUN_BY_PAGE = true;
	public static final boolean ASSIGN_RANK_MODE = true;
	public static final boolean SMOOTHED_UMM = false;
	public static final int ASSIGN_METHOD = 4; // 1- KLDiv, 2- KS, 3- Bhat, 4- Chi
	
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
        String model = "2";
        // 1 for LDA, 2 for kmeans, 3 for unigram tm
        String tw = "0"; // 0- tf (default), 1- td, 2- tds
        String startK="0", startIter="100", startAlpha="1.0", startBeta="260"; 
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
        	stopBeta = "560";
        	stepBeta = "50";
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
        		        		"/home/sumanta/Documents/new_research/unh/test200-v1.4results/detective_results"});
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
