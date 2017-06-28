package edu.unh.cs.treccar.playground;

import java.util.ArrayList;

import com.trolltech.qt.gui.*;

public class RunExperiment {
	
	public static final String RAND_RESULT_FILENAME = "rand_k_umm";
	public static final String TRECEVAL_ASSIGN_FILENAME = "trec_k_umm";
	
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
        String tw = "0"; // currently ignored
        String startK="0", startIter="30", startAlpha="1.3", startBeta="0.3"; // treat Aplpha as AplphaSum
        String stopK=startK, stopIter=startIter, stopAlpha=startAlpha, stopBeta=startBeta;
        String stepK="1", stepIter="1", stepAlpha="1", stepBeta="1";
        boolean[] isVar = {true, false, false, false};
        if(isVar[0]){
        	stopK = "21";
        	stepK = "1";
        }
        if(isVar[1]){
        	stopIter = "100";
        	stepIter = "10";
        }
        if(isVar[2]){
        	stopAlpha = "4.5";
        	stepAlpha = "0.4";
        }
        if(isVar[3]){
        	stopBeta = "1";
        	stepBeta = "0.1";
        }
        for(int a=Integer.parseInt(startK); !(a>Integer.parseInt(stopK)); a+=Integer.parseInt(stepK)){
        	for(int b=Integer.parseInt(startIter); !(b>Integer.parseInt(stopIter)); b+=Integer.parseInt(stepIter)){
        		for(double c=Double.parseDouble(startAlpha); !(c>Double.parseDouble(stopAlpha)); c+=Double.parseDouble(stepAlpha)){
        			for(double d=Double.parseDouble(startBeta); !(d>Double.parseDouble(stopBeta)); d+=Double.parseDouble(stepBeta)){
        				System.out.println("K="+a+",Iter="+b+",Alpha="+c+",Beta="+d);
        				SingleRun sr = new SingleRun(new String[]{Integer.toString(a),Integer.toString(b),model,tw,Double.toString(c),Double.toString(d),
        		        		"/home/sumanta/Documents/new_research/unh/test200-v1.4/all.test200.cbor.paragraphs",
        		        		"/home/sumanta/Documents/new_research/unh/test200-v1.4/all.test200.cbor.outlines",
        		        		"/home/sumanta/Documents/new_research/unh/test200-v1.4/all.test200.cbor.toplevel.qrels",
        		        		"/home/sumanta/Documents/new_research/unh/test200-v1.4/all.test200.cbor.article.qrels",
        		        		"/home/sumanta/Documents/new_research/unh/test200-v1.4results/custom_lda_and_km_results"});
        		        sr.runExperiment();
        			}
        		}
        	}
        }
        
	}

}
