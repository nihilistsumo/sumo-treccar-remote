package edu.unh.cs.treccar.playground;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;

import org.apache.commons.lang3.builder.CompareToBuilder;

import co.nstant.in.cbor.CborException;
import edu.unh.cs.treccar.Data;
import edu.unh.cs.treccar.read_data.DeserializeData;

public class MeasurePerformanceMain {
	//Arguments: ground truth file, outline file, paragraph file
	public static final int TOPICNUMBER=100;
	public static final int CLUSTERNUMBER=100;
	public static final int ITERATION=2;
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		try {			
			FileOutputStream fos = new FileOutputStream(new File(args[3]), true);
			PrintStream ps = new PrintStream(fos);
			ps.println("Topic Number: "+TOPICNUMBER+", Cluster Number: "+CLUSTERNUMBER);
			ps.println("############################################");
			ps.println();
			for(int i=0; i<ITERATION; i++){
				MeasurePerformance mp = new MeasurePerformance(args[0], args[1], args[2]);
				BaseClusterer base = new BaseClusterer(args[2]);
				PerformanceMetrics pm = new PerformanceMetrics(mp.getGroundTruth(), mp.getCandidateAssign());
				runModel(base, pm, ps);
				runClustering(base, pm, ps);
				//runBaseline(base, pm);
			}
			ps.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}	
	}
	private static void runModel(BaseClusterer base, PerformanceMetrics pm, PrintStream ps){
		ps.println("Topic model numbers");
		ps.println(pm.analyzeBasic());
		ps.println(pm.calculateRandIndex());
		ps.println(pm.calculatePurity());
		ps.println("-------------------\n");
	}
	private static void runClustering(BaseClusterer base, PerformanceMetrics pm, PrintStream ps){
		ps.println("Clustering numbers");
		ps.println(pm.calculateRandIndex(base.getClusters(base.clusterNum)));
		ps.println(pm.calculatePurity(base.getClusters(base.clusterNum)));
		ps.println("-------------------\n");
	}
	private static void runBaseline(BaseClusterer base, PerformanceMetrics pm){
		System.out.println("cluster num: "+base.clusterNum);
		System.out.println("All in one:");
		System.out.println(pm.calculateRandIndex(base.getBaselineAllInOne(base.clusterNum)));
		System.out.println(pm.calculatePurity(base.getBaselineAllInOne(base.clusterNum)));
		System.out.println("All correct:");
		System.out.println(pm.calculateRandIndex(base.getBaselineAllCorrect()));
		System.out.println(pm.calculatePurity(base.getBaselineAllCorrect()));
		System.out.println("Random:");
		System.out.println(pm.calculateRandIndex(base.getBaselineRandom(base.clusterNum)));
		System.out.println(pm.calculatePurity(base.getBaselineRandom(base.clusterNum)));
		System.out.println();
	}
}
