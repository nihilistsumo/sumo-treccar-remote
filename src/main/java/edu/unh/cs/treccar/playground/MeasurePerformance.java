package edu.unh.cs.treccar.playground;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;

import edu.unh.cs.treccar.Data;

public class MeasurePerformance {
	static String GTFILE, OUTLINEFILE, PARAFILE;
	// In groundTruth we map sectionID to list of paraIDs from the ground truth file
	private HashMap<String, ArrayList<String>> groundTruth = new HashMap<String, ArrayList<String>>();
	private HashMap<AssignParagraphs.SectionPathID, ArrayList<Data.Paragraph>> candidateAssign = new HashMap<AssignParagraphs.SectionPathID, ArrayList<Data.Paragraph>>();
	public HashMap<String, ArrayList<String>> getGroundTruth() {
		return groundTruth;
	}
	private void setGroundTruth() {
		HashMap<String, ArrayList<String>> localGT = new HashMap<String, ArrayList<String>>();
		BufferedReader br;
		try{
			br = new BufferedReader(new FileReader(MeasurePerformance.GTFILE));
			String line;
			String[] lineData = new String[4];
			while((line = br.readLine()) != null){
				lineData = line.split(" ");
				if(localGT.containsKey(lineData[0])){
					localGT.get(lineData[0]).add(lineData[2]);
				} else{
					ArrayList<String> paraList = new ArrayList<String>();
					paraList.add(lineData[2]);
					localGT.put(lineData[0], paraList);
				}
				
			}
			br.close();
		} catch (Exception e){
			e.printStackTrace();
		}
		this.groundTruth = localGT;
	}
	public HashMap<AssignParagraphs.SectionPathID, ArrayList<Data.Paragraph>> getCandidateAssign() {
		return candidateAssign;
	}
	private void setCandidateAssign() {
		HashMap<AssignParagraphs.SectionPathID, ArrayList<Data.Paragraph>> candidateAssign = new HashMap<AssignParagraphs.SectionPathID, ArrayList<Data.Paragraph>>();
		candidateAssign = (new AssignParagraphs()).assignCandidates(MeasurePerformance.OUTLINEFILE, MeasurePerformance.PARAFILE);
		this.candidateAssign = candidateAssign;
	}
	public MeasurePerformance(String gtFilepath, String outlinePath, String paraPath){
		MeasurePerformance.GTFILE = gtFilepath;
		MeasurePerformance.OUTLINEFILE = outlinePath;
		MeasurePerformance.PARAFILE = paraPath;
		this.setGroundTruth();
		this.setCandidateAssign();
	}
	 
}
