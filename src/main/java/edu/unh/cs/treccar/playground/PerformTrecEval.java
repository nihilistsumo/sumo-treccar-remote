package edu.unh.cs.treccar.playground;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashSet;

public class PerformTrecEval {
	public static final String TRECEVAL_DIR = "/home/sumanta/Documents/trec_eval.9.0";
	//public static final String GT_PATH = "/home/sumanta/Documents/new_research/unh/test200-v1.4/all.test200.cbor.article.qrels";
	public static final String GT_PATH = "/home/sumanta/Documents/new_research/unh/test200-v1.4/all.test200.cbor.toplevel.qrels";
	
	public static void main(String[] args){
		String workDir = "/home/sumanta/Documents/new_research/unh/test200-v1.4results/results4paper/k";
		String assign_filename = "klda_new_trec";
		PerformTrecEval pte = new PerformTrecEval();
		HashSet<String> runids = pte.getRunIDs(workDir+"/"+assign_filename);
		String tempFilepath = workDir+"/temp";
		for(String runid:runids){
			pte.prepareFileForSingleRunid(runid, tempFilepath, workDir+"/"+assign_filename);
			try {
				Process pr = new ProcessBuilder(PerformTrecEval.TRECEVAL_DIR+"/trec_eval", PerformTrecEval.GT_PATH, tempFilepath).start();
				InputStream is = pr.getInputStream();
				InputStream erris = pr.getErrorStream();
				InputStreamReader isr = new InputStreamReader(is);
				InputStreamReader errisr = new InputStreamReader(erris);
				BufferedReader br = new BufferedReader(isr);
				BufferedReader errbr = new BufferedReader(errisr);
				String line, errline;
				while ((line = br.readLine()) != null) {
					System.out.print(line.split("\t")[2]+" ");
				}
				System.out.println();
				while ((errline = errbr.readLine()) != null) {
					System.out.println("Error message from script: "+errline);
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	public HashSet<String> getRunIDs(String assignFilepath){
		HashSet<String> runids = new HashSet<String>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(assignFilepath));
			String line = br.readLine();
			while(line!=null){
				String currRunid = line.split(" ")[5];
				runids.add(currRunid);
				line = br.readLine();
			}
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return runids;
	}
	public void prepareFileForSingleRunid(String runid, String tempFilePath, String assignFilePath){
		try {
			FileWriter tempFileWriter = new FileWriter(tempFilePath, false);
			BufferedReader br = new BufferedReader(new FileReader(assignFilePath));
			String line = br.readLine();
			while(line!=null){
				String currRunid = line.split(" ")[5];
				if(currRunid.equals(runid)){
					tempFileWriter.write(line+"\n");
				}
				line = br.readLine();
			}
			br.close();
			tempFileWriter.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
