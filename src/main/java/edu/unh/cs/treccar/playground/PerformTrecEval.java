package edu.unh.cs.treccar.playground;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashSet;

import edu.unh.cs.treccar.Data;
import edu.unh.cs.treccar.read_data.DeserializeData;
import edu.unh.cs.treccar.read_data.DeserializeData.RuntimeCborException;

public class PerformTrecEval {
	public static final String TRECEVAL_DIR = "/home/sumanta/Documents/trec_eval.9.0";
	public static final String OUTLINE = "/home/sumanta/Documents/new_research/unh/test200-v1.4/all.test200.cbor.outlines";
	//public static final String GT_PATH = "/home/sumanta/Documents/new_research/unh/test200-v1.4/all.test200.cbor.article.qrels";
	public static final String GT_PATH = "/home/sumanta/Documents/new_research/unh/test200-v1.4/all.test200.cbor.toplevel.qrels";
	public static final boolean BY_PAGE = true;
	
	public static void main(String[] args){
		String workDir = "/home/sumanta/Documents/new_research/unh/test200-v1.4results/detective_results2";
		String assign_filename = "random_trec";
		PerformTrecEval pte = new PerformTrecEval();
		HashSet<String> runids = pte.getRunIDs(workDir+"/"+assign_filename);
		String tempFilepath = workDir+"/temp";
		String tempGTpath = workDir+"/tempgt";
		for(String runid:runids){
			if(BY_PAGE){
				ArrayList<Data.Page> pageList = pte.getPageListFromPath(OUTLINE);
				System.out.println();
				for(Data.Page page:pageList){
					pte.prepareGTForPage(GT_PATH, tempGTpath, page);
					pte.prepareFileForPage(runid, tempFilepath, workDir+"/"+assign_filename, page);
					System.out.print(page.getPageId().toString()+" ");
					pte.callTrecEvalScript(tempFilepath, tempGTpath);
				}
			}
			else{
				pte.prepareFileForSingleRunid(runid, tempFilepath, workDir+"/"+assign_filename);
				pte.callTrecEvalScript(tempFilepath, GT_PATH);
			}
		}
	}
	public void callTrecEvalScript(String outfilePath, String gtPath){
		try {
			Process pr = new ProcessBuilder(TRECEVAL_DIR+"/trec_eval", gtPath, outfilePath).start();
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
	public ArrayList<Data.Page> getPageListFromPath(String path){
		ArrayList<Data.Page> pageList = new ArrayList<Data.Page>();
		try {
			FileInputStream fis = new FileInputStream(new File(path));
			for(Data.Page page: DeserializeData.iterableAnnotations(fis))
				pageList.add(page);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (RuntimeCborException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return pageList;
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
	public void prepareFileForPage(String runid, String tempFilePath, String assignFilePath, Data.Page page){
		try {
			FileWriter tempFileWriter = new FileWriter(tempFilePath, false);
			BufferedReader br = new BufferedReader(new FileReader(assignFilePath));
			String line = br.readLine();
			while(line!=null){
				String currRunid = line.split(" ")[5];
				if(currRunid.equals(runid) && line.startsWith(page.getPageId().toString())){
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
	public void prepareGTForPage(String masterGT, String tempGT, Data.Page page){
		try {
			FileWriter tempGTWriter = new FileWriter(tempGT, false);
			BufferedReader br = new BufferedReader(new FileReader(masterGT));
			String line = br.readLine();
			while(line!=null){
				if(line.startsWith(page.getPageId().toString())){
					tempGTWriter.write(line+"\n");
				}
				line = br.readLine();
			}
			br.close();
			tempGTWriter.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
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
