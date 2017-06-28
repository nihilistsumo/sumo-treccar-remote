package edu.unh.cs.treccar.playground.split;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

import co.nstant.in.cbor.CborException;
import edu.unh.cs.treccar.playground.MeasurePerformance;

public class SplitDataset {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String workingFolder = args[0];
		String originalOutline = args[1];
		String originalPara = args[2];
		final String gtFilename = "all.test200.cbor.hierarchical.qrels";
		try {
			SplitOutline so = new SplitOutline(workingFolder, gtFilename, originalOutline);
			SplitParagraph sp = new SplitParagraph(workingFolder, gtFilename, originalPara);
			so.split();
			sp.split();
		}
		catch (IOException | CborException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	public static HashSet<String> getPageidsFromGT(String gtFile) throws IOException{
		ArrayList<String> pageIds = new ArrayList<String>();
		BufferedReader br = new BufferedReader(new FileReader(gtFile));
		String line;
		while((line = br.readLine())!=null){
			pageIds.add(line.split("[' ''/']")[0]);
		}
		HashSet<String> pageIdSet = new HashSet<String>(pageIds);
		br.close();
		return pageIdSet;
	}
	public static HashSet<String> getParaidsFromGT(String gtFile) throws IOException{
		ArrayList<String> paraIds = new ArrayList<String>();
		BufferedReader br = new BufferedReader(new FileReader(gtFile));
		String line;
		while((line = br.readLine())!=null){
			paraIds.add(line.split(" ")[2]);
		}
		HashSet<String> paraIdSet = new HashSet<>(paraIds);
		br.close();
		return paraIdSet;
	}
}
