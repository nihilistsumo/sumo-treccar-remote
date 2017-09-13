package edu.unh.cs.treccar.playground;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Pattern;

import cc.mallet.pipe.CharSequence2TokenSequence;
import cc.mallet.pipe.FeatureSequence2FeatureVector;
import cc.mallet.pipe.Input2CharSequence;
import cc.mallet.pipe.Pipe;
import cc.mallet.pipe.SerialPipes;
import cc.mallet.pipe.TokenSequence2FeatureSequence;
import cc.mallet.pipe.TokenSequenceLowercase;
import cc.mallet.pipe.TokenSequenceRemoveStopwords;
import cc.mallet.types.FeatureVector;
import cc.mallet.types.Instance;
import cc.mallet.types.InstanceList;
import co.nstant.in.cbor.CborException;
import edu.unh.cs.treccar.Data;
import edu.unh.cs.treccar.read_data.DeserializeData;
import edu.unh.cs.treccar.read_data.DeserializeData.RuntimeCborException;

public class ExaminePages {
	static final String TEST200DIR = "/home/sumanta/Documents/new_research/unh/test200-v1.4";
	static final String OUTLINEFILE = "all.test200.cbor.outlines";
	static final String PARAFILE = "all.test200.cbor.paragraphs";
	static final String ARTICLEQRELS = "all.test200.cbor.article.qrels";
	static final String OUTPUT_DIR = "/home/sumanta/Documents/new_research/unh/test200-v1.4results";
	static final String OUTPUT_TEXT = "out_text";
	static final String OUTPUT_VECTOR = "out_vec";
	
	public ArrayList<Data.Page> pagelist;
	public ArrayList<Data.Paragraph> paralist;
	public HashMap<Data.Page, ArrayList<Data.Paragraph>> pageParaMap;
	
	public ExaminePages(){
		this.pagelist = this.getPageListFromPath(TEST200DIR+"/"+OUTLINEFILE);
		this.paralist = this.getParaListFromPath(TEST200DIR+"/"+PARAFILE);
		this.pageParaMap = this.getPageParaObjMap(this.getGTMapFromPath(TEST200DIR+"/"+ARTICLEQRELS));
	}
	
	// ------Existing methods from SingleRun------- //
	private ArrayList<Data.Paragraph> getParaListFromPath(String path){
		ArrayList<Data.Paragraph> paraList = new ArrayList<Data.Paragraph>();
		try {
			FileInputStream fis = new FileInputStream(new File(path));
			for(Data.Paragraph para:DeserializeData.iterableParagraphs(fis))
				paraList.add(para);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (CborException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return paraList;
	}
	private ArrayList<Data.Page> getPageListFromPath(String path){
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
	private HashMap<String, ArrayList<String>> getGTMapFromPath(String path){
		HashMap<String, ArrayList<String>> gtMap = new HashMap<String, ArrayList<String>>();
		BufferedReader br;
		try{
			br = new BufferedReader(new FileReader(path));
			String line;
			String[] lineData = new String[4];
			while((line = br.readLine()) != null){
				lineData = line.split(" ");
				if(gtMap.containsKey(lineData[0])){
					gtMap.get(lineData[0]).add(lineData[2]);
				} else{
					ArrayList<String> paraList = new ArrayList<String>();
					paraList.add(lineData[2]);
					gtMap.put(lineData[0], paraList);
				}	
			}
			br.close();
		} catch (Exception e){
			e.printStackTrace();
		}
		return gtMap;
	}
	private InstanceList convertParasToIList(ArrayList<Data.Paragraph> paraObjs){
		InstanceList iListPara = new InstanceList(ExaminePages.buildPipe());
		for(Data.Paragraph paraObj:paraObjs){
			Instance paraIns = new Instance(paraObj.getTextOnly(), null, paraObj.getParaId(), paraObj.getTextOnly());
			iListPara.addThruPipe(paraIns);
		}
		return iListPara;
	}
	// ----------------------------------------- //
	
	private HashMap<Data.Page, ArrayList<Data.Paragraph>> getPageParaObjMap(
			HashMap<String, ArrayList<String>> pageParaMap){
		HashMap<Data.Page, ArrayList<Data.Paragraph>> objMap = new HashMap<Data.Page, ArrayList<Data.Paragraph>>();
		HashMap<String, ArrayList<String>> idMap = this.getGTMapFromPath(ExaminePages.TEST200DIR+"/"+ExaminePages.ARTICLEQRELS);
		String pageid;
		for(Data.Page page:this.pagelist){
			ArrayList<Data.Paragraph> paralistForPage = new ArrayList<Data.Paragraph>();
			pageid = page.getPageId();
			for(Data.Paragraph para:this.paralist){
				if(idMap.get(pageid).contains(para.getParaId()))
					paralistForPage.add(para);
			}
			objMap.put(page, paralistForPage);
		}
		return objMap;
	}
	private void printBasicAnalysis() throws IOException{
		ArrayList<Data.Paragraph> parasInPage;
		InstanceList iListPara;
		FileWriter fwt = new FileWriter(new File(OUTPUT_DIR+"/"+OUTPUT_TEXT), false);
		FileWriter fwv = new FileWriter(new File(OUTPUT_DIR+"/"+OUTPUT_VECTOR), false);
		for(Data.Page page:this.pageParaMap.keySet()){
			parasInPage = this.pageParaMap.get(page);
			iListPara = this.convertParasToIList(parasInPage);
			fwt.write(page.getPageId()+" has "+parasInPage.size()+" paras\n");
			fwv.write(page.getPageId()+" has "+parasInPage.size()+" paras\n");
			for(Data.Paragraph para:parasInPage){
				fwt.write(para.getParaId()+"\n");
				fwt.write(para.getTextOnly()+"\n");
			}
			for(Instance paraIns:iListPara){
				fwv.write(paraIns.getName().toString()+"\n");
				fwv.write(((FeatureVector)paraIns.getData()).toString(true)+"\n");
			}
			fwt.write("\n");
			fwv.write("\n");
		}
		fwt.close();
		fwv.close();
	}
	public static Pipe buildPipe(){
		ArrayList pipeList = new ArrayList();

        // Read data from File objects
        pipeList.add(new Input2CharSequence("UTF-8"));

        // Regular expression for what constitutes a token.
        //  This pattern includes Unicode letters, Unicode numbers, 
        //   and the underscore character. Alternatives:
        //    "\\S+"   (anything not whitespace)
        //    "\\w+"    ( A-Z, a-z, 0-9, _ )
        //    "[\\p{L}\\p{N}_]+|[\\p{P}]+"   (a group of only letters and numbers OR
        //                                    a group of only punctuation marks)
        Pattern tokenPattern = Pattern.compile("[\\p{L}\\p{N}_]+");

        // Tokenize raw strings
        pipeList.add(new CharSequence2TokenSequence(tokenPattern));

        // Normalize all tokens to all lowercase
        pipeList.add(new TokenSequenceLowercase());

        // Remove stopwords from a standard English stoplist.
        //  options: [case sensitive] [mark deletions]
        pipeList.add(new TokenSequenceRemoveStopwords(false, false));

        // Rather than storing tokens as strings, convert 
        //  them to integers by looking them up in an alphabet.
        pipeList.add(new TokenSequence2FeatureSequence());

        // Do the same thing for the "target" field: 
        //  convert a class label string to a Label object,
        //  which has an index in a Label alphabet.
        //pipeList.add(new Target2Label());

        // Now convert the sequence of features to a sparse vector,
        //  mapping feature IDs to counts.
        pipeList.add(new FeatureSequence2FeatureVector());
        
        // Print out the features and the label
        //pipeList.add(new PrintInputAndTarget());
        
        return (new SerialPipes(pipeList));
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		ExaminePages ep = new ExaminePages();
		try {
			ep.printBasicAnalysis();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
