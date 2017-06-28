package edu.unh.cs.treccar.playground.split;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.HashSet;

import co.nstant.in.cbor.CborDecoder;
import co.nstant.in.cbor.CborEncoder;
import co.nstant.in.cbor.CborException;
import co.nstant.in.cbor.model.DataItem;
import edu.unh.cs.treccar.Data;
import edu.unh.cs.treccar.read_data.DeserializeData;

public class SplitParagraph {
	String folder, gtFile, master;
	final String pFile = "all.test200.cbor.paragraphs";
	HashSet<String> paraIdSet = new HashSet<String>();
	public SplitParagraph(String folder, String gtFile, String paraFile) throws IOException{
		this.folder = folder;
		this.gtFile = gtFile;
		this.master = paraFile;
		this.paraIdSet = SplitDataset.getParaidsFromGT(folder+"/"+gtFile);
	}
	public void split() throws IOException, CborException{
		System.out.println();
		System.out.println("Para ids:");
		for(String id:this.paraIdSet){
			System.out.println(id);
		}
		File trainDir = new File(this.folder+"/train");
		File testDir = new File(this.folder+"/test");
		trainDir.mkdir();
		testDir.mkdir();
		File train = new File(trainDir.getAbsolutePath()+"/"+this.pFile);
		File test = new File(testDir.getAbsolutePath()+"/"+this.pFile);
		train.createNewFile();
		test.createNewFile();
		File input = new File(this.master);
		FileInputStream fis = new FileInputStream(input);
		FileOutputStream fosTr = new FileOutputStream(train);
		FileOutputStream fosTt = new FileOutputStream(test);
		CborDecoder decoder = new CborDecoder(fis);
		CborEncoder encoderTrain = new CborEncoder(fosTr);
		CborEncoder encoderTest = new CborEncoder(fosTt);
		DataItem curr = decoder.decodeNext();
		while(curr!=null){
			//check curr here
			if(checkItem(curr))
				encoderTrain.encode(curr);
			else
				encoderTest.encode(curr);
			curr = decoder.decodeNext();
		}
	}
	private boolean checkItem(DataItem item){
		boolean doEncode = false;
		Data.Paragraph itemPara = DeserializeData.paragraphFromCbor(item);
		doEncode = paraIdSet.contains(itemPara.getParaId());
		return doEncode;
	}
}
