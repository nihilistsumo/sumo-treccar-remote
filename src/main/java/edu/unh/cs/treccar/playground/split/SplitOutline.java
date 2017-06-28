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

public class SplitOutline {
	String folder, gtFile, master;
	final String oFile = "all.test200.cbor.outlines";
	HashSet<String> pageIdSet = new HashSet<String>();
	public SplitOutline(String folder, String gtFile, String outlineFile) throws IOException{
		this.folder = folder;
		this.gtFile = gtFile;
		this.master = outlineFile;
		this.pageIdSet = SplitDataset.getPageidsFromGT(folder+"/"+gtFile);
	}
	public void split() throws IOException, CborException{
		System.out.println();
		System.out.println("Page ids:");
		for(String id:this.pageIdSet){
			System.out.println(id);
		}
		File trainDir = new File(this.folder+"/train");
		File testDir = new File(this.folder+"/test");
		trainDir.mkdir();
		testDir.mkdir();
		File train = new File(trainDir.getAbsolutePath()+"/"+this.oFile);
		File test = new File(testDir.getAbsolutePath()+"/"+this.oFile);
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
		Data.Page itemPage = DeserializeData.pageFromCbor(item);
		//check whether pageid of itemPage matches any from gtfile
		doEncode = this.pageIdSet.contains(itemPage.getPageId());
		return doEncode;
	}
}
