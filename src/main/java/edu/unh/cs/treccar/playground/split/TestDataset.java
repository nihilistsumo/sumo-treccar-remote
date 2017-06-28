package edu.unh.cs.treccar.playground.split;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import co.nstant.in.cbor.CborDecoder;
import co.nstant.in.cbor.CborException;
import co.nstant.in.cbor.model.DataItem;
import edu.unh.cs.treccar.read_data.DeserializeData;

public class TestDataset {
//Arguments: ground truth, outline, paragraph
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String gtFile = args[0];
		String outlineFile = args[1];
		String paraFile = args[2];
		try {
			FileInputStream pinput = new FileInputStream(new File(paraFile));
			FileInputStream oinput = new FileInputStream(new File(outlineFile));
			FileInputStream gtinput = new FileInputStream(new File(gtFile));
			CborDecoder paraDecoder = new CborDecoder(pinput);
			CborDecoder outlineDecoder = new CborDecoder(oinput);
			ArrayList<String> paraList = new ArrayList<String>();
			ArrayList<String> paraListGT = new ArrayList<String>();
			System.out.println("Starting...");
			for(DataItem item : paraDecoder.decode()){
				paraList.add(DeserializeData.paragraphFromCbor(item).getParaId());
			}
			BufferedReader br = new BufferedReader(new FileReader(gtFile));
			String line;
			while((line = br.readLine())!=null){
				paraListGT.add(line.split(" ")[2]);
			}
			for(String p:paraListGT){
				System.out.println(p);
			}
		} catch (CborException | IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
