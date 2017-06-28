package edu.unh.cs.treccar.playground;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import edu.unh.cs.treccar.Data;
import edu.unh.cs.treccar.Data.Section;
import edu.unh.cs.treccar.read_data.DeserializeData;

public class TestOutlineFile {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String outlineFile = args[0];
		FileInputStream fileInputStream;
		String parent;
		try {
			fileInputStream = new FileInputStream(new File(outlineFile));
			for(Data.Page page: DeserializeData.iterableAnnotations(fileInputStream)) {
				System.out.println("Page id: "+page.getPageId());
				parent = page.getPageId();
				System.out.println("page skeleton: "+page.getSkeleton());
				System.out.println("Flat sectionpaths of page: "+page.flatSectionPaths());
				System.out.println("flatSectionPathsParagraphs of page: "+page.flatSectionPathsParagraphs());
				ArrayList<Section> sections = page.getChildSections();
				for(Data.Section sec : sections){
					TestOutlineFile.getSection(sec, parent);
				}
				//System.out.println(Data.sectionPathId(page.getPageId(), sections));
	        }
			fileInputStream.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        
	}
	
	private static void getSection(Data.Section section, String parent){
		System.out.println(parent+"\\"+section.getHeadingId());
		//System.out.println("Inside section...........");
		//System.out.println("heading "+section.getHeading());
		//System.out.println("heading ID: "+section.getHeadingId());
		
		//System.out.println("to string: "+section.toString());
		parent = parent+"\\"+section.getHeadingId();
		if(section.getChildSections() != null){
			for(Data.Section child : section.getChildSections()){
				getSection(child, parent);
			}
		}
		
	}

}
