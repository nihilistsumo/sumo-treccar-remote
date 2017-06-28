package edu.unh.cs.treccar.playground;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;
import java.util.Set;

import org.apache.commons.lang3.builder.EqualsBuilder;
import org.apache.commons.lang3.builder.HashCodeBuilder;

import co.nstant.in.cbor.CborException;
import edu.unh.cs.treccar.Data;
import edu.unh.cs.treccar.read_data.DeserializeData;

public class AssignParagraphs {
	static class SectionPathID{
		private String sectionPathID;

		public String getSectionPathID() {
			return sectionPathID;
		}
		public SectionPathID(String id){
			this.sectionPathID = id;
		}
		public int hashCode(){
			return new HashCodeBuilder(11, 19).append(sectionPathID).toHashCode();
		}
		public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            SectionPathID obj = (SectionPathID) o;
            return new EqualsBuilder().append(sectionPathID, obj.getSectionPathID()).isEquals();
		}
	}
	public HashMap<AssignParagraphs.SectionPathID, ArrayList<Data.Paragraph>> assignCandidates(String outlinePath, String paraPath){
		HashMap<AssignParagraphs.SectionPathID, ArrayList<Data.Paragraph>> candidate = new HashMap<AssignParagraphs.SectionPathID, ArrayList<Data.Paragraph>>();
		ArrayList<AssignParagraphs.SectionPathID> sectionList = new ArrayList<AssignParagraphs.SectionPathID>();
		Iterable<Data.Paragraph> paraList = new ArrayList<Data.Paragraph>();
		Iterable<Data.Page> pageList = new ArrayList<Data.Page>();
		try {
			FileInputStream oulineStream = new FileInputStream(new File(outlinePath));
			FileInputStream paraStream = new FileInputStream(new File(paraPath));
			pageList = DeserializeData.iterableAnnotations(oulineStream);
			String parent;
			for(Data.Page page : pageList){
				parent = page.getPageId();
				ArrayList<Data.Section> sectionsOfPage = page.getChildSections();
				for(Data.Section section : sectionsOfPage){
					sectionList.add((new AssignParagraphs.SectionPathID(parent)));
					getChildren(section, sectionList, parent);
				}
			}
			paraList = DeserializeData.iterableParagraphs(paraStream);
			
			//*** This part will be replaced by clustering/topic modeling etc. ***//
			//*** sectonList, paraList ---> candidate                          ***//
			/*
			Random rand = new Random();
			AssignParagraphs.SectionPathID currentSection;
			for(Data.Paragraph para : paraList){
				currentSection = sectionList.get(rand.nextInt(sectionList.size()));
				if(candidate.containsKey(currentSection)){
					candidate.get(currentSection).add(para);
				} else{
					ArrayList<Data.Paragraph> currentParaList = new ArrayList<Data.Paragraph>();
					currentParaList.add(para);
					candidate.put(currentSection, currentParaList);
				}	
			}
			*/
			//*** ############################################################ ***//
			/*
			FileInputStream pStream = new FileInputStream(new File(paraPath));
			Iterable<Data.Paragraph> pList = new ArrayList<Data.Paragraph>();
			pList = DeserializeData.iterableParagraphs(pStream);
			*/
			TopicModel tm = new TopicModel(sectionList, paraList);
			tm.modelTopics();
			candidate = tm.getCandidate();
			
		} catch (FileNotFoundException | CborException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return candidate;
	}
	private static void getChildren(Data.Section section, ArrayList<AssignParagraphs.SectionPathID> sectionList, String parent){
		sectionList.add(new AssignParagraphs.SectionPathID(parent+"/"+section.getHeadingId()));
		parent = parent+"/"+section.getHeadingId();
		if(section.getChildSections() != null){
			for(Data.Section child : section.getChildSections()){
				getChildren(child, sectionList, parent);
			}
		}
		
	}

}
