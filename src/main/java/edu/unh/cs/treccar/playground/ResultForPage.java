package edu.unh.cs.treccar.playground;

import java.util.ArrayList;
import java.util.HashMap;

public class ResultForPage {
	private HashMap<String, ArrayList<String>> queryParaAssignment;
	private ArrayList<ArrayList<String>> paraClusters;
	private HashMap<ArrayList<String>, Double> queryParaRank;
	public HashMap<ArrayList<String>, Double> getQueryParaRank() {
		return queryParaRank;
	}
	public void setQueryParaRank(HashMap<ArrayList<String>, Double> queryParaRank) {
		this.queryParaRank = queryParaRank;
	}
	public HashMap<String, ArrayList<String>> getQueryParaAssignment() {
		return queryParaAssignment;
	}
	public void setQueryParaAssignment(
			HashMap<String, ArrayList<String>> queryParaAssignment) {
		this.queryParaAssignment = queryParaAssignment;
	}
	public ArrayList<ArrayList<String>> getParaClusters() {
		return paraClusters;
	}
	public void setParaClusters(ArrayList<ArrayList<String>> paraClusters) {
		this.paraClusters = paraClusters;
	}
}
