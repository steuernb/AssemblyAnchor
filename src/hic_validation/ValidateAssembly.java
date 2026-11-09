package hic_validation;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;

import support.BlastIPKReader;
//import support.BlastIPKReader;
import support.BlastIteration;
import support.BlastReader;
import support.BlastXMLReader;

public class ValidateAssembly {
	
	HashMap<String,HiCMarker> markers;
	HashMap<String, Boolean> chromosomeDirection;
	
	
	
	
	public ValidateAssembly(File blastFile)throws IOException {
		
		this.markers = new HashMap<String,HiCMarker>();
		
		int counta = 0;
		
		
		
		BlastReader reader;
		
		
		BufferedReader in = new BufferedReader(new FileReader(blastFile));
		String inputline = in.readLine();
		if( inputline.startsWith("<?xml")) {
			reader = new BlastXMLReader(blastFile);
		}else {
			reader = new BlastIPKReader(blastFile);
		}
		
		in.close();
		
		
	
		
		
		
		for (BlastIteration it = reader.readIteration(); it != null; it = reader.readIteration()) {
			
			HiCMarker marker = new HiCMarker(it, 50);
			markers.put(marker.markerName, marker);
			counta++;
			if(counta%100==0) {
				System.out.print(".");
			}
			if(counta%10000==0) {
				System.out.println("");
			}
		}
		reader.close();
		System.out.println("Blast output read. " + markers.size() + " markers");
		
		this.determineContigDirection();
	}

	
	
	public void writeInconsistencyReport(File outputFile, File fastaIndex)throws IOException{
		
		HashMap<String, Long> offsets = this.getOffsets(fastaIndex);
		
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		
		HashSet<String> incompatibleMarkerNames = new HashSet<String>();
		
		Vector<HiCMarker> v = new Vector<HiCMarker>();
		for(Iterator<String> iterator = markers.keySet().iterator(); iterator.hasNext();) {
			String markerName = iterator.next();
			v.add(markers.get(markerName));
		}
		
		HashMap<String, HashMap<String, Integer>> incompatible = new HashMap<String, HashMap<String, Integer>>();
		
		for(int i = 0; i<v.size(); i++) {
			HiCMarker marker1 = v.get(i);
			for(int j = i+1; j<v.size(); j++ ) {
				HiCMarker marker2 = v.get(j);
				
				int scoring = marker1.checkConsistencyWith(marker2, chromosomeDirection);
				
				if( scoring != 0) {
					if(!incompatible.containsKey(marker1.markerName)) {
						incompatible.put(marker1.markerName, new HashMap<String, Integer> ());
					}
					if(!incompatible.containsKey(marker2.markerName)) {
						incompatible.put(marker2.markerName, new HashMap<String, Integer>  ());
					}
					incompatible.get(marker1.markerName).put(marker2.markerName, scoring);
					incompatible.get(marker2.markerName).put(marker1.markerName, scoring);
					
				}
				
			}
		}
		
		while( incompatible.size() > 0) {
			
			//step 1: find marker with max incompatibilities.
			String max = "";
			int num = 0;
			for(Iterator<String> iterator = incompatible.keySet().iterator(); iterator.hasNext(); ) {
				String markerName = iterator.next();
				int numIncompatible = incompatible.get(markerName).size();
				if( numIncompatible > num ) {
					max = markerName;
					num = numIncompatible;
				}
			}
			
			if(num > 0) {
				//step 2: remove that marker and increment total number of inconsistencies
				
				
				int num1=0; //unanchored
				int num2=0; //wrong LG
				int num3=0; //wrong order
				
				for( Iterator<String> iterator = incompatible.get(max).keySet().iterator(); iterator.hasNext();) {
					String otherMarker = iterator.next();
					int scoring = incompatible.get(max).get(otherMarker);
					if(scoring ==1) {num1++;}
					if(scoring ==2) {num2++;}
					if(scoring ==3) {num3++;}
				}
				String color = "green";
				if(num1 >=num2 && num1>=num3) { color = "red";}
				if(num2 > num1 && num2 >=num3) {color = "blue";}
				
				
				
				HiCMarker marker = markers.get(max);
				out.write(max + "\t" + marker.geneticMapChromosome + "\t" + marker.geneticMapCentiMorgan 
						  + "\t" + marker.primaryContig + "\t" + marker.primaryPosition + "\t"  
						  + (offsets.get(marker.primaryContig).longValue() + marker.primaryPosition)  +"\t" + color );
				out.newLine();
				
				incompatibleMarkerNames.add(max);
				
				
				incompatible.remove(max);
				HashSet<String> remove = new HashSet<String>();
				for(Iterator<String> iterator = incompatible.keySet().iterator(); iterator.hasNext(); ) {
					String markerName = iterator.next();
					incompatible.get(markerName).remove(max);
					if(incompatible.get(markerName).size() == 0) {
						remove.add(markerName);
					}
				}
				for(Iterator<String > iterator = remove.iterator(); iterator.hasNext();) {
					incompatible.remove(iterator.next());
				}
				
				
				
				
				
				
			}
			
			
		}
		
		
		for(Iterator<String> iterator = markers.keySet().iterator();iterator.hasNext();) {
			String markerName = iterator.next();
			if(!incompatibleMarkerNames.contains(markerName)) {
				HiCMarker marker = markers.get(markerName);
				out.write(markerName + "\t" + marker.geneticMapChromosome + "\t" + marker.geneticMapCentiMorgan 
						  + "\t" + marker.primaryContig + "\t" + marker.primaryPosition + "\t"  
						  + (offsets.get(marker.primaryContig).longValue() + marker.primaryPosition)  +"\t" + "black" );
				out.newLine();
			}
		}
		
		
		
		
		
		
		out.close();
	}
	
	
	
	private HashMap<String, Long> getOffsets(File faidx)throws IOException{
		
		
		
		HashMap<String, Long> h = new HashMap<String, Long>();
		
		
		BufferedReader in = new BufferedReader(new FileReader(faidx));

		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			String[] split = inputline.split("\t");
			h.put(split[0], Long.parseLong(split[2]));
		}

		in.close();
		
		
		return h;
		
	}
	
	
	public int getNumInconsistent() {
		
		
		
		
		int totalNumberOfInconsistentMarkers = 0;
		
		
		
		
		Vector<HiCMarker> v = new Vector<HiCMarker>();
		for(Iterator<String> iterator = markers.keySet().iterator(); iterator.hasNext();) {
			String markerName = iterator.next();
			
			v.add(markers.get(markerName));
		}
		
		HashMap<String, HashSet<String>> incompatible = new HashMap<String, HashSet<String>>();
		
		for(int i = 0; i<v.size(); i++) {
			HiCMarker marker1 = v.get(i);
			for(int j = i+1; j<v.size(); j++ ) {
				HiCMarker marker2 = v.get(j);
				
				if( !marker1.isConsistent(marker2, chromosomeDirection)) {
					if(!incompatible.containsKey(marker1.markerName)) {
						incompatible.put(marker1.markerName, new HashSet<String> ());
					}
					if(!incompatible.containsKey(marker2.markerName)) {
						incompatible.put(marker2.markerName, new HashSet<String> ());
					}
					incompatible.get(marker1.markerName).add(marker2.markerName);
					incompatible.get(marker2.markerName).add(marker1.markerName);
				}
				
			}
		}
		
		System.out.println("initial incompatibilities: " + incompatible.size());
		
		while( incompatible.size() > 0) {
			
			//step 1: find marker with max incompatibilities.
			String max = "";
			int num = 0;
			for(Iterator<String> iterator = incompatible.keySet().iterator(); iterator.hasNext(); ) {
				String markerName = iterator.next();
				int numIncompatible = incompatible.get(markerName).size();
				if( numIncompatible > num ) {
					max = markerName;
					num = numIncompatible;
				}
			}
			
			if(num > 0) {
				//step 2: remove that marker and increment total number of inconsistencies
				totalNumberOfInconsistentMarkers++;
				incompatible.remove(max);
				HashSet<String> remove = new HashSet<String>();
				for(Iterator<String> iterator = incompatible.keySet().iterator(); iterator.hasNext(); ) {
					String markerName = iterator.next();
					incompatible.get(markerName).remove(max);
					if(incompatible.get(markerName).size() == 0) {
						remove.add(markerName);
					}
				}
				for(Iterator<String > iterator = remove.iterator(); iterator.hasNext();) {
					incompatible.remove(iterator.next());
				}
				
			}
			
			
		}
		
		return totalNumberOfInconsistentMarkers;
		
	}
	
	public void determineContigDirection() {
		System.out.println("Determining contig directions");
		HashMap<String, int[]> mychromosomeDirection = new HashMap<String,int[]>();
		
		Vector<HiCMarker> v = new Vector<HiCMarker>();
		for(Iterator<String> iterator = markers.keySet().iterator(); iterator.hasNext();) {
			v.add(markers.get(iterator.next()));
		}
		
		for( int i = 0; i< v.size(); i++) {
			HiCMarker marker1 = v.get(i);
			
			if(!mychromosomeDirection.containsKey(marker1.primaryContig)) {
				mychromosomeDirection.put(marker1.primaryContig, new int[2]);
			}
			
			for( int j = i+1; j< v.size();j++) {
				HiCMarker marker2 = v.get(j);
				
				if(marker1.primaryContig.equalsIgnoreCase(marker2.primaryContig)) {
					if(marker1.geneticMapChromosome.equalsIgnoreCase(marker2.geneticMapChromosome)) {
						
						if( marker1.primaryPosition > marker2.primaryPosition && marker1.geneticMapCentiMorgan< marker2.geneticMapCentiMorgan ) {
							mychromosomeDirection.get(marker1.primaryContig)[1]++;
						}
						if( marker1.primaryPosition > marker2.primaryPosition && marker1.geneticMapCentiMorgan> marker2.geneticMapCentiMorgan ) {
							mychromosomeDirection.get(marker1.primaryContig)[0]++;
						}
						if( marker1.primaryPosition < marker2.primaryPosition && marker1.geneticMapCentiMorgan< marker2.geneticMapCentiMorgan ) {
							mychromosomeDirection.get(marker1.primaryContig)[0]++;
						}
						if( marker1.primaryPosition < marker2.primaryPosition && marker1.geneticMapCentiMorgan> marker2.geneticMapCentiMorgan ) {
							mychromosomeDirection.get(marker1.primaryContig)[1]++;
						}
					
						
					}
					
				}
				
				
			}
			
		}
		this.chromosomeDirection = new HashMap<String, Boolean>();
		
		for(Iterator<String> iterator = mychromosomeDirection.keySet().iterator(); iterator.hasNext();) {
			String chr = iterator.next();
			int[] a = mychromosomeDirection.get(chr);
			boolean forward = true;
			if(a[1] > a[0]) {
				forward = false;
			}
			
			
			
			if(forward) {
				//System.out.println(chr + " direction is forward");
			}else {
				//System.out.println(chr + " direction is reverse");
			}
			
			
			chromosomeDirection.put(chr, forward);
		}
		
	}
	
	public void writeMarkerPairs(File outputFile)throws IOException{
		
		Vector<HiCMarker> v = new Vector<HiCMarker>();
		for(Iterator<String> iterator = markers.keySet().iterator(); iterator.hasNext();) {
			String markerName = iterator.next();
			v.add(markers.get(markerName));
		}
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		
		for(int i = 0; i < v.size(); i++) {
			HiCMarker marker1 = v.get(i);
			
			for( int j = 0; j< v.size(); j++) {
				HiCMarker marker2 = v.get(j);
				
				int score = marker1.checkConsistencyWith(marker2, chromosomeDirection);
				if(score !=0) {
					out.write(i + "\t" + j + "\t" + score);
					out.newLine();
				}
				
			}
			
			
			
		}
		
		out.close();
		
	}
	
	
	public HashMap<String, HiCMarker> getMarkers(){
		return this.markers;
	}
	
	public boolean getChromosomeDirection(String chromosome) {
		return this.chromosomeDirection.get(chromosome);
	}
}
