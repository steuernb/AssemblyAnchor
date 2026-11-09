package hic_validation;


import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;

import support.BlastHSP;
import support.BlastHit;
import support.BlastIteration;

public class HiCMarker implements Comparable<HiCMarker> {

	String markerName;
	
	String geneticMapChromosome;
	double geneticMapCentiMorgan;
	
	String primaryContig;
	int primaryPosition;
	
	Vector<BlastHSP> hsps;
	
	
	
	
	
	
	public HiCMarker(BlastIteration it, int numIdenticalThreshold) {
		
		this.markerName = it.getQueryID();
		this.geneticMapChromosome = it.getQueryDescription().split(":")[0];
		this.geneticMapCentiMorgan = Double.parseDouble(it.getQueryDescription().split(":")[1].trim());
		
		BlastHSP primaryHSP = it.getFirstHSP();
		this.primaryContig = primaryHSP.getHitId();
		this.primaryPosition = primaryHSP.getHitStart();
		this.hsps = new Vector<BlastHSP>();
		for(Iterator<BlastHit> iterator = it.getHits().iterator(); iterator.hasNext();) {
			BlastHit hit = iterator.next();
			
			for(Iterator<BlastHSP> iterator2 = hit.getHSPs().iterator(); iterator2.hasNext();) {
				BlastHSP hsp = iterator2.next();
				if( hsp.getNumIdentical()>= numIdenticalThreshold) {
					hsps.add(hsp);
				}else {
					break; // too dangerous? hsps are ordered by score not numIdentical
				}
				
			}
			
			
		}
		
		
	}
	
	
	public String toString() {
		return markerName + "\t" + geneticMapChromosome + "\t" + geneticMapCentiMorgan + "\t" + primaryContig +"\t" + primaryPosition;
	}
	
	
	
	/**
	 * 
	 * Give a value for consistency with another contig. 
	 * 	0 means they match
	 * 	1 means same LG but different contigs/scaffolds
	 * 	2 means same contig/scaffold but different LG
	 *  3 means inconsistent order
	 * 
	 * 
	 * 
	 * @param marker2
	 * @param chromosomeDirection
	 * @return
	 */
	public int checkConsistencyWith(HiCMarker marker2 , HashMap<String, Boolean> chromosomeDirection) {
		if( this.geneticMapChromosome.equalsIgnoreCase(marker2.geneticMapChromosome)) {
			
			if(this.primaryContig.equalsIgnoreCase(marker2.primaryContig)) {
			
				if(chromosomeDirection.get(this.primaryContig).booleanValue()) {
					
					if( this.primaryPosition > marker2.primaryPosition && this.geneticMapCentiMorgan< marker2.geneticMapCentiMorgan) {
						return 3;
					}
					if(this.primaryPosition < marker2.primaryPosition && this.geneticMapCentiMorgan > marker2.geneticMapCentiMorgan) {
						return 3;
					}
					
				}else {
					
					
					if( this.primaryPosition > marker2.primaryPosition && this.geneticMapCentiMorgan> marker2.geneticMapCentiMorgan) {
						return 3;
					}
					if(this.primaryPosition < marker2.primaryPosition && this.geneticMapCentiMorgan < marker2.geneticMapCentiMorgan) {
						return 3;
					}
					
				}
				
				
				
			}else{
				//not sure if this makes sense. It will always flag the unanchored contigs
				return 1;
			}
			
			
			
			
			
		}else {
			// if they are on different genetic map chromosomes but on same contig it is considered inconsistent
			if( this.primaryContig.equalsIgnoreCase(marker2.primaryContig)) {
				return 2;
			}
		}
		
		
		
		return 0;
	}
	
	
	public boolean isConsistent(HiCMarker marker2 , HashMap<String, Boolean> chromosomeDirection) {
		
		if( this.geneticMapChromosome.equalsIgnoreCase(marker2.geneticMapChromosome)) {
			
			if(this.primaryContig.equalsIgnoreCase(marker2.primaryContig)) {
			
				if(chromosomeDirection.get(this.primaryContig).booleanValue()) {
					
					if( this.primaryPosition > marker2.primaryPosition && this.geneticMapCentiMorgan< marker2.geneticMapCentiMorgan) {
						return false;
					}
					if(this.primaryPosition < marker2.primaryPosition && this.geneticMapCentiMorgan > marker2.geneticMapCentiMorgan) {
						return false;
					}
					
				}else {
					
					
					if( this.primaryPosition > marker2.primaryPosition && this.geneticMapCentiMorgan> marker2.geneticMapCentiMorgan) {
						return false;
					}
					if(this.primaryPosition < marker2.primaryPosition && this.geneticMapCentiMorgan < marker2.geneticMapCentiMorgan) {
						return false;
					}
					
				}
				
				
				
			}else{
				//not sure if this makes sense. It will always flag the unanchored contigs
				return false;
			}
			
			
			
			
			
		}else {
			// if they are on different genetic map chromosomes but on same contig it is considered inconsistent
			if( this.primaryContig.equalsIgnoreCase(marker2.primaryContig)) {
				return false;
			}
		}
		
		
		
		return true;
	}
	
	
	
	
	public int compareTo(HiCMarker marker2) { // sort markers according to their position on contig. This comes from best blast hit.
		if( !this.primaryContig.equalsIgnoreCase(marker2.primaryContig)) {
			return this.primaryContig.compareTo(marker2.primaryContig);
		}else {
			
			if( this.primaryPosition <  marker2.primaryPosition) {
				return -1;
			}else if(this.primaryPosition > marker2.primaryPosition) {
				return 1;
			}else {
				return 0;
			}
		}
	}
	
	
	public String getPrimaryContig() {
		return this.primaryContig;
	}
	public int getPrimaryPosition() {
		return this.primaryPosition;
	}


	public String getMarkerName() {
		return markerName;
	}



	public String getGeneticMapChromosome() {
		return geneticMapChromosome;
	}


	public double getGeneticMapCentiMorgan() {
		return geneticMapCentiMorgan;
	}



	
	
}
