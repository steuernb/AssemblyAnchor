package AnchorAssembly;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;

import hic_validation.HiCMarker;
import hic_validation.ValidateAssembly;

import support.BioSequence;
import support.FastaReader;
import support.commandLineInterface.CLI;

public class AnchorAssembly extends ValidateAssembly {

	HashMap<String, Vector<String>> contigOrder; // chromosome<contigID>
	HashMap<String, HashSet<HiCMarker>> contigs;
	HashMap<String, String> contig2chr;
	
	
	public static void main(String[] args) {
		
		mainCLI(args);
		
		
	}
	

	
	public static void mainCLI(String[] args) {
		
		CLI cli = new CLI();
		cli.parseOptions(args);
		
		
		
	
		
		
		
		if(!cli.hasOption("i")) {
			System.out.println("please specify input blast file with option -i");
		}
		
		if( !cli.hasOption("f")) {
			System.out.println("please specify input genome using option -f");
		}
		
		
		if(! cli.hasOption("o")) {
			System.out.println("please specify output file using option -o");
		}
		File blastFile = new File(cli.getArg("i"));
		File fastaFile = new File(cli.getArg("f"));
		File outputFile = new File(cli.getArg("o"));
		
		try {
			
		
			
			AnchorAssembly aa = new AnchorAssembly(blastFile);
			
			aa.writeAssembly(outputFile, fastaFile);
			
			if( cli.hasOption("m")) {
				File outputMarkerFile = new File(cli.getArg("m"));
				aa.writeMarkers(outputMarkerFile, fastaFile);
			}
			
		
		
		} catch (Exception e) {
			
		}
		
		
	}
	
	
	
	
	public AnchorAssembly(File blastFile)throws IOException{
		super(blastFile);
		
		this.collectContigInformation();
		
		
	}
	
	
	/**
	 * This already assumes both contigs are known to be on the same chromosome.
	 * 
	 * 
	 * @param contig1
	 * @param contig2
	 * @return 
	 *			-1 if the first contig is first
	 *			 1 if the second contig is first
	 * 			 
	 */
	private int compareTwoContigs(String contig1, String contig2) {
		
		int before = 0;
		int after = 0;
		
		String chr = contig2chr.get(contig1);
		if(!chr.equalsIgnoreCase(contig2chr.get(contig2))){
			return 0;
		}
		
		
		for(Iterator<HiCMarker> iterator1 = contigs.get(contig1).iterator(); iterator1.hasNext();) {
			HiCMarker marker1 = iterator1.next();
			
			if( !marker1.getGeneticMapChromosome().equalsIgnoreCase(chr)) {
				continue;
			}
			
			for(Iterator<HiCMarker> iterator2 = contigs.get(contig2).iterator(); iterator2.hasNext();) {
				HiCMarker marker2 = iterator2.next();
				
				if(!marker2.getGeneticMapChromosome().equalsIgnoreCase(chr)) {
					continue;
				}
				if( marker1.getGeneticMapCentiMorgan() < marker2.getGeneticMapCentiMorgan()) {
					before++;
				}else if(marker1.getGeneticMapCentiMorgan() > marker2.getGeneticMapCentiMorgan()) {
					after++;
				}
				
			}
			
		}
		
		if(before>after) {
			return -1;
		}
		if(after>before) {
			return 1;
		}
		
		return 0;
	}
	
	private String getChromosomeForContig(String contig) {
		
		HashMap<String, Integer> chromosomeCounts = new HashMap<String, Integer>();
		
		for(Iterator<HiCMarker> iterator = contigs.get(contig).iterator(); iterator.hasNext();) {
			HiCMarker marker = iterator.next();
			String chr = marker.getGeneticMapChromosome();
			int count = 0;
			if( chromosomeCounts.containsKey(chr)) {
				count = chromosomeCounts.get(chr);
			}
			count++;
			chromosomeCounts.put(chr,count);
		}
		String winner = "";
		int largest = 0;
		for( Iterator<String> iterator = chromosomeCounts.keySet().iterator(); iterator.hasNext();) {
			String chr = iterator.next();
			int count = chromosomeCounts.get(chr);
			if( count> largest) {
				winner = chr;
				largest = count;
			}
		}
		
		if(!contigOrder.containsKey(winner)) {
			contigOrder.put(winner, new Vector<String>() );
			
		}
		contigOrder.get(winner).add(contig);
		
		return winner;
	}
	
	
	public void collectContigInformation() {
		
		System.out.println();
		
		
		this.contigs = new HashMap<String, HashSet<HiCMarker>>();
		this.contigOrder = new HashMap<String, Vector<String>> ();
		this.contig2chr = new HashMap<String, String>();
		
		
		for(Iterator<String> iterator = this.getMarkers().keySet().iterator(); iterator.hasNext();) {
			String markerID = iterator.next();
			
			String contig = this.getMarkers().get(markerID).getPrimaryContig();
			if(!contigs.containsKey(contig)) {
				contigs.put(contig, new HashSet<HiCMarker>());
			}
			contigs.get(contig).add(this.getMarkers().get(markerID));
		}
		
		for(Iterator<String> iterator = contigs.keySet().iterator(); iterator.hasNext();) {
			String contig = iterator.next();
			String chr = getChromosomeForContig(contig); //this also adds the contig to contig order
			//System.out.println(contig + "\t" + chr);
			contig2chr.put(contig, chr);
		}
		
		for( Iterator<String> iterator = contigOrder.keySet().iterator(); iterator.hasNext();) {
			String chr = iterator.next();
			Collections.sort(this.contigOrder.get(chr), new Comparator<String>() {public int compare(String contig1, String contig2) {return compareTwoContigs(contig1,contig2);}} );
		}
		
	}
	
	public void writeAssembly(File outputFile, File inputFasta)throws IOException {
		System.out.println("writing assembly");
		
		HashMap<String, BioSequence> inputAssembly = new HashMap<String, BioSequence>();
		FastaReader fastaReader = new FastaReader(inputFasta);
		for (BioSequence seq = fastaReader.readEntry(); seq != null; seq = fastaReader.readEntry()) {
			inputAssembly.put(seq.getIdentifier(), seq);
		}
		fastaReader.close();
		
		HashSet<String> anchored= new HashSet<String>();
		
		Vector<String > chromosomes = new Vector<String>();
		chromosomes.addAll(contigOrder.keySet());
		Collections.sort(chromosomes);
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
		for(Iterator<String> iterator = chromosomes.iterator();iterator.hasNext();) {
			String chromosome = iterator.next();
			
			StringBuffer buffer = new StringBuffer();
			
			
			Vector<String> contigs = contigOrder.get(chromosome);
			for(int i = 0; i< contigs.size(); i++) {
				String contig = contigs.get(i);
				anchored.add(contig);
				
				BioSequence seq = inputAssembly.get(contig);
				if( super.getChromosomeDirection(contig)) {
					buffer.append(seq.getSequence());
				}else {
					buffer.append(seq.getReverseComplementarySequence());
				}
				
				inputAssembly.remove(contig);
				if( i+1 < contigs.size()) {
					for( int j = 0; j< 50; j++) {
						buffer.append("N");
					}
				}
				
				
				
			}
			
			BioSequence newSequence= new BioSequence(chromosome, buffer.toString());
			newSequence.writeFasta(out, 100);
			
		}
		
		Vector<String> v = new Vector<String>();
		v.addAll(inputAssembly.keySet());
		Collections.sort(v);
		
		for(Iterator<String> iterator = v.iterator(); iterator.hasNext();) {
			String contig = iterator.next();
			inputAssembly.get(contig).writeFasta(out, 100);
		}

		
		
		
		
		out.close();
	}
	
	public void writeMarkers(File outputFile, File inputFasta)throws IOException{
		System.out.println("writing marker comparison");
		HashMap<String, BioSequence> inputAssembly = new HashMap<String, BioSequence>();
		FastaReader fastaReader = new FastaReader(inputFasta);
		for (BioSequence seq = fastaReader.readEntry(); seq != null; seq = fastaReader.readEntry()) {
			inputAssembly.put(seq.getIdentifier(), seq);
		}
		fastaReader.close();
		
		HashSet<String> anchored= new HashSet<String>();
		
		Vector<String > chromosomes = new Vector<String>();
		chromosomes.addAll(contigOrder.keySet());
		Collections.sort(chromosomes);
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
		for(Iterator<String> iterator = chromosomes.iterator();iterator.hasNext();) {
			String chromosome = iterator.next();
			
			int posOffset = 0;
			
			Vector<String> contigs = contigOrder.get(chromosome);
			for(int i = 0; i< contigs.size(); i++) {
				String contig = contigs.get(i);
				anchored.add(contig);
				
				boolean forward = super.getChromosomeDirection(contig);
				
				HashSet<HiCMarker> markers = this.contigs.get(contig);
				
				for(Iterator<HiCMarker> iterator3 = markers.iterator(); iterator3.hasNext();) {
					HiCMarker marker = iterator3.next();
					
					String name = marker.getMarkerName();
					String markerChromosome = marker.getGeneticMapChromosome();
					if(!markerChromosome.equalsIgnoreCase(chromosome)) {
						continue;
					}
					double markerPos = marker.getGeneticMapCentiMorgan();
					int positionOnContig = marker.getPrimaryPosition();
					
					int positionOnChromosome = posOffset + positionOnContig;
					if(!forward) {
						positionOnChromosome = posOffset + inputAssembly.get(contig).getLength()-positionOnContig;
					}
					
					
					out.write(name+"\t" + markerChromosome + "\t" + markerPos + "\t" + chromosome +"\t" + positionOnChromosome + "\t" + contig + "\t" + positionOnContig + "\n");
				}
				
				posOffset = posOffset + inputAssembly.get(contig).getLength() + 50;
			}
			
			
			
		}
		
		
		
		
		
		
		out.close();
	}
}
