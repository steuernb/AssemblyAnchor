package support;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.zip.GZIPInputStream;

public class CleanFasta {

	public static void main(String[] args) {
		try {
			
			File inputFile = new File("/Volumes/group-scratch/Matthew-Hartley/steuernb/RachelWells/anchor_napus_assemblies/BnapusAPG.anchoredZs11.fasta");
			File outputFile = new File("/Volumes/group-scratch/Matthew-Hartley/steuernb/RachelWells/anchor_napus_assemblies/BnapusAPG.anchoredZs11_clean.fasta");
			format(inputFile, outputFile, 100);
			
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	
	
	
	
	
	public static void format(File inputFile, File outputFile, int lineLength)throws IOException{
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		
		BufferedReader in;
		
		FileInputStream fis = new FileInputStream(inputFile);
		byte[] bytes = new byte[2];
		fis.read(bytes);
		int head = ((int) bytes[0] & 0xff) | ((bytes[1] << 8) & 0xff00);
		boolean gzip = GZIPInputStream.GZIP_MAGIC == head;
		fis.close();
		if(gzip){
			in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inputFile))));
		}else{
			in = new BufferedReader((new FileReader(inputFile)));
		}		
		
		long counta = 0;
		long numSeqs = 0;
		
		Hashtable<Long, Long> lengths = new Hashtable<Long, Long>();
		Hashtable<String, Long> sizes = new Hashtable<String, Long>();
		
		long seqCounta = 0;
		
		int baseCharacter=in.read();
		
		
		StringBuffer sequence = new StringBuffer("");
		String id= "";
		
		while(baseCharacter >=0){
			char c = (char)baseCharacter;
			
			if( c == '>'){
				
				long num = 0;
				if( lengths.containsKey(seqCounta)){
					num = lengths.get(seqCounta);
				}
				num++;
				lengths.put(seqCounta, num);
				if(id.length()>0) {
					sizes.put(id, seqCounta);
				}
				numSeqs++;
				if(sequence.length()>0){
					out.write(sequence.toString());
					out.newLine();
					sequence = new StringBuffer("");
					
				}
				String header = in.readLine();
				out.write(">" +  header);
				id = header.split("\\s")[0];
				out.newLine();
			}else{
				if(baseCharacter != '\n'){
					counta++;
					sequence.append(c);
					if(sequence.length()==lineLength){
						out.write(sequence.toString());
						out.newLine();
						sequence = new StringBuffer("");
					}
				}
			}
			
			
			
			baseCharacter = in.read();
		}
		
		out.write(sequence.toString());
		out.newLine();
		
		
		in.close();
		out.close();
		System.out.println(counta + "\t" + numSeqs);
		
		for(Iterator<String > iterator = sizes.keySet().iterator(); iterator.hasNext();) {
			String chr = iterator.next();
			System.out.println(chr + "\t" + sizes.get(chr).longValue());
		}
		
		
	}
	
	
	
	
}
