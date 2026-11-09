# Assembly Anchor

The concept is to use a known gene annotation and exploit synteny to position contigs of a new assembly into pseudo chromosomes.
All known genes are aligned to contigs using blastn. Only best (first) blast hit is considered.
Contigs get orientation and position from the majority of genes that hit the contig and are consistent to each other with regards of position to each other in their original assembly and their position on the contig.
All contigs are then re-orientated and contigs with the same chromosome are concatenated adding Ns as spaces in between.


## Prerequisites
* blast+ [https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html]
* Java Runtime Environment 1.8 or higher [https://www.java.com/en/download/]


# Usage

## Step1 create blastable database of your assembly

```
makeblastdb -in assembly.fasta -db nucl
```

## Step 2 prepare reference gene set
Prepare your input gene set by adding the position of the gene to the description of each sequence.
Your fasta file of genes must have the format "`geneID chromosome:position`". The position may be the start position of the gene, taken from a gff file. 

Example
```
>gene1 chr1:50325
ATGCATGCATGCTGA
```


## Step3 run blast

Downstream pipeline requires XML output from blast (`-outfmt 5`)
Since only first blast hit is processed downstream, blast will run much faster by restricting number of reported hits (`-max_target_seqs 1`). This will also produce much smaller output blast files. 

```
blastn -query genes.fasta -db assembly.fasta -out blast.out.xml -outfmt 5 -max_target_seqs 1 -max_hsps 3
```

## Step 4 run AssemblyAnchor.jar

AssemblyAnchor.jar has three mandatory arguments and one optional: The blast output from step3, the assembly file and the path to the output file. An optional output can be generated for qc purpose.

-i blast.out.xml
-f assembly.fasta
-o output.fasta
-m gene_positions.txt

Example

```
java -jar AssemblyAnchor.jar -i blast.out.xml -f assembly.fasta -o assembly.anchored.fasta -m markers.txt
```
