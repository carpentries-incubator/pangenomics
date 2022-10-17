---
title: "Annotating Genomic Data"
teaching: 30
exercises: 15
questions:
- "How to download NCBI genomic data from the command line?"
- "How to annotate genome FASTA files?"
objectives:
- "Explore ncbi-genome-download as a tool for genomic data fetching from the NCBI."
- "Learn how to use the Prokka genome annotation utility."
keypoints:
- "The `ncbi-genome-download` package is a set of scripts to download genomes from the NCBI."
- "Prokka is a command line utility that provides rapid prokaryotic genome annotation."
---

## Prokka: Annotating Genomes

[Annotation](https://en.wikipedia.org/wiki/DNA_annotation) is the process of 
identifying the coordinates of genes and all the coding regions 
in a genome and determining what those genes are for. In order to do this, an unknown 
sequence is enriched with information relating genomic position, regulatory
sequences, repeats, gene name and protein products. This information
is stored in genomic databases to help future analysis processing new data.

[Prokka](https://github.com/tseemann/prokka) 
is a command-line software tool created in Perl to annotate bacterial, 
archaeal and viral genomes and reproduce standards-compliant output files.
It requires a preassembled genomic DNA sequences in FASTA format as input 
file, which is the only mandatory parameter to the software.
For annotation, Prokka relies on external features and databases to 
identify the genomic features within the contigs.

| Tool(reference) | Features predicted |
| --------- | ----------- |
|Prodigal (Hyatt 2010 )   | Coding Sequence (CDS) |
| RNAmmer ( Lagesen et al. , 2007 )  | Ribosomal RNA genes (rRNA) |
| Aragorn ( Laslett and Canback, 2004 )  | Transfer RNA genes |
| SignalP ( Petersen et al. , 2011 )  | Signal leader peptides|
| Infernal ( Kolbe and Eddy, 2011 )  | Non-coding RNA|

Proteins coding genes are annotated in two stages. Prodigal identifies 
the coordinates of candidate genes, but does not describe the putative 
gene product. Usually, in order to predict what a gene encodes 
for, it is compared with a large
database of known sequences, usually at a protein level, 
and transfer the annotation of the best significant match.
Prokka uses this method, but in a hierarchical manner. It starts 
with a small trustworthy database, it then moves to medium
sized but domain specific databases and finally to curated 
models of protein families.

> ## Notes
> [Environment variables](https://opensource.com/article/19/8/what-are-environment-variables) are special variables that contain information about your loggin session. These can be useful when you want to [manage](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#updating-an-environment) default or new settings that your system usually ignores. For instance, you can download a specific package with a downgraded version of perl if needed. 
{: .callout}

Next, we need to change to the directory where we have the assembly (FASTA) 
files of interest. As a simple initial example of execution, we can annotate 
a FASTA file and define names for our output directory and files like this:

~~~
$ mkdir -p ~/gm_workshop/results/annotated/
$ cd ~/gm_workshop/results/annotated/
$ conda deactivate
$ conda activate Prokka_Global
~~~
{: .language-bash}

Now you must be inside the environment  
~~~
(Prokka_Global) $
~~~
{: .language-bash}

You are ready to run your first annotation.  
~~~
$ prokka --prefix thermophilusLMD9_prokka --outdir thermophilusLMD9_prokka --kingdom Bacteria --genus Streptococcus --strain LMD9 --usegenus --addgenes ~/gm_workshop/data/thermophilusLMD9/GCF_000014485.1_ASM1448v1_genomic.fna
~~~
{: .language-bash}

In this example, we have told prokka to:

| code   | meaning |
| ------- | ---------- |
| --prefix | Filename output prefix [auto] (default '') |
| --outdir | Output folder [auto] (default '') |
| --kingdom | Annotation mode: Archaea|Bacteria|Mitochondria|Viruses (default 'Bacteria')|
| --genus | Genus name (default 'Genus') |
| --strain | Strain name (default 'strain') |
| --usegenus | Use genus-specific BLAST databases (needs --genus) (default OFF) |
| --addgens |Add 'gene' features for each 'CDS' feature (default OFF) |


Now prokka has generated a new folder. Lets get in and if you run the `tree` command 
in the current directory, you can preview the system of files created by Prokka:

~~~
exdir/
├── thermophilusLMD9_prokka.err
├── thermophilusLMD9_prokka.faa
├── thermophilusLMD9_prokka.ffn
├── thermophilusLMD9_prokka.fna
├── thermophilusLMD9_prokka.fsa
├── thermophilusLMD9_prokka.gbk
├── thermophilusLMD9_prokka.gff
├── thermophilusLMD9_prokka.log
├── thermophilusLMD9_prokka.sqn
├── thermophilusLMD9_prokka.tbl
├── thermophilusLMD9_prokka.tsv
└── thermophilusLMD9_prokka.txt
~~~
{: .output}

We encourage you to explore each output nevertheless, the following table describes the contents of each output file:

| Extension | Description |
| --------- | ----------- |
| .gff | This is the master annotation in GFF3 format, containing both sequences and annotations. It can be viewed directly in Artemis or IGV. |
| .gbk | This is a standard Genbank file derived from the master .gff. If the input to prokka was a multi-FASTA, then this will be a multi-Genbank, with one record for each sequence. |
| .fna | Nucleotide FASTA file of the input contig sequences. |
| .faa | Protein FASTA file of the translated CDS sequences. |
| .ffn | Nucleotide FASTA file of all the prediction transcripts (CDS, rRNA, tRNA, tmRNA, misc_RNA) |
| .sqn | An ASN1 format "Sequin" file for submission to Genbank. It needs to be edited to set the correct taxonomy, authors, related publication etc. |
| .fsa | Nucleotide FASTA file of the input contig sequences, used by "tbl2asn" to create the .sqn file. It is mostly the same as the .fna file, but with extra Sequin tags in the sequence description lines. |
| .tbl | Feature Table file, used by "tbl2asn" to create the .sqn file. |
| .err | Unacceptable annotations - the NCBI discrepancy report. |
| .log | Contains all the output that Prokka produced during its run. This is a record of what settings you used, even if the --quiet option was enabled. |
| .txt | Statistics relating to the annotated features found. |
| .tsv | Tab-separated file of all features: locus_tag,ftype,len_bp,gene,EC_number,COG,product |

You can also modify parameters as much as you need regarding the organism, gene and even locus tag you are looking for. 


> ## Exercise 1: Extracting tRNAs with Prokka
>
> Suppose you are now asked to annotate the FASTA file you downloaded in Exercise 1 and output the results to a subdirectory called `annotated` within the `thermophilusLMG18311_prokka` directory. Then, a research team asks you to provide them a TSV file titled `trnas.tsv` that only contains *S. thermophilus'*s tRNAs. This file must contain the same headers as the original TSV file, followed by the rows that correspond to tRNAs. Complete the following sequence of commands to perform this actions:
> 
> ~~~
> $ prokka --outdir thermophilusLMG18311_prokka --prefix thermophilusLMG18311_prokka ../../data/thermophilusLMG18311/__________ --kingdom Bacteria --genus Streptococcus --species thermophilus --usegenus --addgenes 
> ~~~
> {: .language-bash}
> 
> ~~~
> $ cd __________
> ~~~
> {: .language-bash}
> 
> ~~~
> $ __________ -n 1 thermophilusLMG18311_prokka.tsv > trnas.tsv  # Get column headers
> $ grep __________ thermophilusLMG18311_prokka.tsv >> trnas.tsv # Append all lines that contain tRNAs
> ~~~
> {: .language-bash}
>
> > ## Solution
> >
> > First, we perform the annotation with Prokka and save all files as `thermophilusLMG18311_prokka`.
> >
> > ~~~
> > $ prokka --prefix thermophilusLMG18311_prokka --outdir thermophilusLMG18311_prokka --kingdom Bacteria --genus Streptococcus --strain LMG18311 --usegenus --addgenes ../../data/thermophilusLMG18311/GCF_000011825.1_ASM1182v1_genomic.fna
> > ~~~
> > {: .language-bash}
> >
> > After switching to the `thermophilusLMG18311_prokka` directory, we shall now filter the data we need and save the outputs to a file named `trnas.tsv`. To do so, we use the `head` command with the `-n 1` argument to get the first line (the headers of the columns). We then append the lines that correspond to tRNAs, which is done with the code `$'\t'tRNA$'\t'`(this means that the program will search for lines that contain the word `tRNA` with tab spaces at the beginning and the end of the word).
> >
> > ~~~
> > $ cd thermophilusLMG18311_prokka
> > $ head -n 1 thermophilusLMG18311_prokka.tsv > trnas.tsv # Get column headers
> > $ grep $'\t'tRNA$'\t' thermophilusLMG18311_prokka.tsv >> trnas.tsv # Append all lines that contain tRNA
> > ~~~
> > {: .language-bash}
> >
> {: .solution}
{: .challenge}

> ## Discussion 1: Number of tRNAs
> 
> Inside the `annotated` directory, run the command `wc -l trnas.tsv` to get the number of lines in the file. Excluding the first line (which contains the header), observe that there are 67 tRNAs, whereas the [standard codon table](https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables#Translation_table_1) contains 61 coding codons (i.e. the ones that are not stop codons). This means that there are several genes that produce the same tRNA. How would you explain this fact?
> 
> > ## Solution
> > 
> > The existence of many genes producing the same product (such as the same tRNA) can happen due to duplications of genes during the evolutionary history of a taxon.
> {: .solution}
{: .discussion}




Lets first obtain the strain name of each fasta.  
~~~
$ cd ~/gm_workshop/data
$ ls */*fasta |while read line; do strain=$(echo $line|cut -d'_' -f3 |cut -d'.' -f1); echo $strain; done
~~~
{: .language-bash}

~~~
18RS21
515
A909
CJB111
COH1
H36B
~~~
{: .output}

You are ready to run all your annotations.  
~~~
$ ls */*fasta | while read line
> do 
> strainName=$(echo $line|cut -d'_' -f3 |cut -d'.' -f1)
> echo prokka $line --kingdom Bacteria --genus Streptococcus --species agalactie --strain $strainName --usegenus --addgenes --prefix Streptococcus_agalactie_${strainName}\.prokka --outdir ~/gm_workshop/results/annotated/Streptococcus_agalactie_${strainName}\_prokka
> done
~~~
{: .language-bash}

~~~
...
prokka COH1/Streptococcus_agalactiae_COH1.fasta --kingdom Bacteria --genus Streptococcus --species agalactie --strain COH1 --usegenus --addgenes --prefix Streptococcus_agalactie_COH1.prokka --outdir /home/alumno17/gm_workshop/results/annotated/Streptococcus_agalactie_COH1_prokka
prokka H36B/Streptococcus_agalactiae_H36B.fasta --kingdom Bacteria --genus Streptococcus --species agalactie --strain H36B --usegenus --addgenes --prefix Streptococcus_agalactie_H36B.prokka --outdir /home/alumno17/gm_workshop/results/annotated/Streptococcus_agalactie_H36B_prokka
~~~
{: .output}

~~~
$ cd ~/gm_workshop/data 
$ ls */*fasta | while read line
> do 
> strainName=$(echo $line|cut -d'_' -f3 |cut -d'.' -f1)
> prokka $line --kingdom Bacteria --genus Streptococcus --species agalactie --strain $strainName --usegenus --addgenes --prefix Streptococcus_agalactie_${strainName}\.prokka --outdir ~/gm_workshop/results/annotated/Streptococcus_agalactie_${strainName}\_prokka
> done
~~~
{: .language-bash}

~~~
$ mv ~/gm_workshop/results/annotated/*/*gbk ~/gm_workshop/results/annotated/
$ ls ~/gm_workshop/results/annotated/
~~~
{: .language-bash}
~~~
Streptococcus_agalactiae_18RS21.gbk_prokka  Streptococcus_agalactiae_COH1.gbk_prokka  
Streptococcus_agalactiae_515.gbk_prokka     Streptococcus_agalactiae_H36B.gbk_prokka  
Streptococcus_agalactiae_A909.gbk_prokka    thermophilusLMD9_prokka  
Streptococcus_agalactiae_CJB111.gbk_prokka  thermophilusLMG18311_prokka
~~~
{: .output}


> ## Genome annotation services
> To know more about prokka you can read [Seemann T. 2014](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517). Other valuable web-based genome annotation services are [RAST](https://rast.nmpdr.org/) and [PATRIC](https://www.patricbrc.org/). Both provides a web-user interface where you can storage your private genomes and share them with your colleagues. To use RAST as a command-line tool you can try the docker container [myRAST](https://github.com/nselem/myrast). 
{: .callout}
