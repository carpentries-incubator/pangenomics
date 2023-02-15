---
title: "Annotating Genomic Data"
teaching: 30
exercises: 15
questions:
- "How can I identify the genes in a genome?"
objectives:
- "Annotate bacterial genomes using Prokka."
keypoints:
- "Prokka is a command line utility that provides rapid prokaryotic genome annotation."
---

## Annotating Genomes

[Annotation](https://en.wikipedia.org/wiki/DNA_annotation) is the process of
identifying the coordinates of genes and all the coding regions
in a genome and determining what proteins are produced from them. In order to do this, an unknown
sequence is enriched with information relating genomic position, regulatory
sequences, repeats, gene name and protein products. This information
is stored in genomic databases to help future analysis processing new data.

[Prokka](https://github.com/tseemann/prokka)
is a command-line software tool created in Perl to annotate bacterial,
archaeal and viral genomes and reproduce standards-compliant output files.
It requires preassembled genomic DNA sequences in FASTA format as input
file, which is the only mandatory parameter to the software.
For annotation, Prokka relies on external features and databases to
identify the genomic features within the contigs.

| Tool (reference) | Features predicted |
| --------- | ----------- |
|Prodigal (Hyatt 2010 )   | Coding Sequence (CDS) |
| RNAmmer ( Lagesen et al. , 2007 )  | Ribosomal RNA genes (rRNA) |
| Aragorn ( Laslett and Canback, 2004 )  | Transfer RNA genes |
| SignalP ( Petersen et al. , 2011 )  | Signal leader peptides|
| Infernal ( Kolbe and Eddy, 2011 )  | Non-coding RNA|

Protein coding genes are annotated in two stages. Prodigal identifies
the coordinates of candidate genes, but does not describe the putative
gene product. Usually, in order to predict what a gene encodes
for, it is compared with a large
database of known sequences, usually at the protein level,
and transferred the annotation of the best significant match.
Prokka uses this method, but in a hierarchical manner. It starts
with a small trustworthy database, it then moves to medium
sized but domain specific databases and finally to curated
models of protein families.

> ## Notes
> [Environment variables](https://opensource.com/article/19/8/what-are-environment-variables) are special variables that contain information about your login session. These can be useful when you want to [manage](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#updating-an-environment) default or new settings that your system usually ignores. For instance, you can download a specific package with a downgraded version of perl if needed.
{: .callout}

First, we need to create a new directory where our annotated genomes will be.

~~~
$ mkdir -p ~/gm_workshop/results/annotated/
$ cd ~/gm_workshop/results/annotated/
$ conda deactivate
$ conda activate Prokka_Global
~~~
{: .language-bash}

Once inside the environment, we are ready to run our first annotation.  
~~~
(Prokka_Global) $
~~~
{: .language-bash}

In this example, we will use the following options:

| Code   | Meaning |
| ------- | ---------- |
| --prefix | Filename output prefix [auto] (default '') |
| --outdir | Output folder [auto] (default '') |
| --kingdom | Annotation mode: Archaea Bacteria Mitochondria Viruses (default 'Bacteria')|
| --genus | Genus name (default 'Genus') |
| --strain | Strain name (default 'strain') |
| --usegenus | Use genus-specific BLAST databases (needs --genus) (default OFF) |
| --addgens |Add 'gene' features for each 'CDS' feature (default OFF) |

~~~
$ prokka --prefix thermophilusLMD9.prokka --outdir thermophilusLMD9_prokka --kingdom Bacteria --genus Streptococcus --strain LMD9 --usegenus --addgenes ~/gm_workshop/data/thermophilusLMD9/GCF_000014485.1_ASM1448v1_genomic.fna
~~~
{: .language-bash}

Now Prokka has generated a new folder. If you run the `tree` command inside the new directory, you can view the set of files created by Prokka:

~~~
.
agalactiae_A909_prokka/
├── agalactiae_A909.prokka.err
├── agalactiae_A909.prokka.faa
├── agalactiae_A909.prokka.ffn
├── agalactiae_A909.prokka.fna
├── agalactiae_A909.prokka.fsa
├── agalactiae_A909.prokka.gbk
├── agalactiae_A909.prokka.gff
├── agalactiae_A909.prokka.log
├── agalactiae_A909.prokka.sqn
├── agalactiae_A909.prokka.tbl
├── agalactiae_A909.prokka.tsv
└── agalactiae_A909.prokka.txt

0 directories, 12 files
~~~
{: .output}

We encourage you to explore each output. The following table describes the contents of each output file:

| Extension | Description |
| --------- | ----------- |
| .gff | This is the master annotation in GFF3 format, containing both sequences and annotations. It can be viewed directly in Artemis or IGV. |
| .gbk | This is a standard GenBank file derived from the master `.gff`. If the input to Prokka was a multi-FASTA, then this will be a multi-GenBank, with one record for each sequence. |
| .fna | Nucleotide FASTA file of the input contig sequences. |
| .faa | Protein FASTA file of the translated CDS sequences. |
| .ffn | Nucleotide FASTA file of all the prediction transcripts (CDS, rRNA, tRNA, tmRNA, misc_RNA). |
| .sqn | An ASN1 format "Sequin" file for submission to GenBank. It needs to be edited to set the correct taxonomy, authors, related publication etc. |
| .fsa | Nucleotide FASTA file of the input contig sequences, used by "tbl2asn" to create the `.sqn` file. It is almost the same as the `.fna` file, but with extra Sequin tags in the sequence description lines. |
| .tbl | Feature Table file, used by "tbl2asn" to create the `.sqn` file. |
| .err | Unacceptable annotations - the NCBI discrepancy report. |
| .log | Contains all the output that Prokka produced during its run. This is the record of the used settings, even if the `--quiet` option was enabled. |
| .txt | Statistics related to the found annotated features. |
| .tsv | Tab-separated file of all features: locus_tag,ftype,len_bp,gene,EC_number,COG,product. |

Parameters can be modified as much as needed regarding the organism, the gene and even the locus tag you are looking for.

> ## Exercise 1: Extracting tRNAs with Prokka
>
> Suppose you are asked to annotate the FASTA file you downloaded in the previous episode and
> output the results to a subdirectory called `annotated` within the `thermophilusLMG18311_prokka`
> directory. Then, a research team requests you a TSV file named `trnas.tsv` that contains only
> *S. thermophilus'*s tRNAs. This file must contain the same headers as the original
> TSV file, followed by the rows that correspond to tRNAs.
>
> Complete the following sequence of commands to perform this actions:
>
> ~~~
> $ prokka --outdir thermophilusLMG18311_prokka --prefix thermophilusLMG18311.prokka ../../data/thermophilusLMG18311/__________ --kingdom Bacteria --genus Streptococcus --species thermophilus --usegenus --addgenes
> ~~~
> {: .language-bash}
>
> ~~~
> $ cd __________
> ~~~
> {: .language-bash}
>
> ~~~
> $ __________ -n 1 thermophilusLMG18311.prokka.tsv > trnas.tsv  # Get column headers
> $ grep __________ thermophilusLMG18311.prokka.tsv >> trnas.tsv # Append all lines that contain tRNAs
> ~~~
> {: .language-bash}
>
> > ## Solution
> >
> > Firstly, perform the annotation with Prokka and save all files as `thermophilusLMG18311_prokka`.
> >
> > ~~~
> > $ prokka --prefix thermophilusLMG18311.prokka --outdir thermophilusLMG18311_prokka --kingdom Bacteria --genus Streptococcus --strain LMG18311 --usegenus --addgenes ../../data/thermophilusLMG18311/GCF_000011825.1_ASM1182v1_genomic.fna
> > ~~~
> > {: .language-bash}
> >
> > After switching to the `thermophilusLMG18311_prokka` directory, we should filter the data we need and save the outputs to a file named `trnas.tsv`. To do so, we use the `head` command with the `-n 1` argument to get the first line (the headers of the columns). Next, we add the lines that correspond to tRNAs, which is done with the code `$'\t'tRNA$'\t'`(this means that the program will search for lines that contain the word `tRNA` with tab spaces at the beginning and the end of the word).
> >
> > ~~~
> > $ cd thermophilusLMG18311_prokka
> > $ head -n 1 thermophilusLMG18311.prokka.tsv > trnas.tsv # Get column headers
> > $ grep $'\t'tRNA$'\t' thermophilusLMG18311.prokka.tsv >> trnas.tsv # Append all lines that contain tRNA
> > ~~~
> > {: .language-bash}
> >
> {: .solution}
{: .challenge}

> ## Discussion 1: Number of tRNAs
>
> Inside the `annotated` directory, run the command `wc -l trnas.tsv` to get the number of lines in the file. Excluding the first line (which contains the header), observe that there are 67 tRNAs, whereas the [standard codon table](https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables#Translation_table_1) contains 61 coding codons (i.e. those that are non stop codons). This means that there are several genes that produce the same tRNA. How could you explain this fact?
>
> > ## Solution
> >
> > The existence of many genes producing the same product (such as the same tRNA) can happen due to duplications of genes during the evolutionary history of a taxon.
> {: .solution}
{: .discussion}

## Annotating multiple genomes

Now that we know how to annotate genomes with Prokka we can annotate all of
the *S. agalactiae* in one run.
For this purpose we will use a complex `while` loop that, for each of the *S. agalactiae* genomes,
first extracts the strain name and saves it in a variable, and then uses it inside the
Prokka command.

~~~
$ cd ~/gm_workshop/data
$ ls */*fasta | while read line
> do
> strainName=$(echo $line|cut -d'_' -f3 |cut -d'.' -f1)
> prokka $line --kingdom Bacteria --genus Streptococcus --species agalactiae --strain $strainName --usegenus --addgenes --prefix Streptococcus_agalactiae_${strainName}\.prokka --outdir ~/gm_workshop/results/annotated/Streptococcus_agalactiae_${strainName}\_prokka
> done
~~~
{: .language-bash}

Since we are only using the `.gbk` files, we will move them to the `results/annotated/` directory and remove the subdirectories.
~~~
$ cd ~/gm_workshop/results/annotated/
$ mv */*gbk .
$ rm -r *_prokka
$ ls
~~~
{: .language-bash}
~~~
Streptococcus_agalactiae_18RS21.prokka.gbk  Streptococcus_agalactiae_CJB111.prokka.gbk  thermophilusLMD9.prokka.gbk
Streptococcus_agalactiae_515.prokka.gbk 	Streptococcus_agalactiae_COH1.prokka.gbk	thermophilusLMG18311.prokka.gbk
Streptococcus_agalactiae_A909.prokka.gbk	Streptococcus_agalactiae_H36B.prokka.gbk
~~~
{: .output}


> ## Genome annotation services
> To learn more about Prokka you can read [Seemann T. 2014](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517). Other valuable web-based genome annotation services are [RAST](https://rast.nmpdr.org/) and [PATRIC](https://www.patricbrc.org/). Both provide a web-user interface where you can store your private genomes and share them with your colleagues. If you want to use RAST as a command-line tool you can try the docker container [myRAST](https://github.com/nselem/myrast).
{: .callout}



