---
title: "Downloading and Annotating Genomic Data"
teaching: 30
exercises: 15
questions:
- "How to download NCBI genomic data from the command line?"
- "What is a quick way to annotate a FASTA file and obtain different outputs?"
objectives:
- "Explore ncbi-genome-download as a tool for genomic data fetching from the NCBI."
- "Learn how to use the Prokka genome annotation utility."
keypoints:
- "The `ncbi-genome-download` package is a set of scripts to download genomes from the NCBI FTP servers implemented in Python."
- "Prokka is a command line utility that provides rapid prokaryotic genome annotation written in Perl."
---

## The `ncbi-genome-download` package: Getting Genomic Data from the NCBI

The NCBI Genome Downloading Scripts provide a shell command that allows users to download bacterial and fungal genomes from the NCBI. This [package]((https://github.com/kblin)) is highly useful because we can access to the genomic databases in a very specific way. Then, it will be easier to get the data directly in our work directory.

The full list of parameters you can incorporate in your dowmloads can be obtained by typing:

~~~
$ ncbi-genome-download --help
~~~
{: .source}

~~~
usage:  
 ncbi-genome-download [-h] [-s {refseq,genbank}] [-F FILE_FORMATS]  
                          [-l ASSEMBLY_LEVELS] [-g GENERA] [--genus GENERA]  
                          [--fuzzy-genus] [-S STRAINS] [-T SPECIES_TAXIDS]  
                          [-t TAXIDS] [-A ASSEMBLY_ACCESSIONS]  
                          [-R REFSEQ_CATEGORIES]  
                          [--refseq-category REFSEQ_CATEGORIES] [-o OUTPUT]  
                          [--flat-output] [-H] [-P] [-u URI] [-p N] [-r N]  
                          [-m METADATA_TABLE] [-n] [-N] [-v] [-d] [-V]  
                          [-M TYPE_MATERIALS]
                          groups
    -s {refseq,genbank}, --section {refseq,genbank}  
                        NCBI section to download (default: refseq)  
    -F FILE_FORMATS, --formats FILE_FORMATS  
                        Which formats to download (default: genbank).A comma-
                        separated list of formats is also possible. For
                        example: "fasta,assembly-report". Choose from:
                        ['genbank', 'fasta', 'rm', 'features', 'gff',
                        'protein-fasta', 'genpept', 'wgs', 'cds-fasta', 'rna-
                        fna', 'rna-fasta', 'assembly-report', 'assembly-
                        stats', 'all'] 
    -g GENERA, --genera GENERA  
                        Only download sequences of the provided genera. A
                        comma-seperated list of genera is also possible. For
                        example: "Streptomyces coelicolor,Escherichia coli".
                        (default: [])  
    -S STRAINS, --strains STRAINS   
                        Only download sequences of the given strain(s). A
                        comma-separated list of strain names is possible, as
                        well as a path to a filename containing one name per
                        line.   
    -A ASSEMBLY_ACCESSIONS, --assembly-accessions ASSEMBLY_ACCESSIONS  
                        Only download sequences matching the provided NCBI
                        assembly accession(s). A comma-separated list of
                        accessions is possible, as well as a path to a
                        filename containing one accession per line.   
    -o OUTPUT, --output-folder OUTPUT   
                        Create output hierarchy in specified folder (default:
                        /home/betterlab)  
    -m METADATA_TABLE, --metadata-table METADATA_TABLE  
                        Save tab-delimited file with genome metadata  
    -n, --dry-run         Only check which files to download, don't download
                        genome files.  
~~~
{: .output}

Once we know about the flags we can use, we are ready to use the package but first we need to go to our data directory.

~~~
$ cd dc_workshop/data
~~~
{: .language-bash}

If you list the contents of this directory (using the `ls` command), you'll see several folders, each of which contain the raw data of different strains of *Streptococcus agalactiae* used in Tettelin *et al*., (2005) in `.gbk` and `.fasta` format. 

~~~
$ ls 
~~~
{: .language-bash}

~~~
18RS21/  515/  A909/  COH1/  CJB111/  H36B/ 
~~~
{: .output}

We are interested in download a FASTA format of *Streptococcus thermophilus* LMD-9 strain from on NCBI. The `-n` flag allow us to check for this assembly availability before data downloading.

~~~
$ ncbi-genome-download --formats fasta --genera "Streptococcus thermophilus" -S LMD-9 -n -o thermophilusLMD9 bacteria
~~~
{: .language-bash}

~~~
Considering the following 1 assemblies for download:
GCF_000014485.1 Streptococcus thermophilus LMD-9        LMD-9
~~~
{: .output}

We will add a new entry to this directory by downloading the available assembly in a directory named `thermophilusLMD9`:

~~~
$ ncbi-genome-download --formats fasta --genera "Streptococcus thermophilus" -S LMD-9 -o thermophilusLMD9 bacteria
~~~
{: .language-bash}

After this execution, we will find the newly created directory and, at the bottom of the folder, the compressed assembly of interest named after the NCBI assembly number. We can check this by:

~~~
$ cd thermophilusLMD9/refseq/bacteria/GCF_000014485.1/
$ ls
~~~
{: .language-bash}

~~~
GCF_000014485.1_ASM1448v1_genomic.fna.gz  MD5SUMS
~~~
{: .output}

For a posterior analysis of this sequence, we need to decompress it

~~~
$ gunzip thermophilusLMD9/refseq/bacteria/GCF_000014485.1/GCF_000014485.1_ASM1448v1_genomic.fna.gz
~~~
{: .language-bash}

Now we can explore the file and move it to the main directory `thermophilusLMD9`

```
$ mv GCF_008803015.1_ASM880301v1_genomic.fna ../../..
$ cd ../../..
$ rm -r refseq
```
{: .language-bash}


> ## Notes
> In this example, we used the format FASTA to fetch our assembly. Nevertheless, we can use the flag `--help` to switch for other formats. For example, the `gbk` format files contained information about the coding sequences, their locus, the name of the protein, and the full nucleotide sequence of the assembly. At the same time it is a format also used for double-check the annotation.
{: .callout}

> ## Exercise 1: Downloading data from NCBI with the command line
> 
> Suppose you are asked to perform the following sequence of actions:
> 
> 1. Download a FASTA format assembly of the *Streptococcus thermophilus* with the NCBI assembly number `GCF_000011825.1` and save it to an output directory titled `thermophilusLMG18311`.
> 2. Change to the directory of the FASTA file and unzip it.
> 3. Move the FASTA file all the way back to the `thermophilusLMG18311` directory.
> 4. Change to the `thermophilusLMG18311` directory and delete the `refseq` subdirectory created by the downloading tool.
> 
> Complete the following set of commands to perform the previous steps:
> 
> Step 1.
> 
> ~~~
> ncbi-genome-download -F __________ --genera __________ -A __________ -o __________ bacteria
> ~~~
> {: .source}
> 
> Step 2.
> 
> ~~~
> cd __________/refseq/bacteria/GCF_000011825.1/
> __________ GCF_000011825.1_ASM1182v1_genomic.fna.gz
> ~~~
> {: .source}
> 
> Step 3.
> 
> ~~~
> mv GCF_000011825.1_ASM1182v1_genomic.fna __________
> ~~~
> {: .source}
> 
> Step 4.
> 
> ~~~
> cd __________
> rm -rf refseq
> ~~~
> {: .source}
>
> > ## Solution
> >
> > Step 1. Using the information provided in the step, we completed each blank space with the corresponding word:
> >
> > ~~~
> > ncbi-genome-download -F fasta --genera "Streptococcus " -A GCF_000011825.1 -o thermophilusLMG18311 bacteria
> > ~~~
> > {: .source}
> >
> > Step 2. The previous command creates a `thermophilusLMG18311` subdirectory. To get to the FASTA file we must go through a sequence of subdirectories and then apply the `gunzip` command to unzip the FASTA file:
> >
> > ~~~
> > $ cd thermophilusLMG18311/refseq/bacteria/GCF_000011825.1/
> > $ gunzip GCF_000011825.1_ASM1182v1_genomic.fna.gz
> > ~~~
> > {: .source}
> >
> > Step 3. We now have an unzipped FASTA file. The parent directory `thermophilusLMG18311` is located three directories above the current one. To get to the first parent directory, one would type `..`; if you want to get to the second parent directory, you would use `../..`. Thus, to move the FASTA file to the `thermophilusLMG18311` directory, we need to type:
> >
> > ~~~
> > $ mv GCF_000011825.1_ASM1182v1_genomic.fna ../../..
> > ~~~
> > {: .source}
> >
> > Step 4. Finally, we move back to the `thermophilusLMG18311` directory (in a similar manner as in the previous step) and delete the `refseq` directory.
> >
> > ~~~
> > $ cd ../../..
> > $ rm -rf refseq
> > ~~~
> > {: .source}
> >
> {: .solution}
{: .challenge}

## Prokka: Annotating Genomes

Annotation is a process of identifying the locations of genes and all the coding regions in a genome and determining
what those genes are for. For this, an unknown sequence is enriched with information relating genomic position, regulatory
sequences, repeats, gene name and protein products [1](https://en.wikipedia.org/wiki/DNA_annotation). This information
is stored in genomic databases to help future analysis processing new data.

Prokka is a command-line software tool created in Perl to annotate bacterial, archaeal and viral genomes and reproduce standards-compliant output files[2](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517?login=false).
It expects a preassembled genomic DNA sequences in FASTA format as the input file, which is the only mandatory parameter to the software.
For annotation, Prokka relies on external features and databases to identify the genomic features within the contigs.

| Tool(reference) | Features predicted |
| --------- | ----------- |
|Prodigal (Hyatt 2010 )   | Coding Sequence (CDS) |
| RNAmmer ( Lagesen et al. , 2007 )  | Ribosomla RNA genes (rRNA) |
| Aragorn ( Laslett and Canback, 2004 )  | Transfer RNA genes |
| SignalP ( Petersen et al. , 2011 )  | Signal leader peptides|
| Infernal ( Kolbe and Eddy, 2011 )  | Non-coding RNA|

Proteins coding genes are annotated in two stages. Prodigal identifies the coordinates of candidate genes, but does not
describe the putative gene product. The traditional way to predict what a gene codes for is to compare it with a large
database of known sequences, usually at a protein level, and transfer the annotation of the best significant match.
Prokka uses this method, but in a hierarchical manner, starting with a smaller trustworthy database, moving to medium
sized but domain specific databases and finally to curated models of protein families.

Next, we need to change to the directory where we have the assembly (FASTA) files of interest. As a simple initial example of execution, we can annotate a FASTA file and define names for our output directory and files like this:

~~~
$ cd ../results/annotated/
$ prokka --prefix thermophilusLMD9_prokka --outdir thermophilusLMD9_prokka --kingdom Bacteria --genus Streptococcus --strain LMD9 --usegenus --addgenes ../../thermophilusLMD9/GCF_000014485.1_ASM1448v1_genomic.fna
~~~
{: .source}

Now prokka has generated a new folder. Lets get in and if you run the `tree` command in the current directory, you can preview the system of files created by Prokka:

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
l

> ## Exercise 2: Extracting tRNAs with Prokka
>
> Suppose you are now asked to annotate the FASTA file you downloaded in Exercise 1 and output the results to a subdirectory called `annotated` within the `ratti` directory. Then, a research team asks you to provide them a TSV file titled `trnas.tsv` that only contains *S. ratti'*s tRNAs. This file must contain the same headers as the original TSV file, followed by the rows that correspond to tRNAs. Complete the following sequence of commands to perform this actions:
> 
> ~~~
> $ prokka --outdir thermophilusLMG18311_prokka --prefix thermophilusLMG18311_prokka ../../thermophilusLMG18311/__________ --kingdom Bacteria --genus Streptococcus --species thermophilus --usegenus --addgenes 
> ~~~
> {: .source}
> 
> ~~~
> $ cd __________
> ~~~
> {: .source}
> 
> ~~~
> $ __________ -n 1 thermophilusLMG18311_prokka.tsv > trnas.tsv  # Get column headers
> $ grep __________ thermophilusLMG18311_prokka.tsv >> trnas.tsv # Append all lines that contain tRNAs
> ~~~
> {: .source}
>
> > ## Solution
> >
> > First, we perform the annotation with Prokka and save all files as `thermophilusLMG18311_prokka`.
> >
> > ~~~
> > $ prokka --prefix thermophilusLMD9_prokka --outdir thermophilusLMD9_prokka --kingdom Bacteria --genus Streptococcus --strain LMD9 --usegenus --addgenes 
../../thermophilusLMD9/GCF_000014485.1_ASM1448v1_genomic.fna
> > ~~~
> > {: .source}
> >
> > After switching to the `thermophilusLMD9_prokka` directory, we shall now filter the data we need and save the outputs to a file named `trnas.tsv`. To do so, we use the `head` command with the `-n 1` argument to get the first line (the headers of the columns). We then append the lines that correspond to tRNAs, which is done with the code `$'\t'tRNA$'\t'`(this means that the program will search for lines that contain the word `tRNA` with tab spaces at the beginning and the end of the word).
> >
> > ~~~
> > cd annotated
> > head -n 1 thermophilusLMD9_prokka.tsv > trnas.tsv # Get column headers
> > grep $'\t'tRNA$'\t' thermophilusLMD9_prokka.tsv >> trnas.tsv # Append all lines that contain tRNA
> > ~~~
> > {: .source}
> >
> {: .solution}
{: .challenge}

> ## Discussion 1: Number of tRNAs
> 
> Inside the `annotated` directory, run the command `wc -l trnas.tsv` to get the number of lines in the file. Excluding the first line (which contains the header), observe that there are 66 tRNAs, whereas the [standard codon table](https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables#Translation_table_1) contains 61 coding codons (i.e. the ones that are not stop codons). This means that there are several genes that produce the same tRNA. How would you explain this fact?
> 
> > ## Solution
> > 
> > The existence of many genes producing the same product (such as the same tRNA) can happen due to duplications of genes during the evolutionary history of a taxon.
> {: .solution}
{: .challenge}
