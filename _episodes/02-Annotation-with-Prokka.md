---
title: "Downloading and Annotating Genomic Data"
teaching: 35
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
ncbi-genome-download --help
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
  positional arguments:  
    groups              The NCBI taxonomic groups to download (default: all).
                        A comma-separated list of taxonomic groups is also
                        possible. For example: "bacteria,viral"Choose from:
                        ['all', 'archaea', 'bacteria', 'fungi',
                        'invertebrate', 'metagenomes', 'plant', 'protozoa',
                        'vertebrate_mammalian', 'vertebrate_other', 'viral']  
    optional arguments:  
    -h, --help            show this help message and exit  
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
    -l ASSEMBLY_LEVELS, --assembly-levels ASSEMBLY_LEVELS  
                        Assembly levels of genomes to download (default: all).  
                        A comma-separated list of assembly levels is also
                        possible. For example: "complete,chromosome". Choose
                        from: ['all', 'complete', 'chromosome', 'scaffold',
                        'contig']  
    -g GENERA, --genera GENERA  
                        Only download sequences of the provided genera. A
                        comma-seperated list of genera is also possible. For
                        example: "Streptomyces coelicolor,Escherichia coli".
                        (default: [])  
    --genus GENERA        Deprecated alias of --genera  
    --fuzzy-genus         Use a fuzzy search on the organism name instead of an
                         exact match.  
    -S STRAINS, --strains STRAINS   
                        Only download sequences of the given strain(s). A
                        comma-separated list of strain names is possible, as
                        well as a path to a filename containing one name per
                        line.  
    -T SPECIES_TAXIDS, --species-taxids SPECIES_TAXIDS  
                        Only download sequences of the provided species NCBI
                        taxonomy IDs. A comma-separated list of species taxids
                        is also possible. For example: "52342,12325".
                        (default: [])  
    -t TAXIDS, --taxids TAXIDS  
                        Only download sequences of the provided NCBI taxonomy
                        IDs. A comma-separated list of taxids is also
                        possible. For example: "9606,9685". (default: [])  
    -A ASSEMBLY_ACCESSIONS, --assembly-accessions ASSEMBLY_ACCESSIONS  
                        Only download sequences matching the provided NCBI
                        assembly accession(s). A comma-separated list of
                        accessions is possible, as well as a path to a
                        filename containing one accession per line.  
    -R REFSEQ_CATEGORIES, --refseq-categories REFSEQ_CATEGORIES  
                        Only download sequences of the provided refseq
                        categories (default: all)  
    --refseq-category REFSEQ_CATEGORIES  
                        Deprecated alias for --refseq-categories  
    -o OUTPUT, --output-folder OUTPUT   
                        Create output hierarchy in specified folder (default:
                        /home/betterlab)  
    --flat-output         Dump all files right into the output folder without
                        creating any subfolders.  
    -H, --human-readable  Create links in human-readable hierarchy (might fail
                        on Windows)  
    -P, --progress-bar    Create a progress bar for indicating the download
                        progress  
    -u URI, --uri URI     NCBI base URI to use (default:
                        https://ftp.ncbi.nih.gov/genomes)  
    -p N, --parallel N    Run N downloads in parallel (default: 1)  
    -r N, --retries N     Retry download N times when connection to NCBI fails
                        (default: 0)  
    -m METADATA_TABLE, --metadata-table METADATA_TABLE  
                        Save tab-delimited file with genome metadata  
    -n, --dry-run         Only check which files to download, don't download
                        genome files.  
    -N, --no-cache        Don't cache the assembly summary file in
                        /home/betterlab/.cache/ncbi-genome-download.  
    -v, --verbose         increase output verbosity  
    -d, --debug           print debugging information   
    -V, --version         print version information  
    -M TYPE_MATERIALS, --type-materials TYPE_MATERIALS  
                        Specifies the relation to type material for the
                        assembly (default: any). "any" will include assemblies
                        with no relation to type material value defined, "all"
                        will download only assemblies with a defined value. A
                        comma-separated list of relatons. For example:
                        "reference,synonym". Choose from: ['any', 'all',
                        'type', 'reference', 'synonym', 'proxytype',
                        'neotype'].  
~~~
{: .output}

Once we know about the flags we can use, we are ready to use the package but first we need to go to our data directory.

~~~
cd dc_workshop/data
~~~
{: .bash}

If you list the contents of this directory (using the `ls` command), you'll see several folders, each of which contain the raw data of different strains of Streptococcus agalactiae used in Tettelin et al (2005). We will add a new entry to this directory by downloading a FASTA file of the genomes of all Streptococcus constellatus strains available on NCBI:

~~~
ncbi-genome-download --formats fasta  --genera "Streptococcus constellatus" bacteria
~~~
{: .bash}

The above code will display a new Genbank directory with subdirectories named after their NCBI assembly number.

The `gbk` files contained within this folder are used for the posterior analysis. It provides information about the coding sequences, their locus, the name of the protein, and the full nucleotide sequence of the assembly. By setting the `--formats` parameter to `genbank`, you tell the script to download `gbk` files.

> ## Exercise 1: Downloading data from NCBI with the command line
>
> Using `ncbi-genome-download`, get a FASTA file of the ATCC 31377 strain of the Streptococcus ratti bacterium and save it to an output directory titled `ratti`. Then, unzip the `gz` file and move the FASTA file all they back to the `ratti` directory. After you've done that, delete the `refseq` directory.
>
> > ## Solution
> >
> > First, we run the download utility.
> >
> > ~~~
> > ncbi-genome-download -F fasta -g "Streptococcus ratti" -S "ATCC 31377" -o ratti bacteria
> > ~~~
> > {: .source}
> >
> > Next, we navigate to the downloaded `gz` file and unzip it.
> >
> > ~~~
> > cd ratti/refseq/bacteria/GCF_008803015.1/
> > gunzip GCF_008803015.1_ASM880301v1_genomic.fna.gz
> > ~~~
> > {: .source}
> >
> > Then, we move the unzipped FASTA file to the `ratti` directory.
> >
> > ~~~
> > mv GCF_008803015.1_ASM880301v1_genomic.fna ../../..
> > ~~~
> > {: .source}
> >
> > Finally, we delete the `refseq` directory.
> >
> > ~~~
> > cd ../../..
> > rm -rf refseq
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

To get started with Prokka, we must first install it, [here](https://czirion.github.io/comparative-genomics-workshop/setup.html) you can find the set up page. 

Next, we need to change to the directory where we have the assembly (FASTA) files of interest. As a simple initial example of execution, we can annotate a FASTA file and define names for our output directory and files like this:

~~~
prokka example.fasta --outdir exdir --prefix exf
~~~
{: .source}

If you run the `tree` command in the current directory, you can preview the system of files created by Prokka:

~~~
exdir/
├── exf.err
├── exf.faa
├── exf.ffn
├── exf.fna
├── exf.fsa
├── exf.gbk
├── exf.gff
├── exf.log
├── exf.sqn
├── exf.tbl
├── exf.tsv
└── exf.txt
~~~
{: .output}

The following table describes the contents of each output file:

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

You can also add further details regarding the organism you search and the way the files will be annotated. For instance, if you'd like to annotate an archaeon of the genus *Nitrososphaera*, you would execute the following command:

~~~
prokka example.fasta --kingdom Archaea --genus Nitrososphaera --outdir exdir
~~~
{: .source}

You may specify your queries as much as you like. Type `prokka --help` in the command line to get the complete list of available parameters:

~~~
General:
  --help            This help
  --version         Print version and exit
  --citation        Print citation for referencing Prokka
  --quiet           No screen output (default OFF)
  --debug           Debug mode: keep all temporary files (default OFF)
Setup:
  --listdb          List all configured databases
  --setupdb         Index all installed databases
  --cleandb         Remove all database indices
  --depends         List all software dependencies
Outputs:
  --outdir [X]      Output folder [auto] (default '')
  --force           Force overwriting existing output folder (default OFF)
  --prefix [X]      Filename output prefix [auto] (default '')
  --addgenes        Add 'gene' features for each 'CDS' feature (default OFF)
  --locustag [X]    Locus tag prefix (default 'PROKKA')
  --increment [N]   Locus tag counter increment (default '1')
  --gffver [N]      GFF version (default '3')
  --compliant       Force Genbank/ENA/DDJB compliance: --genes --mincontiglen 200 --centre XXX (default OFF)
  --centre [X]      Sequencing centre ID. (default '')
Organism details:
  --genus [X]       Genus name (default 'Genus')
  --species [X]     Species name (default 'species')
  --strain [X]      Strain name (default 'strain')
  --plasmid [X]     Plasmid name or identifier (default '')
Annotations:
  --kingdom [X]     Annotation mode: Archaea|Bacteria|Mitochondria|Viruses (default 'Bacteria')
  --gcode [N]       Genetic code / Translation table (set if --kingdom is set) (default '0')
  --prodigaltf [X]  Prodigal training file (default '')
  --gram [X]        Gram: -/neg +/pos (default '')
  --usegenus        Use genus-specific BLAST databases (needs --genus) (default OFF)
  --proteins [X]    Fasta file of trusted proteins to first annotate from (default '')
  --hmms [X]        Trusted HMM to first annotate from (default '')
  --metagenome      Improve gene predictions for highly fragmented genomes (default OFF)
  --rawproduct      Do not clean up /product annotation (default OFF)
Computation:
  --fast            Fast mode - skip CDS /product searching (default OFF)
  --cpus [N]        Number of CPUs to use [0=all] (default '8')
  --mincontiglen [N] Minimum contig size [NCBI needs 200] (default '1')
  --evalue [n.n]    Similarity e-value cut-off (default '1e-06')
  --rfam            Enable searching for ncRNAs with Infernal+Rfam (SLOW!) (default '0')
  --norrna          Don't run rRNA search (default OFF)
  --notrna          Don't run tRNA search (default OFF)
  --rnammer         Prefer RNAmmer over Barrnap for rRNA prediction (default OFF)
~~~
{: .output}

Prokka can understand FASTA files with plain or detailed tag formats. A plain tag format looks like this:

~~~
>SeqID product
~~~
{: .output}

The detailed one consists of a special encoded three-part description line. The parts are the `/EC_number`, the `/gene` code, then the `/product` - and they are separated by a special "~~~" sequence:

~~~
>SeqID EC_number~~~gene~~~product~~~COG
~~~
{: .output}

Here are some examples. Note that not all parts need to be present, but the "~~~" should still be there:

~~~
>YP_492693.1 2.1.1.48~~~ermC~~~rRNA adenine N-6-methyltransferase~~~COG1234
MNEKNIKHSQNFITSKHNIDKIMTNIRLNEHDNIFEIGSGKGHFTLELVQRCNFVTAIEI
DHKLCKTTENKLVDHDNFQVLNKDILQFKFPKNQSYKIFGNIPYNISTDIIRKIVF*
>YP_492697.1 ~~~traB~~~transfer complex protein TraB~~~
MIKKFSLTTVYVAFLSIVLSNITLGAENPGPKIEQGLQQVQTFLTGLIVAVGICAGVWIV
LKKLPGIDDPMVKNEMFRGVGMVLAGVAVGAALVWLVPWVYNLFQ*
>YP_492694.1 ~~~~~~transposase~~~
MNYFRYKQFNKDVITVAVGYYLRYALSYRDISEILRGRGVNVHHSTVYRWVQEYAPILYQ
QSINTAKNTLKGIECIYALYKKNRRSLQIYGFSPCHEISIMLAS*
~~~
{: .output}

> ## Exercise 2: Extracting tRNAs with Prokka
>
> Using Prokka, annotate the FASTA format assembly of the Streptococcus ratti ATCC 31377 strain that you have downloaded in Exercise 1. Then explore the output files, focus on the TSV one and extract the tRNAs in a new file.
>
> > ## Solution
> >
> > First, we perform the annotation with Prokka and save all files as `atcc31377` in a directory titled `atcc31377`.
> >
> > ~~~
> > prokka GCF_008803015.1_ASM880301v1_genomic.fna --outdir atcc31377 --prefix atcc31377
> > ~~~
> > {: .source}
> >
> > Now we must filter the data we need and save the outputs to a file named `trnas.tsv`
> >
> > ~~~
> > cd atcc31377
> > head -n 1 atcc31377.tsv > trnas.tsv # Get column headers
> > grep $'\t'tRNA$'\t' atcc31377.tsv >> trnas.tsv # Append all lines that contain tRNA
> > ~~~
> > {: .source}
> >
> {: .solution}
{: .challenge}

> ## Discussion 1: Number of tRNAs
> 
> Inside the `atcc31377` directory, run the command `wc -l trnas.tsv` to get the number of lines in the file. Excluding the first line (which contains the header), observe that there are 66 tRNAs, whereas the total number of coding codons (i.e. the one that are not stop codons) is 61. This means that there are several genes that produce the same tRNA. How would you explain this fact?
> 
> > ## Solution
> > 
> > The existence of many genes producing the same product (such as the same tRNA) can happen due to duplications of genes during the evolutionary history of a taxon.
> {: .solution}
{: .challenge}
