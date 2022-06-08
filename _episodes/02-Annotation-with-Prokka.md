---
title: "Download data and annotation with Prokka"
teaching: 0
exercises: 0
questions:
- "Available tools for NCBI downloading data"
- "What is a Prokka?"
- "How to install Prokka"
- "Which means 'annotation' and how can I perform one?"


objectives:
- "Explain how to download genomic data and perform an annotation analysis"
keypoints:
- 

---

Instalation of NCBI-genome-download package  

NCBI Genome Downloading are a script to download bacterial and fungal genomes from NCBI after they
restructured their FTP a while ago. You can call the package typing:  

```
ncbi-genome-download --help
```
and you will get a full list of parameters you can incorporate for your downloading and use as a guide 
to use the package.  

usage: ncbi-genome-download [-h] [-s {refseq,genbank}] [-F FILE_FORMATS]
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
  groups                The NCBI taxonomic groups to download (default: all).
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
                        'neotype'] .  

Once we know about the flags we can use, we are ready to use the package. Then our first genome download
is going to be for the strains of Streptococcus equinus in FASTA format:

~~~
ncbi-genome-download --formats fasta .""--section genbank"". --genera "Streptococcus equinus" bacteria
~~~
{: .bash}

~~~

~~~
{: .output}

The above code will display a new genbank folder with the folders named with their assembly NCBI number.
We need then to extract these files from each folder, uncompress them and rename the files in a more descriptive way, 
like species and strain name. For these, we can use a for loop.

~~~
for
~~~
{: .output}

Then we can explore the gbk files. It is important to know these files because they will be used for the posterior
analysis. It will give you information about the coding sequences, their locus, the name of 
the protein, and the full nucleotide sequence of the assembly.

------------------Image


# Prokka

Annotation is a process of identifying the locations of genes and all the coding regions in a genome and determining 
what those genes do. For this, an unknown sequence is enriched with information relating genomic position, regulatory 
sequences, repeats, gene name and protein products [2](https://en.wikipedia.org/wiki/DNA_annotation). This information
is stored in genomic databases to help future analysis processing new data.

Prokka is a command-line software tool to annotate bacterial, archaeal and viral genomes and reproduce standards-compliant output files[2](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517?login=false).
It expects a preassembled genomic DNA sequences in FASTA format, which is the only mandatory parameter to the software.
For annotation, Prokka relies on external features and databases to identify the genomic features within the contigs.

| Tool(reference) | Features predicted |
| --------- | ----------- |
|Prodigal (Hyatt 2010 )   | Coding Sequence (CDS) |
| RNAmmer ( Lagesen et al. , 2007 )  | Ribosomla RNA genes (rRNA) |
| Aragorn ( Laslett and Canback, 2004 )  | Transfer RNA genes |
| SignalP ( Petersen et al. , 2011 )  | Signal leader peptides|
| Infernal ( Kolbe and Eddy, 2011 )  | Non-coding RNA|

Proteins coding genes are annotates in two stages. Prodigal identifies the coordinates of candidate genes, but does not 
describe the putative gene product. The traditional way to predict what a gene codes for is to compare it with a large
database of known sequences, usually at a protein level, and transfer the annotation of the best significant match. 
Prokka uses this method, but in a hierarchical manner, starting with a smaller trustworthy database, moving to medium 
sized but domain specific databases and finally to curated models of protein families.  

Beginning with Prokka, we need to set up on the folder where we have the assembly (FASTA) files of interest. As a simple 
example of execution, we can annotate a FASTA file and define names for our output directory and files like this:

~~~
prokka example.fasta --outdir exdir --prefix exf
~~~
{: .bash}

~~~
example.gff example.gbk example.fna example.faa example.ffn example.sqn
example.fsa example.tbl example.err example.log example.txt example.tsv
~~~
{: .output}




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


## Command line options

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
      
     
    The detailed one consists of a special encoded three-part description line. The parts are the `/EC_number`, 
the `/gene` code, then the `/product` - and they are separated by a special "~~~" sequence:
```
>SeqID EC_number~~~gene~~~product~~~COG
```

Here are some examples. Note that not all parts need to be present, but the "~~~" should still be there:
```
>YP_492693.1 2.1.1.48~~~ermC~~~rRNA adenine N-6-methyltransferase~~~COG1234
MNEKNIKHSQNFITSKHNIDKIMTNIRLNEHDNIFEIGSGKGHFTLELVQRCNFVTAIEI
DHKLCKTTENKLVDHDNFQVLNKDILQFKFPKNQSYKIFGNIPYNISTDIIRKIVF*
>YP_492697.1 ~~~traB~~~transfer complex protein TraB~~~
MIKKFSLTTVYVAFLSIVLSNITLGAENPGPKIEQGLQQVQTFLTGLIVAVGICAGVWIV
LKKLPGIDDPMVKNEMFRGVGMVLAGVAVGAALVWLVPWVYNLFQ*
>YP_492694.1 ~~~~~~transposase~~~
MNYFRYKQFNKDVITVAVGYYLRYALSYRDISEILRGRGVNVHHSTVYRWVQEYAPILYQ
QSINTAKNTLKGIECIYALYKKNRRSLQIYGFSPCHEISIMLAS*
```
