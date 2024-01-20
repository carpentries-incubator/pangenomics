---
title: "Annotating Genomic Data"
teaching: 30
exercises: 15
questions:
- "How can I identify the genes in a genome?"
objectives:
- "Annotate bacterial genomes using Prokka."
- "Use scripts to customize output files."
- "Identify genes conferring antibiotic resistance with RGI."
keypoints:
- "Prokka is a command line utility that provides rapid prokaryotic genome annotation."
- "Sometimes we need manual curation of the output files of the software."
- "Specialized software exist to perform annotation of specific genomic elements."
---

## Annotating genomes

[Annotation](https://en.wikipedia.org/wiki/DNA_annotation) is the process of
identifying the coordinates of genes and all the coding regions
in a genome and determining what proteins are produced from them. In order to do this, an unknown
sequence is enriched with information relating to genomic position, regulatory
sequences, repeats, gene names, and protein products. This information
is stored in genomic databases to help future analysis processing of new data.

[Prokka](https://github.com/tseemann/prokka)
is a command-line software tool created in Perl to annotate bacterial,
archaeal and viral genomes and reproduce standards-compliant output files.
It requires preassembled genomic DNA sequences in FASTA format as input
file, which is the only mandatory parameter to the software.
For annotation, Prokka relies on external features and databases to
identify the genomic features within the contigs.

| Tool (reference) | Features predicted |
| --------- | ----------- |
|Prodigal (Hyatt 2010)   | Coding Sequences (CDS) |
| RNAmmer (Lagesen et al., 2007)  | Ribosomal RNA genes (rRNA) |
| Aragorn (Laslett and Canback, 2004)  | Transfer RNA genes |
| SignalP (Petersen et al., 2011)  | Signal leader peptides|
| Infernal (Kolbe and Eddy, 2011)  | Non-coding RNAs|

Protein coding genes are annotated in two stages. Prodigal identifies
the coordinates of candidate genes, but does not describe the putative
gene product. Usually, in order to predict what a gene encodes
for, it is compared with a large
database of known sequences, usually at the protein level,
and transferred the annotation of the best significant match.
Prokka uses this method but in a hierarchical manner. It starts
with a small trustworthy database, it then moves to medium
sized but domain-specific databases and finally curated
models of protein families.

In this lesson, we'll annotate the `FASTA` files we have downloaded in the
previous lesson. First, we need to create a new directory where our annotated
genomes will be.

~~~
$ mkdir -p ~/pan_workshop/results/annotated/
$ cd ~/pan_workshop/results/annotated/
$ conda deactivate
$ conda activate /miniconda3/envs/Prokka_Global
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
$ prokka --prefix agalactiae_515_prokka --outdir agalactiae_515_prokka --kingdom Bacteria --genus Streptococcus --species agalactiae --strain 515 --usegenus --addgenes ~/pan_workshop/data/agalactiae_515/GCF_012593885.1_ASM1259388v1_genomic.fna
~~~
{: .language-bash}

This command takes about a minute to run, printing a lot of information on the screen while doing so. After finishing, Prokka will create a new folder, inside of which, if you run the `tree` command, you will find the following files:
~~~
tree
~~~
{: .language-bash}

~~~
.
└── agalactiae_515_prokka
    ├── agalactiae_515_prokka.err
    ├── agalactiae_515_prokka.faa
    ├── agalactiae_515_prokka.ffn
    ├── agalactiae_515_prokka.fna
    ├── agalactiae_515_prokka.fsa
    ├── agalactiae_515_prokka.gbk
    ├── agalactiae_515_prokka.gff
    ├── agalactiae_515_prokka.log
    ├── agalactiae_515_prokka.sqn
    ├── agalactiae_515_prokka.tbl
    ├── agalactiae_515_prokka.tsv
    └── agalactiae_515_prokka.txt

1 directory, 12 files
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
| .sqn | An ASN1 format "Sequin" file for submission to GenBank. It needs to be edited to set the correct taxonomy, authors, related publications etc. |
| .fsa | Nucleotide FASTA file of the input contig sequences, used by "tbl2asn" to create the `.sqn` file. It is almost the same as the `.fna` file but with extra Sequin tags in the sequence description lines. |
| .tbl | Feature Table file, used by "tbl2asn" to create the `.sqn` file. |
| .err | Unacceptable annotations - the NCBI discrepancy report. |
| .log | Contains all the output that Prokka produced during its run. This is the record of the used settings, even if the `--quiet` option was enabled. |
| .txt | Statistics related to the found annotated features. |
| .tsv | Tab-separated file of all features: locus_tag, ftype, len_bp, gene, EC_number, COG, product. |

Parameters can be modified as much as needed regarding the organism, the gene, and even the locus tag you are looking for.

> ## Exercise 1: Inspecting the GBK
> Open the `gbk` output file and carefully explore the information it contains. Which of the following statements is TRUE?
> 
> a) Prokka translates every single gene to its corresponding protein, even if the gene isn't a coding one.  
> b) Prokka can find all kinds of protein-coding sequences, not just the ones that have been identified or cataloged in a database.  
> c) Prokka identifies tRNA genes but doesn't mention the anticodon located on the tRNAs.  
> d) Prokka doesn't provide the positions in which a feature starts or ends.  
> e) The coding sequences are identified with the CDS acronym in the `FEATURES` section of each `LOCUS`.  
> 
>> ## Solution
>>  
>> a) FALSE. Prokka successfully identifies non-coding sequences and doesn't translate them. Instead, it provides alternative information (e.g. if it's a rRNA gene, it tells if it's 5S, 16S, or 23S).  
>> b) TRUE. Some coding sequences produce proteins that are marked as "hypothetical", meaning that they haven't been yet identified but seem to show properties of a coding sequence.  
>> c) FALSE. Every tRNA feature has a `/note` subsection mentioning between parentheses the anticodon located on the tRNA.  
>> d) FALSE. Right next to each feature, there's a pair of numbers indicating the starting and ending position of the corresponding feature.  
>> e) TRUE. Each coding sequence is identified by the CDS acronym on the left and information such as coordinates, gene name, locus tag, 
>> product description and translation on the right.
>{: .solution}
{: .challenge}

## Annotating multiple genomes

Now that we know how to annotate genomes with Prokka we can annotate all of
the *S. agalactiae* in one run.
For this purpose, we will use a complex `while` loop that, for each of the *S. agalactiae* genomes,
first extracts the strain name and saves it in a variable, and then uses it inside the
Prokka command.

To get the strain names easily we will update our `TettelinList.txt` to add the strain names 
that it does not have and change the problematic name of the strain 2603V/R.
We could just open the file in nano and edit it, but we will do it by coding the changes.  
With `echo` we will add each strain name in a new line, and with `sed` we will remove the 
characters `/R` of the problematic strain name.

~~~
$ cd ~/pan_workshop/data/
$ echo "18RS21" >> TettelinList.txt 
$ echo "H36B" >> TettelinList.txt
$ echo "515" >> TettelinList.txt 
$ sed -i 's/\/R//' TettelinList.txt
$ cat TettelinList.txt
~~~
{: .language-bash}

~~~
A909  
COH1  
CJB111 
NEM316
2603V
18RS21
H36B
515
~~~
{: .output}

We can now run Prokka on each of these strains. Since the following command can take up to 8 minutes to run we will use a `screen` session to run it.
The screen session will not have the conda environment activated, so let's activate it again.
~~~
screen -R prokka
conda activate /miniconda3/envs/Prokka_Global
~~~
{: .language-bash}  
~~~
$ cat TettelinList.txt | while read line
do 
prokka agalactiae_$line/*.fna --kingdom Bacteria --genus Streptococcus --species agalactiae \
--strain $line --usegenus --addgenes --prefix Streptococcus_agalactiae_${line}_prokka \
--outdir ~/pan_workshop/results/annotated/Streptococcus_agalactiae_${line}_prokka
done
~~~
{: .language-bash}  
Click `Ctrl`+ `a` + `d` to detach from the session and wait until it finishes the run.


> ## Genome annotation services
> To learn more about Prokka you can read [Seemann T. 2014](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517). Other valuable web-based genome annotation services
> are [RAST](https://rast.nmpdr.org/) and [PATRIC](https://www.patricbrc.org/). Both provide a web-based user interface where you can store your private genomes and share them 4
> with your colleagues. If you want to use RAST as a command-line tool you can try the docker container [myRAST](https://github.com/nselem/myrast).
{: .callout}


## Curating Prokka output files

Now that we have our genome annotations, let's take a look at one of them. Fortunately, the `gbk` files are human-readable and we can 
look at a lot of the information in the first few lines:

~~~
$ cd ../results/annotated/
$ head Streptococcus_agalactiae_18RS21_prokka/Streptococcus_agalactiae_18RS21_prokka.gbk
~~~
{: .language-bash}

~~~
LOCUS       AAJO01000553.1           259 bp    DNA     linear       22-FEB-2023
DEFINITION  Streptococcus agalactiae strain 18RS21.
ACCESSION   
VERSION
KEYWORDS    .
SOURCE      Streptococcus agalactiae
  ORGANISM  Streptococcus agalactiae
            Unclassified.
COMMENT     Annotated using prokka 1.14.6 from
            https://github.com/tseemann/prokka.

~~~
{: .output}

We can see that in the `ORGANISM` field we have the word "Unclassified". If we compare it to the `gbk` file for the same strain, that came with the original data folder (which was obtained from the NCBI) we can see that the strain code should be there. 

~~~
$ head ../../data/agalactiae_18RS21/Streptococcus_agalactiae_18RS21.gbk
~~~
{: .language-bash}

~~~
LOCUS       AAJO01000169.1          2501 bp    DNA     linear   UNK 
DEFINITION  Streptococcus agalactiae 18RS21
ACCESSION   AAJO01000169.1
KEYWORDS    .
SOURCE      Streptococcus agalactiae 18RS21.
  ORGANISM  Streptococcus agalactiae 18RS21
            Bacteria; Terrabacteria group; Firmicutes; Bacilli;
            Lactobacillales; Streptococcaceae; Streptococcus; Streptococcus
            agalactiae.
FEATURES             Location/Qualifiers

~~~
{: .output}


This difference could be a problem since some bioinformatics programs could classify two different strains within the same "Unclassified" group. 
For this reason, Prokka's output files need to be corrected before moving forward with additional analyses.

To do this "manual" curation we will use the script `correct_gbk.sh`. Let's first make a directory for the scripts, and then use of nano text editor to create your file.
~~~
$ mkdir ../../scripts
$ nano ../../scripts/correct_gbk.sh
~~~
{: .language-bash}
Paste the following content in your script:
~~~
#This script will change the word Unclassified from the ORGANISM lines by that of the respective strain code.
# Usage: sh correct_gbk.sh <gbk-file-to-edit>
file=$1   # gbk file annotated with prokka
strain=$(grep -m 1 "DEFINITION" $file |cut -d " " -f6,7) # Create a variable with the columns 6 and 7 from the DEFINITION line.

sed -i '/ORGANISM/{N;s/\n//;}' $file # Put the ORGANISM field on a single line.

sed -i "s/\s*Unclassified./ $strain/" $file # Substitute the word "Unclassified" with the value of the strain variable.
~~~
{: .language-bash}

Press `Ctrl + X` to exit the text editor and save the changes. This script allows us to change the term "Unclassified." from the rows ORGANISM with 
that of the respective strain. 

Now, we need to run this script for all the `gbk` files:
~~~
$ ls */*.gbk | while read file
do 
bash ../../scripts/correct_gbk.sh $file
done
~~~
{: .language-bash}

Finally, let's view the result:
~~~
$ head Streptococcus_agalactiae_18RS21_prokka/Streptococcus_agalactiae_18RS21_prokka.gbk
~~~
{: .language-bash}

~~~
LOCUS       AAJO01000553.1           259 bp    DNA     linear       27-FEB-2023
DEFINITION  Streptococcus agalactiae strain 18RS21.
ACCESSION   
VERSION
KEYWORDS    .
SOURCE      Streptococcus agalactiae
  ORGANISM  Streptococcus agalactiae 18RS21.
COMMENT     Annotated using prokka 1.14.6 from
            https://github.com/tseemann/prokka.
FEATURES             Location/Qualifiers
~~~
{: .output}

Voilà! Our `gbk` files now have the strain code in the `ORGANISM` line.

> ## Exercise 2: Counting coding sequences
> 
> Before we build our pangenome it can be useful to take a quick look at how many coding sequences each of our genomes have. This way we can 
> know if they have a number close to the expected one (if we have some previous knowledge of our organism of study).
> 
> Use your `grep`, looping, and piping abilities to count the number of coding sequences in the `gff` files of each genome.  
> 
> Note: We will use the `gff` file because the `gbk` contains the aminoacid sequences, so it is possible that with the `grep` command 
> we find the string `CDS` in these sequences, and not only in the description of the features. The `gff` files also have the description of the
> features but in a different format.
> 
>
> > ## Solution
> > First inspect a `gff` file to see what you are working with. Open it with `nano` and scroll through the file to see its contents. 
> > ~~~
> > nano Streptococcus_agalactiae_18RS21_prokka/Streptococcus_agalactiae_18RS21_prokka.gff
> > ~~~
> > Now make a loop that goes through every `gff` finding and counting each line with the string "CDS".
> > ~~~ 
> > > for genome in */*.gff
> > > do 
> > > echo $genome # Print the name of the file
> > > grep "CDS" $genome | wc -l # Find the lines with the string "CDS" and pipe that to the command wc with the flag -l to count the lines
> > > done
> > ~~~
> > {: .language-bash}
> > ~~~
> > Streptococcus_agalactiae_18RS21_prokka/Streptococcus_agalactiae_18RS21_prokka.gff
> > 1960
> > Streptococcus_agalactiae_2603V_prokka/Streptococcus_agalactiae_2603V_prokka.gff
> > 2108
> > Streptococcus_agalactiae_515_prokka/Streptococcus_agalactiae_515_prokka.gff
> > 1963
> > Streptococcus_agalactiae_A909_prokka/Streptococcus_agalactiae_A909_prokka.gff
> > 2067
> > Streptococcus_agalactiae_CJB111_prokka/Streptococcus_agalactiae_CJB111_prokka.gff
> > 2044
> > Streptococcus_agalactiae_COH1_prokka/Streptococcus_agalactiae_COH1_prokka.gff
> > 1992
> > Streptococcus_agalactiae_H36B_prokka/Streptococcus_agalactiae_H36B_prokka.gff
> > 2166
> > Streptococcus_agalactiae_NEM316_prokka/Streptococcus_agalactiae_NEM316_prokka.gff
> > 2139
> > ~~~
> > {: .output}
> {: .solution}
{: .challenge}

> ## Annotating your assemblies
>
> If you work with your own assembled genomes, other problems may arise when annotating them. One likely problem is that the name of 
> your contigs is very long, and since Prokka will use those names as the LOCUS names, the LOCUS names may turn out problematic.  
> Example of contig name:
> ~~~
> NODE_1_length_44796_cov_57.856817
> ~~~
> Result of LOCUS name in `gbk` file:
> ~~~
> LOCUS       NODE_1_length_44796_cov_57.85681744796 bp   DNA linear
> ~~~
> Here the coverage and the length of the locus are fused, so this will give problems downstream in your analysis.  
> 
> The tool [anvi-script-reformat-fasta](https://anvio.org/help/main/programs/anvi-script-reformat-fasta/) can help you simplify the names 
> of your assemblies and do other magic, such as removing the small contigs or sequences with too many gaps.  
> ~~~
> anvi-script-reformat-fasta my_new_assembly.fasta -o my_reformated_assembly.fasta --simplify-names
> ~~~
> {: .language-bash}
> This will convert `>NODE_1_length_44796_cov_57.856817` into `>c_000000000001` and the LOCUS name into
> `LOCUS       c_000000000001         44796 bp    DNA     linear`.  
> Problem solved!
{: .callout}

## Annotating antibiotic resistance

Whereas Prokka is useful to identify all kinds of genomic elements, other more
specialized pipelines are also available. For example,
[antiSMASH](https://antismash.secondarymetabolites.org) searches genomes for
biosynthetic gene clusters, responsible for the production of secondary
metabolites. Another pipeline of interest is
[RGI](https://github.com/arpcard/rgi): the Resistance Gene Identifier. This
program allows users to predict genes and SNPs which confer antibiotic
resistance to an organism. It is a very complex piece of software subdivided
into several subprograms; RGI main, for instance, is used to annotate contigs,
and RGI bwt, on the other hand, annotates reads. In this lesson, we'll learn
how to use RGI main. To use it, first activate its virtual environment:

~~~
$ conda activate rgi
~~~
{: .language-bash}

You can type `rgi --help` to list all subprograms that RGI provides. In order
to get a manual of a specific subcommand, type `rgi [command] --help`,
replacing `[command]` with your subprogram of interest. Before you do anything
with RGI, however, you must download [CARD](https://card.mcmaster.ca/)
(the Comprehensive Antibiotic Resistance Database), which is used by RGI as
reference. To do so, we will use `wget` and one of the subcommands of RGI,
`rgi load`, as follows:

~~~
$ cd ~/pan_workshop/data/
$ wget -O card_archive.tar.gz https://card.mcmaster.ca/latest/data
$ tar -xf card_archive.tar.gz ./card.json
$ rgi load --local -i card.json
$ rm card_archive.tar.gz card.json
~~~
{: .language-bash}

After performing this sequence of commands, you'll find a directory called
`localDB/` in your current working directory. Its location and name are
extremely important: **you must always run RGI inside the parent directory of
`localDB/`** (which, in our case, is `~/pan_workshop/data/`), and **you shall
not rename `localDB/` to anything else**. RGI will fail if you don't follow
these rules.

As we'll be using RGI main, write `rgi main --help` and take a moment to read
through the help page. The parameters we'll be using in this lesson are:

- `-i` or `--input_sequence`. Sets the genomic sequence (in `fasta` or
`fasta.gz` format) we want to annotate.
- `-o` or `--output_file`. Specifies the basename for the two output files RGI
produces; for example, if you set this option to `outputs`, you'll get two
files: `outputs.json` and `outputs.txt`.
- `--include_loose`. When not using this option, RGI will only return
hits with strict boundaries; on the other hand, if provided, RGI will also
include hits with loose boundaries.
- `--local`. Tells RGI to use the database stored in `localDB/`.
- `--clean`. Removes temporary files created by RGI.

We are now going to create a new directory for RGI main's outputs:

~~~
$ mkdir -p ../results/resistomes/
~~~
{: .language-bash}

Next, let's see how we would find the resistance genes in the 18RS21 strain of
*S. agalactiae*:

~~~
$ rgi main --clean --local --include_loose \
> -i agalactiae_18RS21/Streptococcus_agalactiae_18RS21.fna \
> -o ../results/resistomes/agalactiae_18RS21
~~~
{: .language-bash}

Recall that RGI produces two output files; let's have a look at them: 

~~~
$ cd ../results/resistomes/
$ ls
~~~
{: .language-bash}

~~~
agalactiae_18RS21.json
agalactiae_18RS21.txt
~~~
{: .output}

The `JSON` file stores the complete output whereas the `.TXT` file contains a
subset of this information. However, the former isn't very human-readable, and
is mostly useful for downstream analyses with RGI; the latter, on the contrary,
has everything we might need in a "friendlier" format. This file is
tab-delimited, meaning that it is a table file which uses the tab symbol as
separator. Have a look at the file by running `less -S agalactiae_18RS21.txt`;
use the arrow keys to move left and right. A detailed description of the
meaning of each column can be found in the table below (taken from RGI's
documentation):

|Column|Field                                  | Contents                                       |
|:-----|:--------------------------------------|:-----------------------------------------------|
|   1  |ORF_ID                                 | Open Reading Frame identifier (internal to RGI)|
|   2  |Contig                                 | Source Sequence                                |
|   3  |Start                                  | Start co-ordinate of ORF                       |
|   4  |Stop                                   | End co-ordinate of ORF                         |
|   5  |Orientation                            | Strand of ORF                                  |
|   6  |Cut_Off                                | RGI Detection Paradigm (Perfect, Strict, Loose)|
|   7  |Pass_Bitscore                          | Strict detection model bitscore cut-off        |
|   8  |Best_Hit_Bitscore                      | Bitscore value of match to top hit in CARD     |
|   9  |Best_Hit_ARO                           | ARO term of top hit in CARD                    |
|  10  |Best_Identities                        | Percent identity of match to top hit in CARD   |
|  11  |ARO                                    | ARO accession of match to top hit in CARD      |
|  12  |Model_type                             | CARD detection model type                      |
|  13  |SNPs_in_Best_Hit_ARO                   | Mutations observed in the ARO term of top hit in CARD (if applicable)|
|  14  |Other_SNPs                             | Mutations observed in ARO terms of other hits indicated by model id (if applicable)|
|  15  |Drug Class                             | ARO Categorization                             |
|  16  |Resistance Mechanism                   | ARO Categorization                             |
|  17  |AMR Gene Family                        | ARO Categorization                             |
|  18  |Predicted_DNA                          | ORF predicted nucleotide sequence              |
|  19  |Predicted_Protein                      | ORF predicted protein sequence                 |
|  20  |CARD_Protein_Sequence                  | Protein sequence of top hit in CARD            |
|  21  |Percentage Length of Reference Sequence| (length of ORF protein / length of CARD reference protein)|
|  22  |ID                                     | HSP identifier (internal to RGI)               |
|  23  |Model_id                               | CARD detection model id                        |
|  24  |Nudged                                 | TRUE = Hit nudged from Loose to Strict         |
|  25  |Note                                   | Reason for nudge or other notes                |
|  26  |Hit_Start                              | Start co-ordinate for HSP in CARD reference    |
|  27  |Hit_End                                | End co-ordinate for HSP in CARD reference      |
|  28  |Antibiotic                             | ARO Categorization                             |

When viewing wide tab-delimited files like this one, it might be useful to look
at them one column at a time, which can be accomplished with the `cut` command.
For instance, if we wanted to look at 