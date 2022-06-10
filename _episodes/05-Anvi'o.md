---
title: "Interactive pangenome graphics"
teaching: 15
exercises: 40
questions:
- "How can I obtain an interactive pangenome plot?"
- "How can I measure the homogeneity of the gene families?"
- "How to obtain an enrichement analysis of the gene families?"
- "How to compute the ANI values between the genomes of the pangenome?"

objectives:
- "Construct a pangenome following the Anvi'o workflow"
- "Visualize and interact with the pangenome graph"
- "Compute and visualize the ANI values of the genomes from the pangenome"
- "Perform a functional enrichment analysis on a group of genomes from the pangenome"

keypoints:
- "Anvi’o can build a pangenome starting from genomes or metagenomes, or a combination of both"
- "Anvi'o allows to interactively visualize your pangenomes"
- "Anvi'o platform includes additional scripts to explore the geometric and biochemical homogeneity of the gene clusters, to compute and visualize the ANI values of the genomes, to conduct a functional enrichment analysis in a group of genomes, among others"
---

<a href="../fig/01-03-01.png">
  <img src="../fig/01-03-01.png" alt="Anvi'o network representation" />
</a>

## Anvi'o

Anvi’o is an open-source, community-driven analysis and visualization platform for microbial ‘omics.
It brings together many aspects of today's cutting-edge strategies including **genomics, metagenomics, metatranscriptomics, phylogenomics, microbial population genetis, pangenomics and metapangenomics** in an *integrated* and *easy-to-use* fashion thorugh extensive interactive visualization capabilities. 


### Get all ready to start the Anvi'o workflow to build a pangenome

To start using Anvi'o, activate the conda environment `Pangenomics_Global` 

~~~
conda activate Pangenomics_Global
~~~
{: .laguage-bash}

~~~
(Pangenomics) betterlab@betterlabub:~$

~~~
{: .output}

Move into the directory named `results` and create a new directory called `anvi-o` for the Anvi'o analysis
~~~
cd ~/gm_workshop/results/pangenome
mkdir anvi-o

cd anvi-o
~~~
{: .laguage-bash}

In order to better organize our Anvi'o results, create a new directory named `genome-db` that will be used to storage the genome database needed for the Anvi'o pangenome worflow
~~~
mkdir genome-db
~~~
{: .laguage-bash}

**Note:** The bacterial genomes that will be used in this practice come from the Prokka annotation analysis. We will use the .gbk files as input for the Anvi'o workflow. The .gbk files can be found in ~/gm_workshop/results/annotated


Let's do it!


## Ten steps guide to build a Pangenome in Anvi'o

### Step 1

Process the genome files (.gbk) with the `anvi-script-process-genbank` script

~~~
ls ~/gm_workshop/results/annotated/agalactiae* | cut -d'/' -f7 | cut -d '.' -f1 | while read line; do anvi-script-process-genbank -i GENBANK --input-genbank ~/gm_workshop/results/annotated/$line.gbk -O genome-db/$line; done
~~~
{: .laguage-bash}

~~~
cd genome-db
ls
~~~
{: .laguage-bash}

~~~
agalactiae_18RS21_prokka-contigs.fa               agalactiae_CJB111_prokka-contigs.fa
agalactiae_18RS21_prokka-external-functions.txt   agalactiae_CJB111_prokka-external-functions.txt
agalactiae_18RS21_prokka-external-gene-calls.txt  agalactiae_CJB111_prokka-external-gene-calls.txt
agalactiae_515_prokka-contigs.fa                  agalactiae_COH1_prokka-contigs.fa
agalactiae_515_prokka-external-functions.txt      agalactiae_COH1_prokka-external-functions.txt
agalactiae_515_prokka-external-gene-calls.txt     agalactiae_COH1_prokka-external-gene-calls.txt
agalactiae_A909_prokka-contigs.fa                 agalactiae_H36B_prokka-contigs.fa
agalactiae_A909_prokka-external-functions.txt     agalactiae_H36B_prokka-external-functions.txt
agalactiae_A909_prokka-external-gene-calls.txt    agalactiae_H36B_prokka-external-gene-calls.txt

~~~
{: .output}

### Step 2

Reformat the fasta files using the `anvi-script-reformat-fasta` script

~~~
ls *fa |while read line; do anvi-script-reformat-fasta --seq-type NT $line -o $line\.fasta; done
ls
~~~
{: .laguage-bash}

~~~
agalactiae_18RS21_prokka-contigs.fa               agalactiae_CJB111_prokka-contigs.fa
agalactiae_18RS21_prokka-contigs.fa.fasta         agalactiae_CJB111_prokka-contigs.fa.fasta
agalactiae_18RS21_prokka-external-functions.txt   agalactiae_CJB111_prokka-external-functions.txt
agalactiae_18RS21_prokka-external-gene-calls.txt  agalactiae_CJB111_prokka-external-gene-calls.txt
agalactiae_515_prokka-contigs.fa                  agalactiae_COH1_prokka-contigs.fa
agalactiae_515_prokka-contigs.fa.fasta            agalactiae_COH1_prokka-contigs.fa.fasta
agalactiae_515_prokka-external-functions.txt      agalactiae_COH1_prokka-external-functions.txt
agalactiae_515_prokka-external-gene-calls.txt     agalactiae_COH1_prokka-external-gene-calls.txt
agalactiae_A909_prokka-contigs.fa                 agalactiae_H36B_prokka-contigs.fa
agalactiae_A909_prokka-contigs.fa.fasta           agalactiae_H36B_prokka-contigs.fa.fasta
agalactiae_A909_prokka-external-functions.txt     agalactiae_H36B_prokka-external-functions.txt
agalactiae_A909_prokka-external-gene-calls.txt    agalactiae_H36B_prokka-external-gene-calls.txt

~~~
{: .output}

### Step 3

Create a database per genome with the `anvi-gen-contigs-database` script

~~~
ls *fasta | while read line; do anvi-gen-contigs-database -T 4 -f $line -o $line-contigs.db; done
ls
~~~
{: .laguage-bash}

~~~
agalactiae_18RS21_prokka-contigs.fa                   agalactiae_CJB111_prokka-contigs.fa
agalactiae_18RS21_prokka-contigs.fa.fasta             agalactiae_CJB111_prokka-contigs.fa.fasta
agalactiae_18RS21_prokka-contigs.fa.fasta-contigs.db  agalactiae_CJB111_prokka-contigs.fa.fasta-contigs.db
agalactiae_18RS21_prokka-external-functions.txt       agalactiae_CJB111_prokka-external-functions.txt
agalactiae_18RS21_prokka-external-gene-calls.txt      agalactiae_CJB111_prokka-external-gene-calls.txt
agalactiae_515_prokka-contigs.fa                      agalactiae_COH1_prokka-contigs.fa
agalactiae_515_prokka-contigs.fa.fasta                agalactiae_COH1_prokka-contigs.fa.fasta
agalactiae_515_prokka-contigs.fa.fasta-contigs.db     agalactiae_COH1_prokka-contigs.fa.fasta-contigs.db
agalactiae_515_prokka-external-functions.txt          agalactiae_COH1_prokka-external-functions.txt
agalactiae_515_prokka-external-gene-calls.txt         agalactiae_COH1_prokka-external-gene-calls.txt
agalactiae_A909_prokka-contigs.fa                     agalactiae_H36B_prokka-contigs.fa
agalactiae_A909_prokka-contigs.fa.fasta               agalactiae_H36B_prokka-contigs.fa.fasta
agalactiae_A909_prokka-contigs.fa.fasta-contigs.db    agalactiae_H36B_prokka-contigs.fa.fasta-contigs.db
agalactiae_A909_prokka-external-functions.txt         agalactiae_H36B_prokka-external-functions.txt
agalactiae_A909_prokka-external-gene-calls.txt        agalactiae_H36B_prokka-external-gene-calls.txt

~~~
{: .output}

### Step 4

When using external genomes in anvi'o, a list of the genome ids and their corresponding genome database is required. This list tells Anvi'o which genomes will be processed to construct the pangenome. 
~~~
ls *.fa | cut -d '-' -f1 | while read line; do echo $line$'\t'$line-contigs.db >>external-genomes.txt; done
head external-genomes.txt
~~~
{: .laguage-bash}

~~~
agalactiae_18RS21_prokka        agalactiae_18RS21_prokka-contigs.db
agalactiae_515_prokka   agalactiae_515_prokka-contigs.db
agalactiae_A909_prokka  agalactiae_A909_prokka-contigs.db
agalactiae_CJB111_prokka        agalactiae_CJB111_prokka-contigs.db
agalactiae_COH1_prokka  agalactiae_COH1_prokka-contigs.db
agalactiae_H36B_prokka  agalactiae_H36B_prokka-contigs.db
~~~
{: .output}

### Step 5

Modify the headers of the list `external-genomes.txt`
~~~
nano external-genomes.txt
~~~
{: .laguage-bash}

~~~
  GNU nano 4.8                                                             external-genomes.txt                                                                       
agalactiae_18RS21_prokka        agalactiae_18RS21_prokka-contigs.db
agalactiae_515_prokka   agalactiae_515_prokka-contigs.db
agalactiae_A909_prokka  agalactiae_A909_prokka-contigs.db
agalactiae_CJB111_prokka        agalactiae_CJB111_prokka-contigs.db
agalactiae_COH1_prokka  agalactiae_COH1_prokka-contigs.db
agalactiae_H36B_prokka  agalactiae_H36B_prokka-contigs.db



^G Get Help     ^O Write Out    ^W Where Is     ^K Cut Text     ^J Justify      ^C Cur Pos      M-U Undo        M-A Mark Text   M-] To Bracket  M-Q Previous
^X Exit         ^R Read File    ^\ Replace      ^U Paste Text   ^T To Spell     ^_ Go To Line   M-E Redo        M-6 Copy Text   ^Q Where Was    M-W Next
~~~
{: .output}

~~~
head external-genomes.txt
~~~
{: .laguage-bash}

~~~
name    contigs_db_path
agalactiae_18RS21_prokka        agalactiae_18RS21_prokka-contigs.db
agalactiae_515_prokka   agalactiae_515_prokka-contigs.db
agalactiae_A909_prokka  agalactiae_A909_prokka-contigs.db
agalactiae_CJB111_prokka        agalactiae_CJB111_prokka-contigs.db
agalactiae_COH1_prokka  agalactiae_COH1_prokka-contigs.db
agalactiae_H36B_prokka  agalactiae_H36B_prokka-contigs.db

~~~
{: .output}

### Step 6

Rename the `.db` files

~~~
rename s'/.fa.fasta-contigs.db/.db/' *db
ls *.db
~~~
{: .laguage-bash}

~~~
agalactiae_18RS21_prokka-contigs.db  agalactiae_A909_prokka-contigs.db    agalactiae_COH1_prokka-contigs.db
agalactiae_515_prokka-contigs.db     agalactiae_CJB111_prokka-contigs.db  agalactiae_H36B_prokka-contigs.db

~~~
{: .output}

### Step 7
Execute HMM analysis with the `anvi-run-hmms` script to identify matching genes in each contigs database file

~~~
ls *contigs.db | while read line; do anvi-run-hmms -c $line; done
~~~
{: .laguage-bash}

~~~

HMM Profiling for Ribosomal_RNA_16S
===============================================
Reference ....................................: Seeman T, https://github.com/tseemann/barrnap
Kind .........................................: Ribosomal_RNA_16S
Alphabet .....................................: RNA
Context ......................................: CONTIG
Domain .......................................: N/A
HMM model path ...............................: /tmp/tmpjguiut54/Ribosomal_RNA_16S.hmm
Number of genes in HMM model .................: 3
Noise cutoff term(s) .........................: --cut_ga
Number of CPUs will be used for search .......: 1
HMMer program used for search ................: nhmmscan
Temporary work dir ...........................: /tmp/tmp9q6mevny
Log file for thread 0 ........................: /tmp/tmp9q6mevny/RNA_contig_sequences.fa.0_log
Done 🎊

Number of raw hits in table file .............: 7
Number of weak hits removed by HMMER parser ..: 0
Number of hits in annotation dict  ...........: 7
Pruned .......................................: 3 out of 7 hits were removed due to redundancy
Gene calls added to db .......................: 4 (from source "Ribosomal_RNA_16S")

HMM Profiling for Ribosomal_RNA_23S
===============================================
Reference ....................................: Seeman T, https://github.com/tseemann/barrnap
Kind .........................................: Ribosomal_RNA_23S
Alphabet .....................................: RNA
Context ......................................: CONTIG
Domain .......................................: N/A
HMM model path ...............................: /tmp/tmpjguiut54/Ribosomal_RNA_23S.hmm
Number of genes in HMM model .................: 2
Noise cutoff term(s) .........................: --cut_ga
Number of CPUs will be used for search .......: 1
HMMer program used for search ................: nhmmscan
Temporary work dir ...........................: /tmp/tmp9q6mevny
Log file for thread 0 ........................: /tmp/tmp9q6mevny/RNA_contig_sequences.fa.0_log
Done 🎊

Number of raw hits in table file .............: 7
Number of weak hits removed by HMMER parser ..: 0
Number of hits in annotation dict  ...........: 7
Pruned .......................................: 3 out of 7 hits were removed due to redundancy
Gene calls added to db .......................: 4 (from source "Ribosomal_RNA_23S")

~~~
{: .output}


### Step 8

Create the genome database `genomes-storage-db` using the `anvi-gen-genomes-storage` script. In this case, we named this `genomes-storage-db` as **AGALACTIAE_GENOMES.db**, which will be used downstream as input in other process.

~~~
anvi-gen-genomes-storage -e external-genomes.txt -o AGALACTIAE_GENOMES.db
ls *.db
~~~
{: .laguage-bash}

~~~
agalactiae_18RS21_prokka-contigs.db  agalactiae_CJB111_prokka-contigs.db  agalactiae_H36B_prokka-contigs.db
agalactiae_515_prokka-contigs.db     agalactiae_COH1_prokka-contigs.db
agalactiae_A909_prokka-contigs.db    AGALACTIAE_GENOMES.db
~~~
{: .output}

### Step 9

Construct the pangenome database `pan-db` with the `anvi-pan-pangenome` script using the `genomes-storage-db` named **AGALACTIAE_GENOMES.db** as input 

~~~
anvi-pan-genome -g AGALACTIAE_GENOMES.db \
                --project-name "PANGENOME-AGALACTIAE" \
                --output-dir AGALACTIAE \
                --num-threads 6 \
                --minbit 0.5 \
                --mcl-inflation 10 \
                --use-ncbi-blast
~~~
{: .laguage-bash}

~~~
WARNING
===============================================
If you publish results from this workflow, please do not forget to cite DIAMOND
(doi:10.1038/nmeth.3176), unless you use it with --use-ncbi-blast flag, and MCL
(http://micans.org/mcl/ and doi:10.1007/978-1-61779-361-5_15)

Functions found ..............................:
Genomes storage ..............................: Initialized (storage hash: hash299bb5bf)
Num genomes in storage .......................: 6
Num genomes will be used .....................: 6
Pan database .................................: A new database,
                                                /home/betterlab/gm_workshop/results/anvi-o/genome-db/AGALACTIAE/PANGENOME-AGALACTIAE-PAN.db,
                                                has been created.
Exclude partial gene calls ...................: False

AA sequences FASTA ...........................: /home/betterlab/gm_workshop/results/anvi-o/genome-db/AGALACTIAE/combined-aas.fa

Num AA sequences reported ....................: 13,548
Num excluded gene calls ......................: 0
Unique AA sequences FASTA ....................: /home/betterlab/gm_workshop/results/anvi-o/genome-db/AGALACTIAE/combined-aas.fa.unique

WARNING
===============================================
You elected to use NCBI's `blastp` for amino acid sequence search. Running
blastp will be significantly slower than DIAMOND, but in some cases, slightly
more sensitive. We are unsure about whether the slight increase in sensitivity
may justify significant increase in run time, but you are the boss.


NCBI BLAST MAKEDB
===============================================
BLAST search db ..............................: /home/betterlab/gm_workshop/results/pangenome/anvi-o/genome-db/AGALACTIAE/combined-aas.fa.unique

NCBI BLAST SEARCH
===============================================
BLAST results ................................: /home/betterlab/gm_workshop/results/pangenome/anvi-o/genome-db/AGALACTIAE/blast-search-results.txt

MCL INPUT
===============================================
Min percent identity .........................: 0.0
Minbit .......................................: 0.5
Filtered search results ......................: 92,757 edges stored
MCL input ....................................: /home/betterlab/gm_workshop/results/pangenome/anvi-o/genome-db/AGALACTIAE/mcl-input.txt

MCL
===============================================
MCL inflation ................................: 10.0
MCL output ...................................: /home/betterlab/gm_workshop/results/pangenome/anvi-o/genome-db/AGALACTIAE/mcl-clusters.txt
Number of MCL clusters .......................: 2,711

CITATION
===============================================
The workflow you are using will likely use 'muscle' by Edgar,
doi:10.1093/nar/gkh340 (http://www.drive5.com/muscle) to align your sequences.
If you publish your findings, please do not forget to properly credit this tool.

* Your pangenome is ready with a total of 2,711 gene clusters across 6 genomes 🎉

~~~
{: .output}


### Step 10

Create the interactive pangenome with the `anvi-display-pan` script using as input the `genomes-storage-db` **AGALACTIAE_GENOMES.db** and the `pan-db` **PANGENOME-AGALACTIAE-PAN.db** (located in AGALACTIAE directory)

~~~
anvi-display-pan -g AGALACTIAE_GENOMES.db \
    -p AGALACTIAE/PANGENOME-AGALACTIAE-PAN.db
~~~
{: .laguage-bash}

~~~
* The server is up and running 🎉

WARNING
===============================================
If you are using OSX and if the server terminates prematurely before you can see
anything in your browser, try running the same command by putting 'sudo ' at the
beginning of it (you will be prompted to enter your password if sudo requires
super user credentials on your system). If your browser does not show up, try
manually entering the URL shown below into the address bar of your favorite
browser. *cough* CHROME *cough*.


Server address ...............................: http://0.0.0.0:8080

* When you are ready, press CTRL+C once to terminate the server and go back to the
command line.


~~~
{: .output}


Whitout disturbing the active terminal, open a new window in your prefered browser (recommended Chrome), copy-paste the following link `http://132.248.196.38:8080` and click on the bottom `Draw` to see your results and start interacting with your pangenome

<a href="../fig/01-03-02.svg">
  <img src="../fig/01-03-02.svg" width="956.5" height="453.5" alt="" />
</a>

{: .output}

  
> ## Exercise 1: The homogeneity of gene clusters.
>Anvi’oo allows us to identify different levels of disagreements between amino acid sequences in different genomes. Amino acid sequences from different genomes in a gene cluster that are almost identical tell us that the gene cluster is highly homogeneous. 
>
> The **geometric homogeneity index** tell us the degree of geometric configuration between the genes of a gene cluster and the **functional homogeneity index** considers aligned residues and quantifies differences across residues in a site.
>
>For more info see [this.](https://merenlab.org/2016/11/08/pangenomics-v2/#inferring-the-homogeneity-of-gene-clusters)
>
>Go to this [page](https://anvio.org/help/main/programs/anvi-get-sequences-for-gene-clusters/) and explore the pangenome graph according to the following homogeneity index.
>
>a) Order the pangenome based on the geometric homogeneity index and inspect a gene cluster with a relatively low score.
>
>b) Filter the gene cluster according to a functional homogeneity index above 0.25. 
>
>Extra) How can you estimate evolutionary relationships between genomes? With the `concatenated-gene-alignment-fasta` produce the phylogenomic tree and explore it.
> >## Solution
>> a) Go to the main settings panel and modify the “items order”.
>>
>> b) 
>>~~~
anvi-get-sequences-for-gene-clusters -g genomes-storage-db \
                                     -p pan-db \
                                     -o genes-fasta \
                                     --min-functional-homogenity-index 0.25
>>~~~
>>{: .laguage-bash}
>>
>>Extra)Explore this [page.](https://anvio.org/help/main/programs/anvi-gen-phylogenomic-tree/)
>{: .solution}
{: .challenge}
  
  
> ## Exercise 2: Splitting the pangenome.
> 1. Read about [anvi-split](https://anvio.org/help/main/programs/anvi-split/) 
> 2. With this program split your pangenome in independent pangenomes that:
> > * Contains only singletons.
> > * Contains only core gene clusters.
>
> Tip: [anvi-display-pan](https://anvio.org/help/main/programs/anvi-display-pan/) can be usefull
{: .challenge}

{% include links.md %}
