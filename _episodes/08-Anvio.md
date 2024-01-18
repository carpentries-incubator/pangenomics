---
title: "Interactive Pangenome Plots"
teaching: 15
exercises: 40
questions:
- "How can I obtain an interactive pangenome plot?"
- "How can I measure the homogeneity of the gene families?"
- "How to obtain an enrichment analysis of the gene families?"
- "How to compute the ANI values between the genomes of the pangenome?"

objectives:
- "Construct a pangenome following the Anvi'o workflow"
- "Visualize and interact with the pangenome graph"
- "Compute and visualize the ANI values of the genomes from the pangenome"
- "Perform a functional enrichment analysis on a group of genomes from the pangenome"

keypoints:
- "Anvi’o can build a pangenome starting from genomes or metagenomes, or a combination of both"
- "Anvi'o allows you to interactively visualize your pangenomes"
- "Anvi'o platform includes additional scripts to explore the geometric and biochemical homogeneity of the gene clusters, to compute and visualize the ANI values of the genomes, to conduct a functional enrichment analysis in a group of genomes, among others"
---

## Anvi'o

Anvi’o is an open-source, community-driven analysis and visualization platform for microbial omics.
It brings together many aspects of today's cutting-edge strategies, including **genomics, metagenomics, metatranscriptomics, phylogenomics, microbial population genetics, pangenomics, and metapangenomics** in an *integrated* and *easy-to-use* fashion through extensive interactive visualization capabilities.

In this episode, we will meet another pangenomics powerful tool. The pangenomics workflow of Anvi'o is not suitable for thousands 
of genomes like the PPanGGOLiN workflow, but it allows for an interactive exploration of smaller pangenomes, providing you with a 
closer look at the pangenome matrix and some interesting characteristics of our gene families.

## Preparing the databases for each genome

To start using Anvi'o, activate the conda environment `Pangenomics_Global`.
~~~
$ conda activate /miniconda3/envs/anvio-7.1
~~~
{: .language-bash}

Move into the directory named `results` and create a new directory called `anvi-o` for the Anvi'o analysis.
~~~
$ cd ~/pan_workshop/results/pangenome
$ mkdir anvi-o
$ cd anvi-o
~~~
{: .language-bash}

In order to better organize our Anvi'o results, create a new directory named `genome-db` that will be used to store the genome database needed for the Anvi'o pangenome workflow. We will use the `.gbk` files that came out of Prokka as input for the Anvi'o workflow. They can be found in `~/pan_workshop/results/annotated`.
~~~
$ mkdir genome-db
~~~
{: .language-bash}

To build a pangenome Anvi'o needs to extract the sequences AND OTHER INFORMATION from the `gbk` files. To do this
let's do a while loop to get the file names of each `gbk` and run the `anvi-script-process-genbank` script for each of them.
~~~
$ ls ~/pan_workshop/results/annotated/Streptococcus_agalactiae_* | cut -d'/' -f7 | cut -d '.' -f1 | while read line
do
anvi-script-process-genbank -i GENBANK --input-genbank ~/pan_workshop/results/annotated/$line.gbk -O genome-db/$line
done
~~~
{: .language-bash}

~~~
$ cd genome-db
$ ls
~~~
{: .language-bash}

~~~
Streptococcus_agalactiae_18RS21_prokka-contigs.fa 
Streptococcus_agalactiae_18RS21_prokka-external-functions.txt 
Streptococcus_agalactiae_18RS21_prokka-external-gene-calls.txt 
...
~~~
{: .output}

Now we have our `genome-db/` with the files that Anvi'o needs. Now we need to reformat the generated `fasta` files so that the GENE NAMES ARE STANDARDIZED. Let's do it
with the `anvi-script-reformat-fasta` script.

~~~
$ ls *fa |while read line
do
anvi-script-reformat-fasta --seq-type NT $line -o $line\.fasta
done
$ ls
~~~
{: .language-bash}

~~~
Streptococcus_agalactiae_18RS21_prokka-contigs.fa 
Streptococcus_agalactiae_18RS21_prokka-contigs.fa.fasta 
Streptococcus_agalactiae_18RS21_prokka-external-functions.txt 
Streptococcus_agalactiae_18RS21_prokka-external-gene-calls.txt 
...
~~~
{: .output}

With these new files now we need to GATHER THAT INFORMATION in an Anvi'o database format. To create a database per genome we need to run the `anvi-gen-contigs-database` script.

~~~
$ ls *fasta | while read line; do anvi-gen-contigs-database -T 4 -f $line -o $line-contigs.db; done
$ ls
~~~
{: .language-bash}

~~~
Streptococcus_agalactiae_18RS21_prokka-contigs.fa 
Streptococcus_agalactiae_18RS21_prokka-contigs.fa.fasta 
Streptococcus_agalactiae_18RS21_prokka-contigs.fa.fasta-contigs.db 
Streptococcus_agalactiae_18RS21_prokka-external-functions.txt 
Streptococcus_agalactiae_18RS21_prokka-external-gene-calls.txt 
...
~~~
{: .output}

The database files have a super long name, so we should replace the extension to only `.db`.

~~~
$ rename s'/.fa.fasta-contigs.db/.db/' *db
$ ls *.db
~~~
{: .language-bash}

~~~
Streptococcus_agalactiae_18RS21_prokka-contigs.db 
...
~~~
{: .output}

When using external genomes (genomes that are not part of the Anvi'o collection), a list of the genome IDs and their corresponding genome database is required. This list tells Anvi'o which genomes will be processed to construct the pangenome.
~~~
$ ls *.fa | cut -d '-' -f1 | while read line
do
echo $line$'\t'$line-contigs.db >> external-genomes.txt
done
$ head external-genomes.txt
~~~
{: .language-bash}

~~~
Streptococcus_agalactiae_18RS21_prokka  Streptococcus_agalactiae_18RS21_prokka-contigs.db
Streptococcus_agalactiae_2603V_prokka   Streptococcus_agalactiae_2603V_prokka-contigs.db
Streptococcus_agalactiae_515_prokka     Streptococcus_agalactiae_515_prokka-contigs.db
Streptococcus_agalactiae_A909_prokka    Streptococcus_agalactiae_A909_prokka-contigs.db
Streptococcus_agalactiae_CJB111_prokka  Streptococcus_agalactiae_CJB111_prokka-contigs.db
Streptococcus_agalactiae_COH1_prokka    Streptococcus_agalactiae_COH1_prokka-contigs.db
Streptococcus_agalactiae_H36B_prokka    Streptococcus_agalactiae_H36B_prokka-contigs.db
Streptococcus_agalactiae_NEM316_prokka  Streptococcus_agalactiae_NEM316_prokka-contigs.db
~~~
{: .output}

Let's add a header to the list that we made.
~~~
$ nano external-genomes.txt
~~~
{: .language-bash}

~~~
name    contigs_db_path
Streptococcus_agalactiae_18RS21_prokka  Streptococcus_agalactiae_18RS21_prokka-contigs.db
Streptococcus_agalactiae_2603V_prokka   Streptococcus_agalactiae_2603V_prokka-contigs.db
Streptococcus_agalactiae_515_prokka     Streptococcus_agalactiae_515_prokka-contigs.db
Streptococcus_agalactiae_A909_prokka    Streptococcus_agalactiae_A909_prokka-contigs.db
Streptococcus_agalactiae_CJB111_prokka  Streptococcus_agalactiae_CJB111_prokka-contigs.db
Streptococcus_agalactiae_COH1_prokka    Streptococcus_agalactiae_COH1_prokka-contigs.db
Streptococcus_agalactiae_H36B_prokka    Streptococcus_agalactiae_H36B_prokka-contigs.db
Streptococcus_agalactiae_NEM316_prokka  Streptococcus_agalactiae_NEM316_prokka-contigs.db
~~~
{: .output}

## Building the pangenome

### HMM
Now we are ready to identify matching genes in each contigs database file, for this, we will execute the HMM analysis with the `anvi-run-hmms` script.

~~~
$ ls *contigs.db | while read line
do
anvi-run-hmms -c $line
done
~~~
{: .language-bash}

~~~
Contigs DB ...................................: Streptococcus_agalactiae_18RS21_prokka-contigs.db
HMM sources ..................................: Ribosomal_RNA_5S, Ribosomal_RNA_12S, Bacteria_71, Ribosomal_RNA_16S, Archaea_76,
                                                Ribosomal_RNA_28S, Ribosomal_RNA_18S, Protista_83, Ribosomal_RNA_23S
Alphabet/context target found ................: AA:GENE
Alphabet/context target found ................: RNA:CONTIG

HMM Profiling for Ribosomal_RNA_5S
===============================================
Reference ....................................: Seeman T, https://github.com/tseemann/barrnap
Kind .........................................: Ribosomal_RNA_5S
Alphabet .....................................: RNA
Context ......................................: CONTIG
Domain .......................................: N/A
HMM model path ...............................: /tmp/tmpigzysqa6/Ribosomal_RNA_5S.hmm
Number of genes in HMM model .................: 5
Noise cutoff term(s) .........................: --cut_ga
Number of CPUs will be used for search .......: 1
HMMer program used for search ................: nhmmscan
Temporary work dir ...........................: /tmp/tmplk81rft0
Log file for thread 0 ........................: /tmp/tmplk81rft0/RNA_contig_sequences.fa.0_log
Done 🎊

Number of raw hits in table file .............: 0

* The HMM source 'Ribosomal_RNA_5S' returned 0 hits. SAD (but it's stil OK).

HMM Profiling for Ribosomal_RNA_12S
===============================================
Reference ....................................: Seeman T, https://github.com/tseemann/barrnap
Kind .........................................: Ribosomal_RNA_12S
Alphabet .....................................: RNA
Context ......................................: CONTIG
Domain .......................................: N/A
HMM model path ...............................: /tmp/tmpigzysqa6/Ribosomal_RNA_12S.hmm
Number of genes in HMM model .................: 1
Noise cutoff term(s) .........................: --cut_ga
Number of CPUs will be used for search .......: 1
HMMer program used for search ................: nhmmscan
Temporary work dir ...........................: /tmp/tmplk81rft0
Log file for thread 0 ........................: /tmp/tmplk81rft0/RNA_contig_sequences.fa.0_log
Done 🎊

Number of raw hits in table file .............: 0

* The HMM source 'Ribosomal_RNA_12S' returned 0 hits. SAD (but it's stil OK).

~~~
{: .output}

> ## Know more
> If you want to read more about HMM .
{: .callout}

### Creating a combined database

Create the genome database `genomes-storage-db` using the `anvi-gen-genomes-storage` script. In this case, we named this `genomes-storage-db` as **STREPTOCOCCUS_AGALACTIAE_GENOMES.db**, which will be used downstream as input in other processes.

~~~
$ anvi-gen-genomes-storage -e external-genomes.txt -o STREPTOCOCCUS_AGALACTIAE_GENOMES.db
$ ls *.db
~~~
{: .language-bash}

~~~
Streptococcus_agalactiae_18RS21_prokka-contigs.db  Streptococcus_agalactiae_COH1_prokka-contigs.db
Streptococcus_agalactiae_2603V_prokka-contigs.db   STREPTOCOCCUS_AGALACTIAE_GENOMES.db
Streptococcus_agalactiae_515_prokka-contigs.db     Streptococcus_agalactiae_H36B_prokka-contigs.db
Streptococcus_agalactiae_A909_prokka-contigs.db    Streptococcus_agalactiae_NEM316_prokka-contigs.db
Streptococcus_agalactiae_CJB111_prokka-contigs.db
~~~
{: .output}

### Making a pangenomic database

Construct the pangenome database `pan-db` with the `anvi-pan-pangenome` script using the `genomes-storage-db` named `STREPTOCOCCUS_AGALACTIAE_GENOMES.db` as input.

The desciption of this script is the next using the flag '-g' indicated that to create a pangenome using the database GENOMES.db, --project-name indicate the name of your preference, '--num-threads' indicate the number of cores tha using your computer for do the procees in this case 6. The flag '--minbit 0.5'. The flag '--mcl-inflation' indicate the sensitive to grouping a genes in the genomes. Whe use this flag '--use-ncbi-blast' is the program that use for align the genes

ABEL EXPLICA LAS FLAGS DE ESTO

~~~
$ anvi-pan-genome -g STREPTOCOCCUS_AGALACTIAE_GENOMES.db \
            	--project-name "PANGENOME-AGALACTIAE" \
            	--output-dir AGALACTIAE \
            	--num-threads 6 \
            	--minbit 0.5 \
            	--mcl-inflation 10 \
            	--use-ncbi-blast
~~~
{: .language-bash}

> ## Know more
> If you want to read more about th flags "minbit" and "MCL".
> Too-weak matches were culled by employing the --minbit criterion with the default value of 0.5, to minimize feeding MCL irrelevant similarities - by calculating all 
possible pairwise similarities, we introduce nonsensical comparisons, e.g., DnaK to TonB, that can be disregarded by their poor alignment quality prior to MCL. CITAR A DUTTER
> MCL uses a hyperparameter, “inflation,” to adjust the clustering sensitivity, i.e., the tendency to split clusters. MCL then uses these pairwise identities to group ORFs into gene clusters, putatively homologous gene groups. MCL uses a hyperparameter (inflation, --mcl-inflation) to adjust the clustering sensitivity, i.e., the tendency to split clusters. The decision for this flag depend of the level of you genomes, for example if you want to construct a pangenome to level Genus you can use a MCL more low and if you want to construct a pangenome level you can use more high. {: .callout}



~~~
WARNING
===============================================
If you publish results from this workflow, please do not forget to cite DIAMOND
(doi:10.1038/nmeth.3176), unless you use it with --use-ncbi-blast flag, and MCL
(http://micans.org/mcl/ and doi:10.1007/978-1-61779-361-5_15)

Functions found ..............................:
Genomes storage ..............................: Initialized (storage hash: hash8a837d50)
Num genomes in storage .......................: 8
Num genomes will be used .....................: 8
Pan database .................................: A new database,
                                                /home/shaday/pan_workshop/results/pangenome/anvi-o/genome-db/AGALACTIAE/PANGENOME-AGALACTIAE-PAN.db,
                                                has been created.
Exclude partial gene calls ...................: False

AA sequences FASTA ...........................: /home/shaday/pan_workshop/results/pangenome/anvi-o/genome-db/AGALACTIAE/combined-aas.fa

Num AA sequences reported ....................: 17,199
Num excluded gene calls ......................: 0
Unique AA sequences FASTA ....................: /home/shaday/pan_workshop/results/pangenome/anvi-o/genome-db/AGALACTIAE/combined-aas.fa.unique

WARNING
===============================================
You elected to use NCBI's `blastp` for amino acid sequence search. Running
blastp will be significantly slower than DIAMOND, but in some cases, slightly
more sensitive. We are unsure about whether the slight increase in sensitivity
may justify significant increase in run time, but you are the boss.


NCBI BLAST MAKEDB
===============================================
BLAST search db ..............................: /home/shaday/pan_workshop/results/pangenome/anvi-o/genome-db/AGALACTIAE/combined-aas.fa.unique

NCBI BLAST SEARCH
===============================================
BLAST results ................................: /home/shaday/pan_workshop/results/pangenome/anvi-o/genome-db/AGALACTIAE/blast-search-results.txt

MCL INPUT
===============================================
Min percent identity .........................: 0.0
Minbit .......................................: 0.5
Filtered search results ......................: 140,540 edges stored
MCL input ....................................: /home/shaday/pan_workshop/results/pangenome/anvi-o/genome-db/AGALACTIAE/mcl-input.txt

MCL
===============================================
MCL inflation ................................: 10.0
MCL output ...................................: /home/shaday/pan_workshop/results/pangenome/anvi-o/genome-db/AGALACTIAE/mcl-clusters.txt
Number of MCL clusters .......................: 2,842

CITATION
===============================================
The workflow you are using will likely use 'muscle' by Edgar,
doi:10.1093/nar/gkh340 (http://www.drive5.com/muscle) to align your sequences.
If you publish your findings, please do not forget to properly credit this tool.

* Your pangenome is ready with a total of 2,842 gene clusters across 8 genomes 🎉

~~~
{: .output}


## Creating an interactive plot

Create the interactive pangenome with the `anvi-display-pan` script using as input the `genomes-storage-db`  `STREPTOCOCCUS_AGALACTIAE_GENOMES.db` and the `pan-db`  `PANGENOME-AGALACTIAE-PAN.db` (located in `AGALACTIAE` directory)

~~~
$ anvi-display-pan -g STREPTOCOCCUS_AGALACTIAE_GENOMES.db \
	-p AGALACTIAE/PANGENOME-AGALACTIAE-PAN.db
~~~
{: .language-bash}

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


Without disturbing the active terminal, open a new window in your preferred browser (recommended Chrome), copy-paste the following link `http://bioinformatica.matmor.unam.mx:8080` and click on the bottom `Draw` to see your results and start interacting with your pangenome

<a href="../fig/01-03-02.svg">
  <img src="../fig/01-03-02.svg" width="956.5" height="453.5" alt="Interactive Anvio pan genome analysis of six S. agalactiae genomes.
                                                               	Each circle corresponds to one genome and each radius represents a gene family. " />
</a>

{: .output}
> ## Exercise 1: Explore the interactive plot.
> 
>>## Solution
>>
>>~~~
>>{: .language-bash}
>>
>>
>{: .solution}
{: .challenge}

## Special analyses
Anvi'o allows us to identify different levels of disagreement between amino acid sequences in different genomes. Amino acid sequences from different genomes in a gene cluster that are almost identical tell us that the gene cluster is highly homogeneous. The **geometric homogeneity index** tells us the degree of geometric configuration between the genes of a gene cluster and the **functional homogeneity index** considers aligned residues and quantifies differences across residues in a site. For more info see [this.](https://merenlab.org/2016/11/08/pangenomics-v2/#inferring-the-homogeneity-of-gene-clusters)

Anvi'o can also estimate evolutionary relationships between genomes with the `concatenated-gene-alignment-fasta` to produce the phylogenomic tree. For more information see [this](https://anvio.org/help/main/programs/anvi-gen-phylogenomic-tree/)

 With Anvi'o you can further analyze your pangenome with [`anvi-split`](https://anvio.org/help/main/programs/anvi-split/) to create independent pangenomes that contain only singletons or contain only core gene clusters.

{% include links.md %}
