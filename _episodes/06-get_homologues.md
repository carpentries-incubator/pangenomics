---
title: "Clustering Protein Sequences"
teaching: 30
exercises: 10
questions:
- "Can I cluster my sequences automatically?"
objectives:
- "Cluster proteins from GenBank files with an automatic software"
keypoints:
- "Clustering protein sequences refers to the process of grouping similar sequences into distinct clusters or families."
- "GET_HOMOLOGUES is a software package for microbial pangenome analysis"
- "Three sequence clustering algorithms are supported by GET_HOMOLOGUES; BDBH, COGtriangles, and OrthoMCL"
---

## Diversity of algorithms for clustering

The key step in pangenomics is knowing what to compare between genomes, in other words, what are the gene/protein families in our group of genomes? The clustering of homologous sequences is a very 
complex problem in bioinformatics, and it can be tackled in different ways. This is why there are many clustering programs that focus on different things to join the individual sequences into families. 

1. [OrthoFinder](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1832-y) uses a combination of sequence similarity and phylogenetic tree-based approaches to infer orthology relationships. It is designed for comparative genomics and phylogenetic analysis, not for pangenomics.

2. [Roary](https://academic.oup.com/bioinformatics/article/31/22/3691/240757) clusters proteins based on pairwise protein similarity. It was designed specifically for pangenomics so it is one of the most used programs in the field and it has set the standard format for the files that describe a pangenome.

3. [GET_HOMOLOGUES](https://journals.asm.org/doi/10.1128/AEM.02411-13) offers various algorithms for clustering proteins, including Bidirectional Best-Hits, Markov clustering, and COGtriangles. It also provides additional functionalities, such as the identification of strain-specific genes and visualization of pangenome data. It provides a comparison of different clustering algorithms.

4. [PPanGGOLiN](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007732) uses the program MMseqs2 to cluster proteins based on sequence similarity. It allows users to define the similarity threshold for clustering, enabling customization according to the specific requirements of the analysis. Also provides features for visualizing and exploring pangenome data. 

It's important to acknowledge the specific requirements of your analysis, such as scalability, speed, and the desired output, and evaluate different 
tools to determine which program best suits your needs. 

## What is GET_HOMOLOGUES?

In this episode, we will use the [GET_HOMOLOGUES](https://journals.asm.org/doi/10.1128/AEM.02411-13) suite of tools for pangenome analysis. Its main task is clustering protein and nucleotide sequences in homologous (possibly orthologous) groups. This software identifies orthologous groups of intergenic regions, flanked by orthologous open reading frames (ORFs), conserved across related genomes. The definition of pan- and core-genomes by GET_HOLOGUES is done by calculation of overlapping sets of proteins.

GET_HOMOLOGUES supports three sequence-clustering methods; Bidirectional Best-Hit (BDBH), OrthoMCL (OMCL), and COGtriangles clustering algorithms.

|    	Method   		 |                           	Definition                         		 |
|:---------------------:    |:---------------------------------------------------------------------:    |
| **Bidirectional Best-Hit (BDBH)**     |      	Clusters proteins by identifying reciprocal best hits between genomes.    		 |
|	**OrthoMCL (OMCL)**  	 | Uses graph theory to cluster proteins based on sequence similarity, handling paralogous genes and gene duplications.    |
|	**COGtriangles:**  	 |   	Assigns proteins to predefined functional categories (COGs) based on best matches to the COG database using a triangle inequality-based algorithm.  		 |


## Clustering protein families with GET_HOMOLOGUES

For this lesson, we will cluster all of our genomes with one of the algorithms of GET_HOMOLOGUES, but then, we will use 
PPanGGOLiN to analyze the pangenome, instead of the GET_HOMOLOGUES tools designed for pangenomics.
Before starting to use GET_HOMOLOGUES, we need to activate the `Pangenomics_Global` environment.

~~~
$ conda deactivate
$ conda activate /miniconda3/envs/Pangenomics_Global  
~~~
{: .language-bash}

Now let's display the programs' help to confirm that it is correctly installed and know what are the options to run the clustering.
~~~
$ get_homologues.pl -h #This command display the options
~~~
{: .language-bash}

~~~
-v print version, credits and checks installation
-d directory with input FASTA files ( .faa / .fna ),  		 (overrides -i,
   GenBank files ( .gbk ), 1 per genome, or a subdirectory 	 use of pre-clustered sequences
   ( subdir.clusters / subdir_ ) with pre-clustered sequences   ignores -c, -g)
   ( .faa / .fna ); allows new files to be added later;    
   creates output folder named 'directory_homologues'
-i input amino acid FASTA file with [taxon names] in headers,  (required unless -d is set)
   creates output folder named 'file_homologues'

Optional parameters:
-o only run BLAST/Pfam searches and exit              		 (useful to pre-compute searches)
-c report genome composition analysis                 		 (follows order in -I file if enforced,
                                                       		 ignores -r,-t,-e)
-R set random seed for genome composition analysis    		 (optional, requires -c, example -R 1234,
                                                       		 required for mixing -c with -c -a runs)
-m runmode [local|cluster|dryrun]                     		 (default local)
-n nb of threads for BLAST/HMMER/MCL in 'local' runmode   	 (default=2)
-I file with .faa/.gbk files in -d to be included     		 (takes all by default, requires -d)

Algorithms instead of default bidirectional best-hits (BDBH):
-G use COGtriangle algorithm (COGS, PubMed=20439257)  		 (requires 3+ genomes|taxa)
-M use orthoMCL algorithm (OMCL, PubMed=12952885)

[...]

Options that control clustering:
-t report sequence clusters including at least t taxa 		 (default t=numberOfTaxa,
                                                       		 t=0 reports all clusters [OMCL|COGS])
-a report clusters of sequence features in GenBank files  	 (requires -d and .gbk files,
   instead of default 'CDS' GenBank features           		 example -a 'tRNA,rRNA',
                                                       		 NOTE: uses blastn instead of blastp,
                                                       		 ignores -g,-D)
-g report clusters of intergenic sequences flanked by ORFs     (requires -d and .gbk files)
   in addition to default 'CDS' clusters
-f filter by %length difference within clusters       		 (range [1-100], by default sequence
                                                       		 length is not checked)
-r reference proteome .faa/.gbk file                  		 (by default takes file with
                                                       		 least sequences; with BDBH sets
                                                       		 first taxa to start adding genes)
-e exclude clusters with inparalogues                 		 (by default inparalogues are
                                                       		 included)
[...]
~~~
{: .output}

> ## Best practices
> GET_HOMOLOGUES suggests that the user runs their data inside a directory because in the future they might want to add a new *.gbk* file to the analysis.
{: .callout}

### Create a directory for the output

> ## Note for longer workshops
> Instead of running the following commands on all files, in longer workshops, it should be executed on the mini genomes.
> First, it will be useful to locate the mini genomes:
> ~~~
> $ ls ~/pan_workshop/data/annotated_mini/
> ~~~
> {: .language-bash}
> ~~~
> Streptococcus_agalactiae_2603V_mini.faa  Streptococcus_agalactiae_A909_mini.faa
> Streptococcus_agalactiae_515_mini.faa    Streptococcus_agalactiae_NEM316_mini.faa
> ~~~
> {: .output}

{: .callout}

It's necessary to create a new folder to store all the results.

~~~
$ mkdir -p ~/pan_workshop/results/pangenome/get_homologues/data_gbks #Create directory (option -p create all parent directories)
$ cd  ~/pan_workshop/results/pangenome/get_homologues/data_gbks #Locates you in the directory 'data_gbks'
~~~
{: .language-bash}
We need to create symbolic links to have easy access to all the *.gbk* files of our genomes.
~~~

$ ln -s ~/pan_workshop/results/annotated/Streptococcus_*/Streptococcus_agalactiae_*_prokka.gbk .
$ ls #List the symbolic links
~~~
{: .language-bash}

~~~
Streptococcus_agalactiae_18RS21_prokka.gbk  Streptococcus_agalactiae_CJB111_prokka.gbk
Streptococcus_agalactiae_2603V_prokka.gbk   Streptococcus_agalactiae_COH1_prokka.gbk
Streptococcus_agalactiae_515_prokka.gbk     Streptococcus_agalactiae_H36B_prokka.gbk
Streptococcus_agalactiae_A909_prokka.gbk    Streptococcus_agalactiae_NEM316_prokka.gbk
~~~
{: .output}

### Run the OMCL algorithm

To do the clustering we will use only the OMCL algorithm implemented in GET_HOMOLOGUES. 

Since the following command can take around 8 minutes to run we will use a screen session to run it. The screen session will not have the conda environment activated, so let’s activate it again.
~~~
$ cd ..
$ screen -R clustering
$ conda activate /miniconda3/envs/Pangenomics_Global
~~~
{: .language-bash}
And now let's run our program.
~~~
$ get_homologues.pl -d data_gbks -M -t 0 -c -n 1
~~~
{: .language-bash}
Click `Ctrl`+ `a` + `d` to detach from the session and wait 8 minutes to attach back the screen and check if it has finished.

~~~
# number_of_clusters = 3464
# cluster_list = data_gbks_homologues/Streptococcusagalactiae18RS21prokka_f0_0taxa_algOMCL_e0_.cluster_list
# cluster_directory = data_gbks_homologues/Streptococcusagalactiae18RS21prokka_f0_0taxa_algOMCL_e0_

# runtime: 479 wallclock secs (62.16 usr  1.81 sys + 1347.29 cusr 43.47 csys = 1454.73 CPU)
# RAM use: 89.0 MB
~~~
{: .output}

> ## Making single-copy clusters
> If the option `-e` is added, the resulting clusters will contain only single-copy genes from each taxon, i.e. orthologues. This flag excludes
> clusters within paralogues. This is useful for making genome-level phylogenetic analyses in only single-copy genes.
> However, it's important to note that using only orthologues may overlook genes that have undergone significant divergence or genes
> with species-specific functions. Including paralogues, which are duplicated genes within a genome, or genes with no clear orthologues
> can also be informative in certain pangenomic analyses, such as studying gene family expansions or specific adaptations within a species.
{: .callout}

## Describe your gene families in one table

GET_HOMOLOGUES gave us one FASTA file for each gene family, with the sequences of the genes included in the family. These files look like this:
~~~
$ head data_gbks_homologues/Streptococcusagalactiae18RS21prokka_f0_0taxa_algOMCL_e0_/3491_IS30_family_transpos...faa
~~~
{: .language-bash}
~~~
>ID:GBPINHCM_01628 |[Streptococcus agalactiae]|2603V|Streptococcus_agalactiae_2603V_prokka.gbk|IS30 family transpos..|804|NC_004116.1(2160267):1582200-1583003:1 ^COG:COG2826^ Streptococcus agalactiae strain 2603V.|neighbours:ID:GBPINHCM_01627(-1),ID:GBPINHCM_01629(-1)|neighbour_genes:tmk,IMPDH_2|
MGVKKGQRIYHILKTNDLEVSSSTVYRHIKKGYLSITPIDLPRAVKFKKRRKSTLPPIPKAIKEGRRYEDFIEHMNQSELNSWLEMDTVIGRIGGKVLLTFNVAFCNFIFAKLMDSKTAIETAKHIQVIKRTLYDNKRDFFELFPVILTDNGGEFARVDDIEIDVCGQSQLFFCDPNRSDQKARIEKNHTLVRDILPKGTSFDNLTQEDINLALSHINSVKRQALNGKTAYELFSFTYGKDIASILGIEEITAEDVCQSPKLLKDKI
~~~
{: .output}

We want to create a file that summarizes the information of the clustering by showing only which genes correspond to which family. 
We will need that file in the next episode to explore our pangenome with another program.

To obtain this file, which we will name `gene_families.tsv`, we will extract the IDs of the genes from the FASTA headers (in the FASTA header we see the ID of the gene after `ID:`) and the name of the families from the file names. For this, we will use the following short script.  

Copy the contents and paste them into a file:
~~~
$ nano ~/pan_workshop/scripts/get_gene_families.sh
~~~
{: .language-bash}

~~~
# Location to use: ~/pan_workshop/results/pangenome/get_homologues/
# Usage: bash ~/pan_workshop/scripts/get_gene_families.sh
# Output: 	One text file per gene family with the IDs of the genes it contains in the directory: 
# 				~/pan_workshop/results/pangenome/get_homologues/families/
# 			A tsv file with the name of the gene family in the first column and the name of the gene in the second column in:
# 				~/pan_workshop/results/pangenome/get_homologues/gene_families.tsv

mkdir families

# Obtain the gene IDs from the FASTA headers of the gene families FASTAs and put them in a text file:
for i in data_gbks_homologues/Streptococcusagalactiae18RS21prokka_f0_0taxa_algOMCL_e0_/*.faa
do 
base=$(basename $i .faa)
grep ">" $i | cut -d':' -f2 | cut -d' ' -f1 > families/${base}.txt
done

# Print the name of the family, a tab separator and the name of the gene in a tsv file:
for i in families/*.txt
do
base=$(basename $i .txt)
cat $i | while read line 
do
echo $base$'\t'$line >> gene_families.tsv
done
done
~~~
{: .output}

Now let's run this script to obtain the file.
~~~
$ bash ~/pan_workshop/scripts/get_gene_families.sh
~~~
{: .language.bash}

And let's see the file that was created.
~~~
$ head gene_families.tsv
~~~
{: .language-bash}

~~~
10003_Int-Tn_1	IGCLFMIO_00101
10003_Int-Tn_1	MGPKLEAL_00603
10004_hypothetical_protein	IGCLFMIO_00102
10004_hypothetical_protein	MGPKLEAL_00604
10005_hypothetical_protein	IGCLFMIO_00103
10005_hypothetical_protein	MGPKLEAL_00605
10006_hypothetical_protein	IGCLFMIO_00104
10006_hypothetical_protein	MGPKLEAL_00606
1000_zinT	GMBKAPON_01606
10010_upp_2	IGCLFMIO_01195
~~~
{: .output}

Now we have in only one file the description of our clustering results!

## Obtaining a pangenomic matrix for a shell genome database

We are going to obtain a file with our gene families that are found in the core genome, or well, perhaps a bit beyond the core, the Shell genome. This will be our database and will serve us to search for enzymes with expansions close to specialized metabolism with EvoMining.

Well, first, let's create our pangenomic matrix. This matrix has the presence and absence of genes from each gene family in the genomes.

~~~
$ compare_clusters.pl -o sample_intersection -m -d data_gbks_homologues/Streptococcusagalactiae18RS21prokka_f0_0taxa_algOMCL_e0_
~~~
{: .language-bash}

~~~
[...]
# intersection size = 3464 clusters

# intersection list = sample_intersection/intersection_t0.cluster_list

# pangenome_file = sample_intersection/pangenome_matrix_t0.tab transposed = sample_intersection/pangenome_matrix_t0.tr.tab
# pangenome_genes = sample_intersection/pangenome_matrix_genes_t0.tab transposed = sample_intersection/pangenome_matrix_genes_t0.tr.tab
# pangenome_phylip file = sample_intersection/pangenome_matrix_t0.phylip
# pangenome_FASTA file = sample_intersection/pangenome_matrix_t0.fasta
# pangenome CSV file (Scoary) = sample_intersection/pangenome_matrix_t0.tr.csv

# WARNING: Venn diagrams are only available for 2 or 3 input cluster directories
~~~
{: .output}

|	gene families / genomes	|                	Streptococcus_agalactiae_18RS21_prokka.gbk                         		 |	Streptococcus_agalactiae_2603V_prokka.gbk	|... |
|:-----------------:    |:---------------------------------------------------------------------:    | :------------------------------------------------:	| :---: |
| **1_rnmV_1.faa**  |						1 			   		 |			0				| ... |
| **2_scaC.faa**    |						 1				    	|			0				| ... |
| **3_hypothetical_protein.faa**  	 |   			1  					 |			0				| ... |
| **4_yecS.faa**  	 |   			1  					 |			1				| ... |
| 	...		|					...					|			...				|    ...    | ... |

Now we are going to use the table to count the gene families that are present in more than half of the genomes (>50%). The Shell genome consists of the genes shared by the majority of genomes (10-95% occurrence).


To achieve this, we will use a function in R that is available in the following [Github]([inst/extdata/search_shell_enzymes_DB_panworkshop.R](https://github.com/andrespan/MetaEvoMining/blob/14e1b45af93f53c6e2a716700ad9c7442d87a3e2/inst/extdata/search_shell_enzymes_DB_panworkshop.R)) repository. This function selects the protein sequences present in more than half of the genomes in the pangenomic matrix. It utilizes the sequence files obtained from Get_homologues and concatenates them into a FASTA file with the appropriate headers for use in EvoMining. Remember to import the necessary libraries which are in the '@import' section of the function documentation.

To facilitate the use of this lesson, let's download the *Streptococcus agalactiae* shell database which is available at this [GitHub](https://github.com/andrespan/MetaEvoMining/tree/e2945cd4555622c9677ca23d9ed1e161fe2c44fa/inst/extdata/DATABASES) link.

> ## References:
> Go to the GET_HOMOLOGUES [GitHub](https://github.com/eead-csic-compbio/get_homologues) for the complete collection of instructions and posibilities.  
> And read the original GET_HOMOLOGUES article to understand the details:
>
> Contreras-Moreira, B., Vinuesa, P. (2013) GET_HOMOLOGUES, a Versatile Software Package for Scalable and Robust Microbial Pangenome Analysis. Applied and Environmental Microbiology 79(24) [https://doi.org/10.1128/AEM.02411-13](https://journals.asm.org/doi/10.1128/aem.02411-13).
{: .callout}
