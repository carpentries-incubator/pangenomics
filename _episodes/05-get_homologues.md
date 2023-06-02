---
title: "Clustering Protein Families"
teaching: 15
exercises: 10
questions:
- "What software is recommended for clustering protein families?"
- "What is GET_HOMOLOGUES?"
objectives:
- "Cluster orthologous proteins from GenBank (gbk) files"
- "Explore clusters using GET_HOMOLOGUES suit of tools"
- "Understand basic pangenomics metrics"
keypoints:
- "Clustering protein families refers to the process of grouping proteins that share similar characteristics or functions into distinct clusters or families."
- "GET_HOMOLOGUES is a software package for microbial pangenome analysis"
- "Three sequence clustering algorithms are supported by Get_Homologues; BDBH, COGtriangles and OrthoMCL"

---

## What software is recommended for clustering protein families?

When it comes to pangenome analysis, which involves analyzing the complete set of genes in a given species or a group of related organisms, there are several software tools that can be used for clustering protein families:

1. [OrthoFinder](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1832-y): It clusters proteins based on orthology relationships, identifying orthologous protein families across multiple genomes. OrthoFinder uses a combination of sequence similarity and phylogenetic tree-based approaches to infer orthology relationships.

2. [Roary](https://academic.oup.com/bioinformatics/article/31/22/3691/240757): It clusters proteins based on pairwise protein similarity. It utilizes a fast algorithm to construct clusters and determines core and accessory genes in the pangenome. Roary is known for its speed and scalability, making it suitable for large-scale pangenome analyses.

3. [GET_HOMOLOGUES](https://journals.asm.org/doi/10.1128/AEM.02411-13): It offers various algorithms for clustering proteins, including bidirectional best hit, Markov clustering, and COGtriangles. Also provides additional functionalities, such as identification of strain-specific genes and visualization of pangenome data.

4. [PPanGGOLiN](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007732): It uses the CD-HIT algorithm to cluster proteins based on sequence similarity. It allows users to define the similarity threshold for clustering, enabling customization according to the specific requirements of the analysis. Also provides features for visualizing and exploring pangenome data. 

> ## Considerations when choosing
> It's important to acknowledge the specific requirements of your analysis, such as scalability, speed, and the desired output, and evaluate different 
> tools to determine which one best suits your needs. 
{: .callout}

## What is GET_HOMOLOGUES?

In this episode we will use get [GET_HOMOLOGUES](https://journals.asm.org/doi/10.1128/AEM.02411-13) suite of tools for pangenome analysis.

Its main task is clustering protein and nucleotide sequences in homologous (possibly orthologous) groups. This software identifies orthologous groups of intergenic regions, flanked by orthologous open reading frames (ORFs), conserved across related genomes. The definition of pan- and core-genomes by Get_Homologues is done by calculation of overlapping sets of proteins. It is maintained by Bruno Contreras-Moreira and Pablo Vinuesa.

GET_HOMOLOGUES supports three sequence-clustering methods; bidirectional best-hit (BDBH), OrthoMCL (OMCL) or COGtriangles clustering algorithms (COG).

|    	Method   		 |                           	Definition                         		 |
|:---------------------:    |:---------------------------------------------------------------------:    |
| **Bidirectional Best-Hit (BDBH)**     |      	Clusters proteins by identifying reciprocal best hits between genomes.    		 |
|	**OrthoMCL (OMCL)**  	 | Uses graph theory to cluster proteins based on sequence similarity, handling paralogous genes and gene duplications.    |
|	**COGtriangles:**  	 |   	Assigns proteins to predefined functional categories (COGs) based on best matches to the COG database using a triangle inequality-based algorithm.  		 |

<a href="../fig/GET_HOMOLOGUES_flow_char.jpeg">
  <img src="../fig/GET_HOMOLOGUES_flow_char.jpeg" width="435" height="631" alt="GET_HOMOLOGUES flow chart." />
</a>

> ## Exercise 1: Clustering algorithms
>
> Complete the line blank with the correct clustering algorithms
>
> |------------------------------+------------------------------------------------------------------------------|  
> | **Algorithms**                  		 |     **Information required**                            		 |  
> |------------------------------+------------------------------------------------------------------------------|  
> | ___________________ |  Starting from a reference genome, keep adding genomes stepwise while storing the sequence clusters that result from merging the latest bidirectional best hits                         		 |  
> |------------------------------+------------------------------------------------------------------------------|  
> | ___________________ | Merges triangles of inter-genomic symmetrical best matches |   
> |------------------------------+------------------------------------------------------------------------------|  
> | ___________________ | uses the Markov Cluster Algorithm to group sequences, with inflation (-F) controlling cluster granularity  |  
> |------------------------------+------------------------------------------------------------------------------|
>
>
>> ## Solution
>>
>> |------------------------------+------------------------------------------------------------------------------|  
>> | **algorithms**                  		 |     **Information required**                            		 |  
>> |------------------------------+------------------------------------------------------------------------------|  
>> | BDBH             		 |  Starting from a reference genome, keep adding genomes stepwise while storing the sequence clusters that result from merging the latest bidirectional best hits                         		 |  
>> |------------------------------+------------------------------------------------------------------------------|  
>> | COGS  | Merges triangles of inter-genomic symmetrical best matches |   
>> |------------------------------+------------------------------------------------------------------------------|  
>> | OMCL    | uses the Markov Cluster Algorithm to group sequences, with inflation (-F) controlling cluster granularity  |  
>> |------------------------------+------------------------------------------------------------------------------|
>>
>>
> {: .solution}
{: .challenge}

## Considerations
Please make sure that you are in the Pangenomics environment.

~~~
$ conda deactivate
$ conda activate Pangenomics_Global  
~~~
{: .language-bash}
Make sure that get_homologues is installed

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

> ## Notes
> Get_homologues suggests that you run your data with a directory because you could add a new file *.gbk* in the future, if necessary.
{: .callout}


## Create a directory for get_homologues output
It's necessary to create a new folder to store all the results.

~~~
$ mkdir -p ~/pan_workshop/results/pangenome/get_homologues/data_gbks #Create directory (option -p create all parent directories)
$ cd  ~/pan_workshop/results/pangenome/get_homologues/data_gbks #Locates you in the directory 'data_gbks'
~~~
{: .language-bash}
We need to create symbolic links with all the *.gbk* files created with prokka or downloaded with ncbi
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

## Create the directory to run OMCL algorithm
~~~
$ cd ..
~~~
{: .language-bash}

Generate the clusters with OMCL (OMCL, PubMed=12952885)

Since the following command can take around 8 minutes to run we will use a screen session to run it. The screen session will not have the conda environment activated, so let’s activate it again.
~~~
screen -R clustering
conda activate Pangenomics_Global
~~~
{: .language-bash}
And now let's run our program.
~~~
get_homologues.pl -d data_gbks -M -t 0 -c -n 8
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

> ## Notes
If the option -e is added, the resulting clusters will contain only single-copy genes from each taxon, i.e. orthologues. This flag exclude clusters within paralogues. This is useful to make genome-level phylogenetic analyses in only single copy-genes.

However, it's important to note that using only orthologues may overlook genes that have undergone significant divergence or genes with species-specific functions. Including paralogues, which are duplicated genes within a genome, or genes with no clear orthologues can also be informative in certain pangenomic analyses, such as studying gene family expansions or specific adaptations within a species.

{: .callout}

## Describe your gene families in one table

Get_homologues gave us one FASTA file for each gene family, with the sequences of the genes included in the family. These files look like this:
~~~
$ head data_gbks_homologues/Streptococcusagalactiae18RS21prokka_f0_0taxa_algOMCL_e0_/1_IS30_family_transpos...faa
~~~
{: .language-bash}
~~~
>ID:GMBKAPON_00410 |[Streptococcus agalactiae]|18RS21|Streptococcus_agalactiae_18RS21_prokka.gbk|IS30 family transpos..|1167|AAJO01000202.1(1503):112-1278:1 ^COG:COG2826^ Streptococcus agalactiae strain 18RS21.|neighbours:start(),end()|neighbour_genes:start(),end()| | aligned:1-388 (388)
MTKHKHLTLLDRNDIQSGLDRGETFKAIGLNLLKHPTTIAKEVKRNKQLRESTKDCLDCPLLRKAPYVCNGCPKRRINCGYKKTFYLAKQAQRNYEKLLVESREGIPLNKETFWKIDRVLSNGVKKGQRIYHILKTNDLEVSSSTVYRHIKKGYLSITPIDLPRAVKFKKRRKSTLPPIPKAIKEGRRYEDFIEHMNQSELNSWLEMDTVIGRIGGKVLLTFNVAFCNFIFAKLMDSKTAIETAKHIQVIKRTLYDNKRDFFELFPVILTDNGGEFARVDDIEIDVCGQSQLFFCDPNRSDQKARIEKNHTLVRDILPKGTSFDNLTQEDINLALSHINSVKRQALNGKTAYELFSFTYGKDIASILGIEEITAEDVCQSPKLLKDKI
~~~
{: .output}

We want to create a file that summarizes the information of the clustering by showing only which genes correspond to which family. 
We will need that file in the next episode to explore our pangenome with another program.

To obtain this file, that we will name `gene_families.tsv`, we will extract the IDs of the genes from the FASTA headers(in the FASTA header we see the ID of the gene at the beggining after `ID:`) and the name of the families from the file names. For this we will use the following short script.  

Copy the contents and paste them in a file:
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

Now we have in only one file all the description of our clustering results!
