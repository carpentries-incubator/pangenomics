---
title: "Clustering Protein Families"
teaching: 20
exercises: 5
questions:
- "What is Get_Homologues?"
- "What is Clustering?"
- "Which are the clustering algorithms that use Get_Homologues?"

objectives:
- "Clustering orthologous proteins from GenBank files"
- "Create a Venn diagram using different clustering algorithms"
- "Implement and interpret the evolutionary history using Clustering orthologous proteins"

keypoints:
- "Get_Homologues is a software package for pangenome analyses "
- "Clustering refers to the task of grouping sequences in such a way that sequences in the same group are more similar to each other than to those in other groups "
- "Three sequence-clustering algorithms are supported by Get_Homologues; BDBH, OMCL and COGtriangles "
---

## What is Get_Homologues?

Get_Homologues is a versatile software package for pan-genome analysis, maintained by Bruno Contreras-Moreira and Pablo Vinuesa.
Its main task is grouping or clustering protein and nucleotide sequences in homologous (possibly orthologous) groups,
on the grounds of sequence similarity.
This software identifies orthologous groups of intergenic regions, flanked by orthologous open reading frames (ORFs),
conserved across related genomes.

Get_Homologues supports three sequence-clustering algorithms; bidirectional best-hit (BDBH), OrthoMCL (OMCL) or COGtriangles clustering algorithms.
BDBH algorithm uses one sequence from the reference genome and maintains the clusters growing.
OMCL algorithm groups nodes in a BLAST graph to build clusters.
And the COG algorithm requires a triangle of reciprocal hits and merges them.

The definition of pan- and core-genomes by Get_Homologues is done by calculation of overlapping sets of proteins.


<a href="../fig/GET_HOMOLOGUES_flow_char.jpeg">
  <img src="../fig/GET_HOMOLOGUES_flow_char.jpeg" width="435" height="631" alt="GET_HOMOLOGUES flow chart.
                                                                       		 Input files are either GenBank or FASTA and can produce different outputs.
                                                                       		 BLAST and Pfam searches are optimized for local and cluster computer environments.
                                                                       		 Once this sequences are sorted and indexed, one of the three clustering algorithms
                                                                       		 (BDBH, OMCL and COGtriangles) yields FASTA files of sequence clusters.
                                                                       		 This clusters could be flanked intergene clusters,
                                                                       		 pan/core-genome size estimates, pangenome matrices and syntenic clusters." />
</a>

## Considerations
Please ensure that you are in the Pangenomics environment. You can omit this step if you have already activated the environment.

~~~
$ conda activate Pangenomics_Global  
~~~
{: .language-bash}
Make sure that get_homologues is install

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


## Step 1. Generate a folder get_homologues
It's necessary to create a new folder where all the results will be sent.

~~~
$ mkdir -p ~/gm_workshop/results/pangenome/get_homologues/data_get #Create directory (-p create all parents)
$ cd  ~/gm_workshop/results/pangenome/get_homologues/data_get # Change to the directory 'data_get'
~~~
{: .language-bash}
We need to create a Symbolic link with the file *.gbk*
~~~
$ ln -s ~/gm_workshop/results/annotated/*Streptococcus_agalactiae_*[A-Z]*.prokka.gbk .
$ ls ~/gm_workshop/results/pangenome/get_homologues/data_get #List the symbolic link
~~~
{: .language-bash}

~~~
Streptococcus_agalactiae_18RS21.prokka.gbk  Streptococcus_agalactiae_COH1.prokka.gbk
Streptococcus_agalactiae_A909.prokka.gbk    Streptococcus_agalactiae_H36B.prokka.gbk
Streptococcus_agalactiae_CJB111.prokka.gbk
~~~
{: .output}

## Step 2. Generate the directory clusters
~~~
$ cd  ~/gm_workshop/results/pangenome/get_homologues/
~~~
{: .language-bash}

To generate the directory clusters with BDBH, this option is default. We use the -c flag to generate a report from core and pangenome.
~~~
$ get_homologues.pl -d data_get -c
~~~
{: .language-bash}

> ## Notes
> This script typically takes about 20 minutes to run.
{: .callout}

~~~
# number_of_clusters = 1177
# cluster_list = data_get_homologues/Streptococcusagalactiae18RS21_f0_alltaxa_algBDBH_e0_.cluster_list
# cluster_directory = data_get_homologues/Streptococcusagalactiae18RS21_f0_alltaxa_algBDBH_e0_

# runtime: 581 wallclock secs (19.27 usr  0.49 sys + 303.82 cusr  4.51 csys = 328.09 CPU)
# RAM use: 62.1 MB
~~~
{: .output}

To generate the directory cluster with COG

~~~
$ get_homologues.pl -d data_get -G
~~~
{: .language-bash}

~~~
# number_of_clusters = 1180
# cluster_list = data_get_homologues/Streptococcusagalactiae18RS21_f0_alltaxa_algCOG_e0_.cluster_list
# cluster_directory = data_get_homologues/Streptococcusagalactiae18RS21_f0_alltaxa_algCOG_e0_

# runtime: 15 wallclock secs ( 0.82 usr  0.05 sys +  1.79 cusr  0.52 csys =  3.18 CPU)
# RAM use: 50.9 MB
~~~
{: .output}

To Generate the OMCL cluster directory (OMCL, PubMed=12952885)

~~~
$ get_homologues.pl -d data_get -M
~~~
{: .language-bash}

~~~
# number_of_clusters = 1181
# cluster_list = data_get_homologues/Streptococcusagalactiae18RS21_f0_alltaxa_algOMCL_e0_.cluster_list
# cluster_directory = data_get_homologues/Streptococcusagalactiae18RS21_f0_alltaxa_algOMCL_e0_

# runtime:  5 wallclock secs ( 2.02 usr  0.11 sys +  0.49 cusr  0.24 csys =  2.86 CPU)
# RAM use: 55.4 MB
~~~
{: .output}

> ## Notes
If the option -e is added, the resulting clusters will contain only single-copy genes from each taxon, i.e. the orthologues. This flag forms singleton clusters, which are created when you exclude clusters within paralogues. This is useful to make genome-level phylogenetic analyses in only single copy-genes.
{: .callout}


## Step 3. Compare all clusters from different algorithms

Get_Homologues algorithm, BDBH, takes by default the smallest genome from the set of genomes as reference genome
~~~
$ ls -d data_get_homologues/*alltaxa* #List the genome reference
~~~
{: .language-bash}

> ## Notes
> After typing the comma in the following command line make sure you don't leave a space after the comma. This could cause an error.
{: .callout}

~~~
$ compare_clusters.pl -o alg_intersection -d\
data_get_homologues/Streptococcusagalactiae18RS21_f0_alltaxa_algBDBH_e0_,\
data_get_homologues/Streptococcusagalactiae18RS21_f0_alltaxa_algCOG_e0_,\
data_get_homologues/Streptococcusagalactiae18RS21_f0_alltaxa_algOMCL_e0_
~~~
{: .language-bash}
~~~
# Venn diagram = alg_intersection/venn_t0.pdf alg_intersection/venn_t0.svg
# Venn region file: alg_intersection/unique_Streptococcusagalactiae18RS21_f0_alltaxa_algBDBH_e0_.venn_t0.txt (9)
# Venn region file: alg_intersection/unique_Streptococcusagalactiae18RS21_f0_alltaxa_algCOG_e0_.venn_t0.txt (28)
# Venn region file: alg_intersection/unique_Streptococcusagalactiae18RS21_f0_alltaxa_algOMCL_e0_.venn_t0.txt (14)
# Venn region file: alg_intersection/intersection_Streptococcusagalactiae18RS21_f0_alltaxa_algBDBH_e0__Streptococcusagalactiae18RS21_f0_alltaxa_algCOG_e0_.venn_t0.txt (1)
# Venn region file: alg_intersection/intersection_Streptococcusagalactiae18RS21_f0_alltaxa_algBDBH_e0__Streptococcusagalactiae18RS21_f0_alltaxa_algOMCL_e0_.venn_t0.txt (16)
# Venn region file: alg_intersection/intersection_Streptococcusagalactiae18RS21_f0_alltaxa_algCOG_e0__Streptococcusagalactiae18RS21_f0_alltaxa_algOMCL_e0_.venn_t0.txt (0)
~~~
{: .output}

Download the Venn diagram to your local machine to be able to see it.
~~~
$ scp user@bioinformatica.matmor.unam.mx:~/gm_workshop/results/pangenome/get_homologues/alg_intersection/*.svg .
~~~
{: .language-bash}


<a href="../fig/venn_t0_GET_HOMOLOGUES.svg">
  <img src="../fig/venn_t0_GET_HOMOLOGUES.svg" alt="Venn diagram of core genomes generated by GET_HOMOLOGUES clustering algorithms.
                                           		 1105 BDBH clusters are labeled in red, 1115 COG are marked in green and 1110 OMCL in blue.
                                           		 1077 consensus clusters were detected by three algorithms." />
</a>

> ## Exercise 1: Comparing clustering algorithms
>
> Explore one of the gene clusters that result from the intersection of all algorithms with grep command:
>~~~
> $ ls alg_intersection | grep clpX
>~~~
>{: .language-bash}
>
> * Why are these genes at the intersection?
> * Is this cluster gene essential for living?
> * What other gene could be present in this output folder?
>
>
>> ## Solution
>>clpX is a gene that encodes part of a protease found in mitochondria, which is essential for living. The reason why they are in the intersection folder is that these cluster genes belong to the core genome.
> {: .solution}
{: .challenge}

> ## Exercise 2: Clustering algorithms
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

## Step 4. Plot core- and pan-genome clustering from BDBH algorithm

Use plot_pancore-matrix.pl to plot the core and pan-genome

~~~
$ plot_pancore_matrix.pl -i data_get_homologues/pan_genome_algBDBH.tab
~~~
{: .language-bash}
~~~
# /opt/anaconda3/envs/Pangenomics_Global/bin/plot_pancore_matrix.pl -i data_get_homologues/pan_genome_algBDBH.tab -f core_Tettelin -F 0.80 -a
# outfiles: data_get_homologues/pan_genome_algBDBH.tab_core_Tettelin.log , data_get_homologues/pan_genome_algBDBH.tab_core_Tettelin.png , data_get_homologues/pan_genome_algBDBH.tab_core_Tettelin.pdf , data_get_homologues/pan_genome_algBDBH.tab_core_Tettelin.svg
~~~
{: .output}

Download the pan and core-genome plots to your local machine.
~~~
$ scp user@bioinformatica.matmor.unam.mx:~/gm_workshop/results/pangenome/get_homologues/data_get_homologues/*_genome_algBDBH.tab_core_Tettelin.png
~~~
{: .language-bash}

<a href="../fig/core_genome_algBDBH.tab_core_Tettelin.png">
  <img src="../fig/core_genome_algBDBH.tab_core_Tettelin.png" alt="Aquí va el texto que describe a la imagen." />
</a>

<a href="../fig/pan_genome_algBDBH.tab_core_Tettelin.png">
  <img src="../fig/pan_genome_algBDBH.tab_core_Tettelin.png" alt="Aquí va el texto que describe a la imagen." />
</a>

> ## Exercise 3:
>
> Add another genome of the Streptococcus family. Try with another S. agalactiae genomes which are in the annotated folder.
> It is required to make a symbolic path in our data_get directory:
> ~~~
> $ find ~/gm_workshop/results/annotated/. -name "*Streptococcus_agalactiae_[1-9]*.prokka.gbk*" -exec ln -s {} . ';'
> ~~~
> {: .language-bash}
> Now ask for clustering all gene sequences with the get_homologues.pl default algorithm
> ~~~
> $ get_homologues.pl -d data_get -c
> ~~~
> {: .language-bash}
> What do you think happens to the number of gene clusters?
> Does it increase or decrease?
>> ## Solution
>> As we can check in the output:
>> ~~~
>> # number_of_clusters = 1105
>> # runtime: 247 wallclock secs ( 9.44 usr  0.37 sys + 128.79 cusr  1.90 csys = 140.50 CPU)
>> # RAM use: 56.4 MB
>> ~~~
>> {: .output}
>>    
>> The number of clusters decreases from  1177 to 1105. This is because the number of genes shared by all the genomes, i.e. core genome, decreases as we add another
>> genome, while the pangenome increases. We can see this with the following command :
>> ~~~
>> $ less data_get_homologues/pan_genome_algBDBH.tab
>> ~~~
>> {: .language-bash}
>>
> {: .solution}
{: .challenge}  

## Step 5. Obtaining a pangenome matrix
Firstly, use the -t 0 option with COG ang OMCL algorithms to include all possible clusters, considering those which might not contain sequences from all input genomes (taxa)
~~~
$ get_homologues.pl -d data_get -t 0 -M
~~~
{: .language-bash}

~~~
# number_of_clusters = 3634
# cluster_list = data_get_homologues/Streptococcusagalactiae18RS21_f0_0taxa_algOMCL_e0_.cluster_list
# cluster_directory = data_get_homologues/Streptococcusagalactiae18RS21_f0_0taxa_algOMCL_e0_

# runtime:  5 wallclock secs ( 2.10 usr  0.12 sys +  0.61 cusr  0.32 csys =  3.15 CPU)
# RAM use: 60.3 MB
~~~
{: .output}

~~~
$ get_homologues.pl -d data_get -t 0 -G
~~~
{: .language-bash}
 
~~~
# number_of_clusters = 3632
# cluster_list = data_get_homologues/Streptococcusagalactiae18RS21_f0_0taxa_algCOG_e0_.cluster_list
# cluster_directory = data_get_homologues/Streptococcusagalactiae18RS21_f0_0taxa_algCOG_e0_

# runtime: 17 wallclock secs ( 1.39 usr  0.12 sys +  2.46 cusr  0.71 csys =  4.68 CPU)
# RAM use: 56.0 MB
~~~
{: .output}

Then use the option

~~~~
$ ls -d data_get_homologues/*0taxa* #list cluster directories of COG and OMCL
~~~~
{: .language-bash}

~~~
$ compare_clusters.pl -o alg_intersection -m -T -d\
data_get_homologues/Streptococcusagalactiae18RS21_f0_0taxa_algCOG_e0_,\
data_get_homologues/Streptococcusagalactiae18RS21_f0_0taxa_algOMCL_e0_
~~~
{: .language-bash}

~~~
# Venn diagram = alg_intersection/venn_t0.pdf alg_intersection/venn_t0.svg
# Venn region file: alg_intersection/unique_Streptococcusagalactiae18RS21_f0_0taxa_algCOG_e0_.venn_t0.txt (182)
# Venn region file: alg_intersection/unique_Streptococcusagalactiae18RS21_f0_0taxa_algOMCL_e0_.venn_t0.txt (186)
~~~
{: .output}

## Step 6. Create a cladogram with our data

Now create a cladogram with the file `pangenome_matrix_t0.phylip.ph`, which is one of the different versions of the same pangenome matrix. This version contains a tree in Newick format. Make sure that it is in Newick format and change the extension to visualize in microreact.org. Check with the head command:
~~~
$ head alg_intersection/pangenome_matrix_t0.phylip.ph
~~~
{: .language-bash}

Rename the file by replacing the extension .ph to .nwk with the mv command

~~~
$ mv alg_intersection/pangenome_matrix_t0.phylip.ph alg_intersection/pangenome_matrix_t0.phylip.nwk
~~~
{: .language-bash}

Download the cladogram and see on microreact

~~~
$ scp user@bioinformatica.matmor.unam.mx:~/gm_workshop/results/pangenome/get_homologues/alg_intersection/*nwk .
~~~
{: .language-bash}

<a href="../fig/tree.png">
  <img src="../fig/tree.png" alt="Cladogram as visualized in microreact.org representing the relationship between the six S. agalactiae strains." />
</a>

<a href="../fig/legend.png">
  <img src="../fig/legend.png" alt="Cladogram legend labeling the S. agalactiae strains 18RS21 (green), 515 (yellow), A909 (purple), CJB111 (red), COH1 (blue) and H36B (orange)." />
</a>





{% include links.md %}



