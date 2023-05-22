---
title: "Clustering Protein Families"
teaching: 20
exercises: 5
questions:
- "Which clustering algorithm to choose and why?"
- "What is Get_Homologues?"

objectives:
- "Cluster orthologous proteins from GenBank (.gbk) files"
- "Explore protein clusters from OrthoMCL algorithm"
- "Understand basic pangenomics metrics"

keypoints:
- "Clustering refers to the task of grouping sequences in which the same group are more similar to each other than to those in other groups "
- "Get_Homologues is a software package for pangenomic analysis"
- "Three sequence-clustering algorithms are supported by Get_Homologues; BDBH, OrthoMCL and COGtriangles "
---

## Which clustering algorithm to choose and why?

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
$ ln -s ~/pan_workshop/results/annotated/Streptococcus_agalactiae_*_prokka.gbk .
$ ls ~/pan_workshop/results/pangenome/get_homologues/data_gbks #List the symbolic links
~~~
{: .language-bash}

~~~
Streptococcus_agalactiae_18RS21.prokka.gbk  Streptococcus_agalactiae_COH1.prokka.gbk
Streptococcus_agalactiae_A909.prokka.gbk    Streptococcus_agalactiae_H36B.prokka.gbk
Streptococcus_agalactiae_CJB111.prokka.gbk
~~~
{: .output}

## Create the directory to run OMCL algorithm
~~~
$ cd  ~/pan_workshop/results/pangenome/get_homologues/
~~~
{: .language-bash}

Generate the clusters with OMCL (OMCL, PubMed=12952885)

~~~
$ nano get_homologues_omcl.sh
~~~
{: .language-bash}

~~~
echo "Job started at:" && date
get_homologues.pl -d data_gbks -M -t 0 -c -n 8
echo "Job finished successfully at:" && date
~~~
{: .language-bash}

~~~
$ chmod +x get_homologues_omcl.sh
$ ./get_homologues_omcl.sh #This script can take up to 10 minutes
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

~~~
