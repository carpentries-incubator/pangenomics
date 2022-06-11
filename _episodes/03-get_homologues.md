---
title: "Clustering Proteins and Nucleotides"
teaching: 20
exercises: 5 
questions:
- "What is Get_Homologues?"
- "What is Clustering?"
- "Which are the clustering algorithms that use Get_Homologues?"

objectives:
- "Clustering orthologous proteins from Gen Bank files."
- "Create a Venn diagram using diferents clustering algorithms."
- "Implemented and interpreted the evolutionary history using Clustering orthologous proteins."

keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---

Get_Homologues a versatile software package for pan-genome analysis is maintained by Bruno Contreras-Moreira and Pablo Vinuesa. 
## Introducction 
## Main Task
- Clustering protein and nucleotide sequences in homologous (possibly orthologous) groups, on the grounds of sequence similarity.
- Identification of orthologous groups of intergenic regions, flanked by orthologous open reading frames (ORFs), conserved across related genomes.
- Definition of pan- and core-genomes by calculation of overlapping sets of proteins.

<a href="../fig/GET_HOMOLOGUES_flow_char.jpeg">
  <img src="../fig/GET_HOMOLOGUES_flow_char.jpeg" width="435" height="631" alt="Aquí va el texto que describe a la imagen." />
</a>

## Considerations
Please ensure that you are in the environment of Pangenomics. You can omit his step if you have activated the environment.

~~~
$ conda activate Pangenomics_Global  
~~~
{: .language-bash}
Now, We ensure that get_homologues is install

~~~
$ get_homologues.pl -h #This command display the options
~~~
{: .language-bash}


~~~
-v print version, credits and checks installation
-d directory with input FASTA files ( .faa / .fna ),           (overrides -i,
   GenBank files ( .gbk ), 1 per genome, or a subdirectory      use of pre-clustered sequences
   ( subdir.clusters / subdir_ ) with pre-clustered sequences   ignores -c, -g)
   ( .faa / .fna ); allows for new files to be added later;    
   creates output folder named 'directory_homologues'
-i input amino acid FASTA file with [taxon names] in headers,  (required unless -d is set)
   creates output folder named 'file_homologues'

Optional parameters:
-o only run BLAST/Pfam searches and exit                       (useful to pre-compute searches)
-c report genome composition analysis                          (follows order in -I file if enforced,
                                                                ignores -r,-t,-e)
-R set random seed for genome composition analysis             (optional, requires -c, example -R 1234,
                                                                required for mixing -c with -c -a runs)
-m runmode [local|cluster|dryrun]                              (default local)
-n nb of threads for BLAST/HMMER/MCL in 'local' runmode        (default=2)
-I file with .faa/.gbk files in -d to be included              (takes all by default, requires -d)

Algorithms instead of default bidirectional best-hits (BDBH):
-G use COGtriangle algorithm (COGS, PubMed=20439257)           (requires 3+ genomes|taxa)
-M use orthoMCL algorithm (OMCL, PubMed=12952885)

Options that control sequence similarity searches:
-X use diamond instead of blastp                               (optional, set threads with -n)
-C min %coverage in BLAST pairwise alignments                  (range [1-100],default=75)
-E max E-value                                                 (default=1e-05,max=0.01)
-D require equal Pfam domain composition                       (best with -m cluster or -n threads)
   when defining similarity-based orthology
-S min %sequence identity in BLAST query/subj pairs            (range [1-100],default=1 [BDBH|OMCL])
-N min BLAST neighborhood correlation PubMed=18475320          (range [0,1],default=0 [BDBH|OMCL])
-b compile core-genome with minimum BLAST searches             (ignores -c [BDBH])

Options that control clustering:
-t report sequence clusters including at least t taxa          (default t=numberOfTaxa,
                                                                t=0 reports all clusters [OMCL|COGS])
-a report clusters of sequence features in GenBank files       (requires -d and .gbk files,
   instead of default 'CDS' GenBank features                    example -a 'tRNA,rRNA',
                                                                NOTE: uses blastn instead of blastp,
                                                                ignores -g,-D)
-g report clusters of intergenic sequences flanked by ORFs     (requires -d and .gbk files)
   in addition to default 'CDS' clusters
-f filter by %length difference within clusters                (range [1-100], by default sequence
                                                                length is not checked)
-r reference proteome .faa/.gbk file                           (by default takes file with
                                                                least sequences; with BDBH sets
                                                                first taxa to start adding genes)
-e exclude clusters with inparalogues                          (by default inparalogues are
                                                                included)
-x allow sequences in multiple COG clusters                    (by default sequences are allocated
                                                                to single clusters [COGS])
-F orthoMCL inflation value                                    (range [1-5], default=1.5 [OMCL])
-A calculate average identity of clustered sequences,          (optional, creates tab-separated matrix,
 by default uses blastp results but can use blastn with -a      recommended with -t 0 [OMCL|COGS])
-P calculate percentage of conserved proteins (POCP),          (optional, creates tab-separated matrix,
 by default uses blastp results but can use blastn with -a      recommended with -t 0 [OMCL|COGS])
-z add soft-core to genome composition analysis                (optional, requires -c [OMCL|COGS])
~~~
{: .output}

> ## Notes
> Get_homologues suggests that you run your data with a directory because you could add a new file *.gbk* in the future, if necessary.
{: .callout}


## Step 1. Generate a folder get_homologues
It's necessary that we create a new folder when all results are sent.

~~~
$ mkdir -p ~/gm_workshop/results/pangenome/get_homologues/data_get #Create directory (-p create all parents)
$ cd  ~/gm_workshop/results/pangenome/get_homologues/data_get # Change to the directory 'data_get'
~~~
{: .language-bash}
We need to create a Symbolic link with the file *.gbk*
~~~
$ find ~/gm_workshop/results/annotated/. -name "*aga*_prokka.gbk*" -exec ln -s {} . ';' 
$ ls ~/gm_workshop/results/pangenome/get_homologues/data_get #List the symbolic link
~~~
{: .language-bash}

~~~
agalactiae_18RS21_prokka.gbk  agalactiae_A909_prokka.gbk    agalactiae_COH1_prokka.gbk
agalactiae_515_prokka.gbk     agalactiae_CJB111_prokka.gbk  agalactiae_H36B_prokka.gbk
~~~
{: .output}

## Step 2. Generate the directory clusters
~~~
$ cd  ~/gm_workshop/results/pangenome/get_homologues/
~~~
{: .language-bash}

To generate the directory clusters with BDBH, this option is default.
~~~
$ get_homologues.pl -d data_get
~~~
{: .language-bash}

> ## Notes
> When run the script above typically takes about 20 minutes.
{: .callout}

~~~
# number_of_clusters = 1105
# cluster_list = data_get_homologues/agalactiae18RS21prokka_f0_alltaxa_algBDBH_e0_.cluster_list
# cluster_directory = data_get_homologues/agalactiae18RS21prokka_f0_alltaxa_algBDBH_e0_

# runtime: 771 wallclock secs ( 9.01 usr  0.21 sys + 419.58 cusr  6.00 csys = 434.80 CPU)
# RAM use: 61.1 MB
~~~
{: .output}

To generate the directory cluster with COG 

~~~
$ get_homologues.pl -d data_get -G
~~~
{: .language-bash}

~~~
# number_of_clusters = 1115
# cluster_list = data_get_homologues/agalactiae18RS21prokka_f0_alltaxa_algCOG_e0_.cluster_list
# cluster_directory = data_get_homologues/agalactiae18RS21prokka_f0_alltaxa_algCOG_e0_

# runtime: 16 wallclock secs ( 1.06 usr  0.08 sys +  2.50 cusr  0.68 csys =  4.32 CPU)
# RAM use: 55.4 MB
~~~
{: .output}

To Generate the OMCL cluster directory (OMCL, PubMed=12952885)

~~~
$ get_homologues.pl -d data_get -M
~~~
{: .source}

~~~
# number_of_clusters = 1110
# cluster_list = data_get_homologues/agalactiae18RS21prokka_f0_alltaxa_algOMCL_e0_.cluster_list
# cluster_directory = data_get_homologues/agalactiae18RS21prokka_f0_alltaxa_algOMCL_e0_

# runtime:  6 wallclock secs ( 3.13 usr  0.11 sys +  0.56 cusr  0.34 csys =  4.14 CPU)
# RAM use: 61.9 MB
~~~
{: .output}


## Step 3. Compare all clusters from different algorithms

Get_Homologues for default the algorithm BDBH takes from the set of genomes to the smallest as a reference genome
~~~
$ ls -d data_get_homologues/*alltaxa* #List the genome reference
~~~

> ## Notes
> After typing the comma in the following command line make sure you don't leave a space after the comma. This could cause an error.
{: .callout}

~~~
$ compare_clusters.pl -o alg_intersection -d\
data_get_homologues/agalactiae18RS21prokka_f0_alltaxa_algBDBH_e0_,\
data_get_homologues/agalactiae18RS21prokka_f0_alltaxa_algCOG_e0_,\
data_get_homologues/agalactiae18RS21prokka_f0_alltaxa_algOMCL_e0_
~~~
{: .language-bash}
~~~
# Venn diagram = alg_intersection/venn_t0.pdf alg_intersection/venn_t0.svg
# Venn region file: alg_intersection/unique_agalactiae18RS21prokka_f0_alltaxa_algBDBH_e0_.venn_t0.txt (8)
# Venn region file: alg_intersection/unique_agalactiae18RS21prokka_f0_alltaxa_algCOG_e0_.venn_t0.txt (31)
# Venn region file: alg_intersection/unique_agalactiae18RS21prokka_f0_alltaxa_algOMCL_e0_.venn_t0.txt (16)
# Venn region file: alg_intersection/intersection_agalactiae18RS21prokka_f0_alltaxa_algBDBH_e0__agalactiae18RS21prokka_f0_alltaxa_algCOG_e0_.venn_t0.txt (5)
# Venn region file: alg_intersection/intersection_agalactiae18RS21prokka_f0_alltaxa_algBDBH_e0__agalactiae18RS21prokka_f0_alltaxa_algOMCL_e0_.venn_t0.txt (15)
# Venn region file: alg_intersection/intersection_agalactiae18RS21prokka_f0_alltaxa_algCOG_e0__agalactiae18RS21prokka_f0_alltaxa_algOMCL_e0_.venn_t0.txt (2)
~~~
{: .output}

Use plot_pancore-matrix.pl to plot the core and pan-genome

~~~
$ plot_pancore_matrix.pl -i data_get_homologues/pan_genome_algCOG.tab
~~~
{: .language-bash}

Use the scp protocol in order to see the venn diagram
~~~
$ scp alumno6@bioinformatica.matmor.unam.mx:/home/alumno6/gm_workshop/results/pangenome/get_homologues/alg_intersection/*.svg .
~~~
{: .language-bash}

~~~
$ alumno6@bioinformatica.matmor.unam.mx's password:
~~~
{: .output}

Use the scp protocol in order to see the pan and core-genome size plot
~~~
$ scp alumno6@bioinformatica.matmor.unam.mx:/home/alumno6/gm_workshop/results/pangenome/get_homologues/alg_intersection/pan* .
~~~
{: .language-bash}
~~~
$ alumno6@bioinformatica.matmor.unam.mx's password:
~~~
{: .output}

search file in the file browser on your computer.

<a href="../fig/venn_t0_GET_HOMOLOGUES.svg">
  <img src="../fig/venn_t0_GET_HOMOLOGUES.svg" alt="Aquí va el texto que describe a la imagen." />
</a>

## Step 4. Obtaining a pangenome matrix
first we use the -t 0 option with COG ang OMCL alghortims to include all possible clusters, including those which might not contain sequences from all input genomes (taxa)
~~~
$ get_homologues.pl -d data_get -t 0 -M
~~~
{: .language-bash}

~~~
# number_of_clusters = 3634
# cluster_list = data_get_homologues/agalactiae18RS21prokka_f0_0taxa_algOMCL_e0_.cluster_list
# cluster_directory = data_get_homologues/agalactiae18RS21prokka_f0_0taxa_algOMCL_e0_

# runtime:  4 wallclock secs ( 1.30 usr  0.16 sys +  0.30 cusr  0.04 csys =  1.80 CPU)
# RAM use: 49.9 MB
~~~
{: .output}

~~~
$ get_homologues.pl -d data_get -t 0 -G
~~~
{: .language-bash}
 
~~~
# number_of_clusters = 3632
# cluster_list = data_get_homologues/agalactiae18RS21prokka_f0_0taxa_algCOG_e0_.cluster_list
# cluster_directory = data_get_homologues/agalactiae18RS21prokka_f0_0taxa_algCOG_e0_

# runtime:  4 wallclock secs ( 1.41 usr  0.10 sys +  0.30 cusr  0.03 csys =  1.84 CPU)
# RAM use: 55.5 MB
~~~
{: .output}

then we use the option 
~~~~
$ ls -d data_get_homologues/0taxa #list directories cluster COG and OMCL
~~~~
## Step 5. Compare only clustering algoriths COG and OMCL
~~~
$ compare_clusters.pl -o alg_intersection -m -d\ 
data_get_homologues/agalactiae18RS21prokka_f0_0taxa_algCOG_e0_,\
data_get_homologues/agalactiae18RS21prokka_f0_0taxa_algOMCL_e0_
~~~
{: .language-bash}

~~~ 
# Venn diagram = alg_intersection/venn_t0.pdf alg_intersection/venn_t0.svg
# Venn region file: alg_intersection/unique_agalactiae18RS21prokka_f0_0taxa_algCOG_e0_.venn_t0.txt (182)
# Venn region file: alg_intersection/unique_agalactiae18RS21prokka_f0_0taxa_algOMCL_e0_.venn_t0.txt (186)
~~~
{: .output}

<a href="../fig/venn_t0_GET_HOMOLOGUES.svg">
  <img src="../fig/venn_t0_GET_HOMOLOGUES.svg" alt="Aquí va el texto que describe a la imagen." />
</a>

<a href="../fig/venn_t0_COG_OMCL.svg">
  <img src="../fig/venn_t0_COG_OMCL.svg" alt="Aquí va el texto que describe a la imagen." />
</a>

> ## Exercise 1: 
> 
> What is the interpret the Venn diagrams?
>> ## Solution
>> 
> {: .solution}
{: .challenge} 

> ## Exercise 2: 
> 
> Complete the line blank with the correct clustering algorithms
> 
> |------------------------------+------------------------------------------------------------------------------|  
> | **algorithms**                           |     **Information required**                                     |  
> |------------------------------+------------------------------------------------------------------------------|  
> | ___________________ |  Starting from a reference genome, keep adding genomes stepwise while storing the sequence clusters that result of merging the latest bidirectional best hits                                  |  
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
>> | **algorithms**                           |     **Information required**                                     |  
>> |------------------------------+------------------------------------------------------------------------------|  
>> | BDBH                      |  Starting from a reference genome, keep adding genomes stepwise while storing the sequence clusters that result of merging the latest bidirectional best hits                                  |  
>> |------------------------------+------------------------------------------------------------------------------|  
>> | COGS  | Merges triangles of inter-genomic symmetrical best matches |   
>> |------------------------------+------------------------------------------------------------------------------|  
>> | OMCL    | uses the Markov Cluster Algorithm to group sequences, with inflation (-F) controlling cluster granularity  |  
>> |------------------------------+------------------------------------------------------------------------------| 
>> 
>>
> {: .solution}
{: .challenge} 



{% include links.md %}
