---
title: "Get_Homologues"
teaching: 20 min
exercises: 5 min
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

# Get_Homologues a versatile software package for pan-genome analysis is maintained by Bruno Contreras-Moreira and Pablo Vinuesa. 
## Main Task
- Clustering protein and nucleotide sequences in homologous (possibly orthologous) groups, on the grounds of sequence similarity.
- Identification of orthologous groups of intergenic regions, flanked by orthologous open reading frames (ORFs), conserved across related genomes.
- Definition of pan- and core-genomes by calculation of overlapping sets of proteins.

## Considerations
Please ensure that you are in the environment of Pangenomics. You can omit his step if you have activated the environment.

~~~
conda activate Pangenomics
~~~
{: .source}
Now, We ensure that get_homologues is install
~~~
get_homologues.pl -h
~~~
{: .source}
{: .laguage-bash}

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

## Step 1. Generate a folder get_homologues
It's necessary that we create a new folder when all results are sent.
~~~
mkdir ~dc_workshop/results/pangenome/get_homologues
mkdir ~dc_workshop/results/pangenome/get_homologues/data_get
cd  ~dc_workshop/results/pangenome/get_homologues/data_get
~~~
{: .source}
We need to create a Symbolic link with the file *.gbk*
~~~
find ~/dc_workshop/results/annotated/. -name "*aga*_prokka.gbk*" -exec ln -s {} . ';'
ls ~dc_workshop/results/pangenome/get_homologues/data_get
~~~
{: .source}

~~~
agalactiae_18RS21_prokka.gbk  agalactiae_A909_prokka.gbk    agalactiae_COH1_prokka.gbk
agalactiae_515_prokka.gbk     agalactiae_CJB111_prokka.gbk  agalactiae_H36B_prokka.gbk
~~~
{: .output}

## Step 2. Generate the directory clusters
To generate the directory clusters with BDBH, this option is default.
~~~
get_homologues.pl -d dc_workshop/data/*/*.gbk
~~~
{: .source}
~~~
 /home/betterlab/.conda/envs/Pangenomics/bin/get_homologues.pl -i 0 -d data_get -o 0 -X 0 -e 0 -f 0 -r 0 -t all -c 0 -z 0 -I 0 -m local -n 2 -M 0 -G 0 -p 0 -C 75 -S 1 -E 1e-05 -F 1.5 -N 0 -B 50 -b 0 -s 0 -D 0 -g 0 -a '0' -x 0 -R 0 -A 0 -P 0

 version 28042022
 results_directory=/home/betterlab/dc_workshop/results/pangenome/get_homologues/data_get_homologues
 parameters: MAXEVALUEBLASTSEARCH=0.01 MAXPFAMSEQS=250 BATCHSIZE=100 KEEPSCNDHSPS=1
 diamond job:0

 checking input files...
 agalactiae_18RS21_prokka.gbk 1960
 agalactiae_515_prokka.gbk 2059
 agalactiae_A909_prokka.gbk 2067
 agalactiae_CJB111_prokka.gbk 2044
 agalactiae_COH1_prokka.gbk 2143
 agalactiae_H36B_prokka.gbk 2166

 6 genomes, 12439 sequences

 taxa considered = 6 sequences = 12439 residues = 3548657 MIN_BITSCORE_SIM = 18.3

 mask=agalactiae18RS21prokka_f0_alltaxa_algBDBH_e0_ (_algBDBH)

 running makeblastdb with /home/betterlab/dc_workshop/results/pangenome/get_homologues/data_get_homologues/agalactiae_18RS21_prokka.gbk.fasta

 running makeblastdb with /home/betterlab/dc_workshop/results/pangenome/get_homologues/data_get_homologues/agalactiae_515_prokka.gbk.fasta

 running makeblastdb with /home/betterlab/dc_workshop/results/pangenome/get_homologues/data_get_homologues/agalactiae_A909_prokka.gbk.fasta

 running makeblastdb with /home/betterlab/dc_workshop/results/pangenome/get_homologues/data_get_homologues/agalactiae_CJB111_prokka.gbk.fasta

 running makeblastdb with /home/betterlab/dc_workshop/results/pangenome/get_homologues/data_get_homologues/agalactiae_COH1_prokka.gbk.fasta

 running makeblastdb with /home/betterlab/dc_workshop/results/pangenome/get_homologues/data_get_homologues/agalactiae_H36B_prokka.gbk.fasta

 running BLAST searches ...
 done

 concatenating and sorting BLAST/DIAMOND results...
 sorting _agalactiae_18RS21_prokka.gbk results (2.8MB)
 sorting _agalactiae_515_prokka.gbk results (3MB)
 sorting _agalactiae_A909_prokka.gbk results (3MB)
 sorting _agalactiae_CJB111_prokka.gbk results (3MB)
 sorting _agalactiae_COH1_prokka.gbk results (3.2MB)
 sorting _agalactiae_H36B_prokka.gbk results (3.1MB)
 done


 parsing blast result! (/home/betterlab/dc_workshop/results/pangenome/get_homologues/data_get_homologues/tmp/all.blast , 18MB)
 parsing file finished

 creating indexes, this might take some time (lines=3.15e+05) ...

 construct_taxa_indexes: number of taxa found = 6
 number of file addresses/BLAST queries = 1.2e+04

 clustering orthologous sequences

 clustering inparalogues in agalactiae_18RS21_prokka.gbk (reference)
 162 sequences

 clustering inparalogues in agalactiae_515_prokka.gbk
 13 sequences

 finding BDBHs between agalactiae_18RS21_prokka.gbk and agalactiae_515_prokka.gbk
 1382 sequences

 clustering inparalogues in agalactiae_A909_prokka.gbk
 51 sequences

 finding BDBHs between agalactiae_18RS21_prokka.gbk and agalactiae_A909_prokka.gbk
 1455 sequences

 clustering inparalogues in agalactiae_CJB111_prokka.gbk
 23 sequences

 finding BDBHs between agalactiae_18RS21_prokka.gbk and agalactiae_CJB111_prokka.gbk
 1413 sequences

 clustering inparalogues in agalactiae_COH1_prokka.gbk
 60 sequences

 finding BDBHs between agalactiae_18RS21_prokka.gbk and agalactiae_COH1_prokka.gbk
 1369 sequences

 clustering inparalogues in agalactiae_H36B_prokka.gbk
 71 sequences

 finding BDBHs between agalactiae_18RS21_prokka.gbk and agalactiae_H36B_prokka.gbk
 1390 sequences

 looking for valid ORF clusters (n_of_taxa=6)...


 number_of_clusters = 1105
 cluster_list = data_get_homologues/agalactiae18RS21prokka_f0_alltaxa_algBDBH_e0_.cluster_list
 cluster_directory = data_get_homologues/agalactiae18RS21prokka_f0_alltaxa_algBDBH_e0_

 runtime: 840 wallclock secs (13.44 usr  0.24 sys + 593.88 cusr 10.63 csys = 618.19 CPU)
 RAM use: 65.8 MB
~~~
{: .output}

To generate the directory cluster with COG 

~~~
get_homologues.pl -d dc_workshop/data/*/*.gbk -G
~~~
{: .source}
To Generate the OMCL cluster director y(OMCL, PubMed=12952885)

~~~
get_homologues.pl -d dc_workshop/data/*/*.gbk -M
~~~
{: .source}

## Step 3. Compare all clusters from diferent algoritms
~~~
compare_clusters.pl -o alg_intersection -m -d\
gbk_homologues/A909_f0_alltaxa_algBDBH_e0_,\
gbk_homologues/A909_f0_alltaxa_algCOG_e0_,\
gbk_homologues/A909_f0_alltaxa_algOMCL_e0_
~~~
{: .source}

Use the scp protocol in order to see the venn diagram
~~~
scp user@ip:/path/to/file/venn_t0.pdf .
~~~
{: .source}

~~~
usuario@ip password:
~~~
{: .output}

search file in the file browser on your computer.

## Step 4. Obtaining a pangenome matrix


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
