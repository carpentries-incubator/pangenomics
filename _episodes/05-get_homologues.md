---
title: "Get_Homologues"
teaching: 20 min
exercises: 5 min
questions:
- "What is Get_Homologues?"
- "What is Clustering?"
- "Which are the clustering algorithms tha use Get_Homologues?"

objectives:
- "Clustering orthologous proteins from Gen Bank files."	
- Create a Venn diagram using diferents clustering algorithms."
- "Implemented and interpreted the evolutionary history using Clustering orthologous proteins."

keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---

# Get_Homologues a versatile software package for pan-genome analysis is maintained by Bruno Contreras-Moreira and Pablo Vinuesa . 
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

## Step 1
It's necessary that we create a new folder when all results are sent.
~~~
mkdir ~dc_workshop/results/pangenome/get_homologues
cd ~dc_workshop/results/pangenome/get_homologues
~~~
{: .source}

## Step 2 Generate the directory clusters
To generate the directory clusters with BDBH, this option is default.
~~~
get_homologues.pl -d dc_workshop/data/*/*.gbk
~~~
{: .source}

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

## Step 3 Compare all clusters from diferent algoritms
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
psw:##
search file in the file browser on your computer.


> ## Exercise 2: 
> 
> Complete the lines of code to obtain the required information
>> ## Solution
>> 
> {: .solution}
{: .challenge} 

{% include links.md %}
