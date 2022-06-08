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
- "Implemented and intepreted the evolutionary history using Clustering orthologous proteins."

keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---

# Get_Homologues
GET_HOMOLOGUES: a versatile software package for pan-genome analysis is maintained by Bruno Contreras-Moreira (bcontreras at eead.csic.es) and Pablo Vinuesa (vinuesa at ccg.unam.mx). 
- Clustering protein and nucleotide sequences in homologous (possibly orthologous) groups, on the grounds of sequence similarity.
- Identification of orthologous groups of intergenic regions, flanked by orthologous open reading frames (ORFs), conserved across related genomes.
- Definition of pan- and core-genomes by calculation of overlapping sets of proteins.

## Considerations
Please ensure that you be in the environment Pangenomics
~~~
conda activate Pangenomics
~~~
{: .source}

## Step 1
We need to create a new folder
~~~
mkdir get_homologues_clus
cd get_homologues_clus
~~~
{: .source}

~~~
betterlab@132.248.196.38's password:
~~~
{: .output}

## Step 2
To generate the directory clusters with BDBH, this option is default.
~~~
get_homologues.pl -d datos/gbk
~~~
{: .source}

To generate the directory cluster with COG 
~~~
get_homologues.pl -d datos/gbk -G
~~~
{: .source}

To Generate the OMCL cluster directory
~~~
get_homologues.pl -d datos/gbk -M
~~~
{: .source}

## Step 3
Use orthoMCL algorithm (OMCL, PubMed=12952885)
~~~
get_homologues.pl -d /home/betterlab/GenomeMining/datos/gbk -e -M 
~~~
{: .source}

~~~
betterlab@132.248.196.38's password:
~~~
{: .output}

> ## Exercise 2: 
> 
> Complete the lines of code to obtain the required information
> What is the interpret the Venn diagrams?
> Which do you think the best clustering algorithms?

>> ## Solution
>> 
> {: .solution}
{: .challenge} 

{% include links.md %}
