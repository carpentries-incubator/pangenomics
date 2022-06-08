---
title: "Get_Homologues"
teaching: 0
exercises: 0
questions:
- "What is Get_Homologues?"
- "¿Que es Get_Homologues?"
- "¿Que es Clusterin protein?"
- "¿Que es un diagrama de Venn?"
- "¿Cual los son los algoritmos de clustering?"

objectives:
- "First learning objective. (FIXME)"
keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---

# Get_Homologues
- Clustering protein and nucleotide sequences in homologous (possibly orthologous) groups, on the grounds of sequence similarity.
- Identification of orthologous groups of intergenic regions, flanked by orthologous open reading frames (ORFs), conserved across related genomes.
- Definition of pan- and core-genomes by calculation of overlapping sets of proteins.

## Considerations
Please ensure that you be in the environment Pangenomics
~~~
Conda activate Pangenomics
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
Use clustering COGtriangle algorithm (COGS, PubMed=20439257)
~~~
get_homologues.pl -d /home/betterlab/GenomeMining/datos/gbk -e -G 
~~~
{: .source}

~~~
betterlab@132.248.196.38's password:
~~~
{: .output}
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

{% include links.md %}
