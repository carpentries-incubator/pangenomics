---
title: "Introduction to Pangenomics"
teaching: 10
exercises: 10
questions:
- "What is the origin of the term **pangenome**"
- "What kind of biological data are analyzed in a pangenome?"
- "Which factors affect the structure of a pangenome?"
objectives:
- "Understand the concept of pangenome"
- "Define the different subsets of the pangenome"
- "Learn which biological data are integrated in a pangenome analysis"
- "Discuss applications and questions that could be solved with a pangenome analysis"

keypoints:
- "A pangenome is the entire set of genes from all strains of a clade."
- "The principal subsets of the pangenome are: core, shell, and dispensable genome."
- "An *open pangenome* occurs when the size of the pangenome increases considerably with every added genome."
- "A *closed pangenome* occurs when only few gene families are incorporated to the pangenome when a new genome is added." 
- "A pangenome can be obtained by comparing the complete genome sequences of all members of a clade."
---
## What is a pangenome?

The term **pangenome** derives from the Greek *pan*, meaning 'whole' or 'everything', while *genome* describes 
an organism's complete genetic material. In 2005, Tettelin *et al.* used the term pangenome by the first time 
to describe the entire collection of genes from a group of pathogenic bacteria, named *Streptococcus agalactiae*. 
Nowadays, the pangenome concept has been also applied to other organisms including plants, animals, fungi, virus, etc. 
To summarize, the term pangenome is defined as the entire set of genes from all strains of a clade.

The pangenome can be divided in three principal components: **core genome** which contains genes present 
in all the strains compared, **shell genome** a proportion of genes absent in one or more strains 
and **dispensable genome** comprised by genes that are unique to each strain. 

There are two classes of pangenomes: open and closed. An **open pangenome** occurs when the size of the 
pangenome increases considerably every added genome (e.g. *Escherichia coli*). Whereas, a **closed pangenome** 
results when only few gene families are incorporated to the pangenome when a new member is added.


<a href="{{ page.root }}/fig/Characteristics_of_open_and_closed_pangenomes.png">
   <img src="{{ page.root }}/fig/Characteristics_of_open_and_closed_pangenomes.png" alt="Figure 1. Venn diagram representation of a) a closed pangenome and b) an open pangenome, comparing the sizes of their core and accesory genomes. c) Graphic depicting the differences of closed and open pangenomes regarding their size, total genes in pangenome, in relation with the number of sequenced genomes." />
  </a>

## How to select a group of genomes to construct a pangenome?

Selection of proper data to construct a pangenome analysis is a crucial task of the process. By definition, a pangenome represents 
the entire set of genes from all strains in a clade. Thus, let's define and understand what a clade is. 

A **clade** (*kládos*, 'branch'), also known as a *monophyletic group*, is a group of organisms that share the same common ancestor. 
In a phylogenetic tree, clades can be distinguished by identifying the common ancestor of branch subgroups. For instance, in the 
following figure, the blue and red subgroups are considered clades whereas the green subgroup is not a clade, it is a *parahyletic group*. 
Why is that? Because it excludes the blue clade which has descended from the same common ancestor. Nonetheless, the green and blue subgroups 
together form a clade.

<a href="{{ page.root }}/fig/Cladogram.png">
   <img src="{{ page.root }}/fig/Cladogram.png" alt="Figure 2. Cladogram representation indicating in blue and red two clades or monophyletic groups and in green a paraphyletic group" />
  </a>

## Brief description of *Streptococcus agalactiae* dataset

For this lesson purposes, we selected six *Streptococcus agalactiae* genomes. 


**Selected genomes**


| Organism                | Host    | Serotype   | Country     | SRA accesion number |
|-------------------------|---------|------------|-------------|---------------------|
|*S. agalactiae*  18RS21  | Human   |            |             |                     |
|*S. agalactiae*  515     | Human   |            |             |                     |
|*S. agalactiae*  A909    | Human   |            |             |                     |
|*S. agalactiae*  CJB111  | Human   |            |             |                     |
|*S. agalactiae*  COH1    | Human   |            |             |                     |
|*S. agalactiae*  H36B    | Human   |            |             |                     |



> ## Prepare your genome database 
> Make sure you have the six genomes previously described. If you do not have them yet, you can download them with the following instruction
> 
> ~~~
> $ cd ~
> $ wget https://zenodo.org/record/6636296/files/gm_workshop.zip
> $ unzip gm_workshop.zip
> ~~~
> {: .language-bash}
> 
> ~~~
>  ... 
>   inflating: dc_workshop/data/COH1/Streptococcus_agalactiae_COH1.gbk  
>   creating: dc_workshop/data/H36B/
>  inflating: dc_workshop/data/H36B/Streptococcus_agalactiae_H36B.fasta  
>  inflating: dc_workshop/data/H36B/Streptococcus_agalactiae_H36B.gbk  
> ~~~
> {: .output}
{: .checklist}

{% include links.md %}


