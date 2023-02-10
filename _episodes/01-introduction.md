---
title: "Introduction to Pangenomics"
teaching: 10
exercises: 10
questions:
- "What is a **pangenome**?" 
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
- "A *closed pangenome* occurs when only a few gene families are incorporated to the pangenome when a new genome is added."
- "A pangenome can be obtained by comparing the complete genome sequences of all members of a clade."
---
## What is a pangenome?

The term **pangenome** derives from the Greek *pan*, meaning 'whole' or 'everything', while *genome* describes
an organism's complete genetic material. In 2005, Tettelin *et al.* used the term pangenome by the first time
to describe the entire collection of genes from a group of pathogenic bacteria, named *Streptococcus agalactiae*.
Nowadays, the pangenome concept has been also applied to other organisms including plants, animals, fungi, viruses, etc.
To summarize, the term pangenome is defined as the entire set of genes from all strains of a clade.

The pangenome can be divided in three principal components: **core genome** which contains genes present
in all the strains compared, **shell genome** a proportion of genes absent in one or more strains
and **dispensable genome** comprised by genes that are unique to each strain.

There are two classes of pangenomes: open and closed. An **open pangenome** occurs when the size of the
pangenome increases considerably every added genome (e.g. *Escherichia coli*). Whereas, a **closed pangenome**
results when only few gene families are incorporated to the pangenome when a new member is added.


<a href="{{ page.root }}/fig/01-01-01.png">
   <img src="{{ page.root }}/fig/01-01-01.png" alt="Venn diagram of a) a closed pangenome and b) an open pangenome, comparing the sizes of their core and accessory genomes. c) Graphic depicting the differences of closed and open pangenomes regarding their size, total genes in pangenome, in relation with the number of sequenced genomes." />
  </a>

## How to select a group of genomes to construct a pangenome?

Selection of proper data to construct a pangenome analysis is a crucial task of the process. By definition, a pangenome represents
the entire set of genes from all strains in a clade. Thus, let's define and understand what a clade is.

A **clade** (*kládos*, 'branch'), also known as a *monophyletic group*, is a group of organisms that share the same common ancestor.
In a phylogenetic tree, clades can be distinguished by identifying the common ancestor of branch subgroups. For instance, in the
following figure, the blue and red subgroups are considered clades whereas the green subgroup is not a clade, it is a *paraphyletic group*.
Why is that? Because it excludes the blue clade which has descended from the same common ancestor. Nonetheless, the green and blue subgroups
together form a clade.

<a href="{{ page.root }}/fig/01-01-02.png">
   <img src="{{ page.root }}/fig/01-01-02.png" alt="Cladogram indicating in blue and red two clades or monophyletic groups and in green a paraphyletic group" />
  </a>

## The *Streptococcus agalactiae* dataset

For this lesson, we will be working with the six *Streptococcus agalactiae* strains that were included in the first pangenome made by *Tettelin et al., 2005*. The genomes of the strains 18RS21 and H36B are already in our `pan_workshop/data`, but the other ones will be downloaded and annotated in the following episodes.


**Metadata**


|Strain	| Host	| Serotype   |
|-------------------------|---------|------------|
|*S. agalactiae*  18RS21  | Human   | II       	|
|*S. agalactiae*  515 	| Human   | Ia       	|
|*S. agalactiae*  A909	| Human   | Ia       	|
|*S. agalactiae*  CJB111  | Human   | V       	|
|*S. agalactiae*  COH1	| Human   | III       	|
|*S. agalactiae*  H36B	| Human   | Ib       	|



> ## Prepare your genome database
> Make sure you have the six genomes previously described. If you do not have them yet, you can download them with the following instruction
>
> ~~~
> $ cd ~
> $ wget https://zenodo.org/record/7620704/files/pan_workshop.zip?download=1
> $ unzip 'pan_workshop.zip?download=1'
> $ rm 'pan_workshop.zip?download=1'
> ~~~
> {: .language-bash}
{: .checklist}

{% include links.md %}





