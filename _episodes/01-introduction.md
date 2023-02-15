---
title: "Introduction to Pangenomics"
teaching: 10
exercises: 10
questions:
- "What is a **pangenome**?" 
- "What kind of biological data are analyzed in a pangenome?"
- "Which factors affect the structure of a pangenome?"
objectives:
- "Understand the concept of pangenome."
- "Define the different subsets of the pangenome."
- "Learn which biological data are integrated into a pangenome analysis."
- "Discuss applications and questions that could be solved with a pangenome analysis."

keypoints:
- "A pangenome is the entire set of genes from all strains of a clade."
- "The principal subsets of the pangenome are core, shell, and dispensable genome."
- "An *open pangenome* occurs when the size of the pangenome increases considerably with every added genome."
- "A *closed pangenome* occurs when only a few gene families are incorporated to the pangenome when a new genome is added."
- "A pangenome can be obtained by comparing the complete genome sequences of all clade members."
---
## What is a pangenome?

The term **pangenome** derives from the Greek *pan*, meaning 'whole' or 'everything', while *genome* describes
an organism's complete genetic material. In 2005, Tettelin *et al.* used the name pangenome for the first time
to explain the entire collection of genes from a group of pathogenic bacteria named *_Streptococcus agalactiae_*.
The pangenome concept has been applied to other organisms, including plants, animals, fungi, viruses, etc.
To summarize, a pangenome is the entire set of genes from all clade strains.

The pangenome can be divided into three principal components: **core genome**, which contains genes present
in all the strains compared, **shell genome** is a proportion of genes absent in one or more strains
and **dispensable genome** comprised of genes unique to each strain.

There are two classes of pangenomes: open and closed. An **open pangenome** occurs when the number of gene families increases considerably for every added genome (e.g., *Escherichia coli*). Whereas, a **closed pangenome**
results when only a few gene families are incorporated into the pangenome when a new member is added.


<a href="{{ page.root }}/fig/01-01-01.png">
   <img src="{{ page.root }}/fig/01-01-01.png" alt=" Venn diagram of a) a closed pangenome and b) an open pangenome, comparing the sizes of their core and accessory genomes. c) Graphic depicting the differences between closed and open pangenomes regarding their size, total genes in pangenome, and the number of sequenced genomes." />
  </a>

## How to select a group of genomes to construct a pangenome?

Selecting proper data to construct a pangenome analysis is a crucial process task. A pangenome represents
the entire set of genes from all strains in a clade. Thus, let's define and understand what a clade is.

A **clade** (*kládos*, 'branch'), also known as a *monophyletic group*, is a group of organisms with the same common ancestor.
In a phylogenetic tree, clades can be distinguished by identifying the common ancestor of branch subgroups. For instance, the blue and red subgroups are considered clades in the following figure, whereas the green subgroup is not a clade. It is a *paraphyletic group*.
Why is that? Because it excludes the blue clade, which has descended from the same common ancestor. Nonetheless, the green and blue subgroups
together form a clade.

<a href="{{ page.root }}/fig/01-01-02.png">
   <img src="{{ page.root }}/fig/01-01-02.png" alt=" Cladogram indicating in blue and red two clades or monophyletic groups and green a paraphyletic group" />
  </a>



> ## Exercise 1: Select strains for a genus pangenome  
>  The following figure shows a simplification of the phylogenetic relationships in the genus _Streptococcus_ given in the article [Comparative Genomics of the Bacterial Genus Streptococcus Illuminates Evolutionary Implications of Species Groups](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0101229) by Gao et al. (2014). Observe the figure and discuss which strains you would select to construct a Pangenome that includes _S. agalactiae_ among several other species? 
>  
>  <a href="{{ page.root }}/fig/01-01-03.png"><img src="{{ page.root }}/fig/01-01-03.png" alt="Cladogram indicating _Streptococcus_ phylogeny" /></a>
>  A) Only _S. agalactiae_ can be chosen  
>  B) _S. urinalis_, _S. canis_, _S.poricuns_  
>  C)  _S. urinalis_ ,_S. agalactiae_ 
>  
> > ## Solution
> > There are several ways to include _S. agalactiae_ in a more comprehensive pangenome. We will discuss some of these possibilities. 
> > A) No. A genus pangenome is possible; it must include all genomes available for the species in the genera, not only _S. agalactie_   
> > B) No. Chosing _S. urinalis_, _S. canis_, _S.poricuns_ will leave out many species that we should include to avoid parafileptic  
> > C) Yes. This selection will be a small pangenome but does not contain a paraphyletic group.  
> > 
> {: .solution}
{: .challenge}

## The *Streptococcus agalactiae* dataset

For this lesson, we will work with the six *Streptococcus agalactiae* strains included in the first pangenome made by *Tettelin et al., 2005*. The genomes of strains 18RS21 and H36B are already in our `pan_workshop/data`, but the others will be downloaded and annotated in the following episodes.


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
> Make sure you have the six genomes previously described. If you do not have them, you can download them with the following instructions.
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





