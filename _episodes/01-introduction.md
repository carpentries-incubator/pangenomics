---
title: "Introduction"
teaching: 0
exercises: 0
questions:
- "What is a pangenome?"
- "What are the principal componenets of a pangenome?"
- "Which means 'open pangenome' and 'closed pangenome'?"
- "How can I perform a pangenome analysis?"

objectives:
- "To learn to perform a pangenome analysis by using different bioinformatic tools"
keypoints:
- "The pangenome is the entire set of genes from all strains of a clade."
- "The principal components of a pangenome are: 'core', 'shell', and 'dispensable'"
- "Open pangenome occurs when the size of the pangenome increases considerably every added genome (p. e. *Escherichia coli*)"
- "Closed pangenome occurs when only few gene families are incorporated to the pangenome when a new member of the clade is added." 
- "A pangenome can be obtained by gathering and computing the complete genome sequence of the members of a clade by using specialized bioinformatic tools"

---

## What is a pangenome?

The term **pangenome** derives from the Greek *pan*, meaning 'whole' or 'everything', while *genome* describes an organism's complete genetic material. In 2005, Tettelin et al. used the term pangenome by first time to describe the entire collection of genes from a group of bacteria. Currently, this term has been used for other grups of organisms, including plants, animals, fungi, virus, etc. Today, the pangenome is defined as the entire set of genes from all strains of a clade.

The pangenome could be divided in three principal components: **core genome** which contains genes present in all the strains compared, **shell genome** a proportion of genes absent in one or more strains and **dispensable genome** comprised by genes that are unique to each strain. 

There are two classes of pangenomes: open and closed. An **open pangenome** occurs when the size of the pangenome increases considerably every added genome (p. e. *Escherichia coli*). Whereas, a **closed pangenome** results when only few gene families are incorporated to the pangenome when a new member is added.


![Figure 1. Characteristics of open and closed pangenomes](../fig/Characteristics_of_open_and_closed_pangenomes.png)



## How to select a group of genomes to construct a pangenome?

Selection of proper data to construct a pangenome analysis is a crucial task of the process. By definition, a pangenome represents the entire set of genes from all strains in a clade. Thus, let's define and understand what a clade is. 

A **clade** (*kládos*, 'branch'), also known as *monophyletic group*, is a group of organisms that share the same common ancestor. In a phylogenetic tree, clades can be distinguished by identified the common ancestor of branch subgroups. For instance, in Fig.2, the blue and red subgroups are considered clades but green subgruop is not a clade, it is a *parahyletic group*. Why is that? Becuse it excludes the blue clade which has descended from the same common ancestor. Instead, the green and blue subgroups together form a clade.



![Figure 2. Cladogram representation](../fig/Cladogram.png)





## Brief description of *Mycobacterium tuberculosis* dataset

For this lesson purposes, we selected 20 representative *Mycobacterium tuberculosis* genomes. 

Tuberculosis (TB) is the most lethal infection disease worldwide with 1.4 million causalities annually. Members of the *Mycobacterium tuberculosis* complex (MTBC) are recognized as the causative agents of TB in humans and other mammals. The MTBC is composed by nine *Mycobacterium* species affecting wild and domesticated animals, and two species that infect principally humans *Mycobacterium tuberculosis sensu stricto* and *Mycobacterium africanum*. Interestingly, all MTBC members form a clade when compared with other *Mycobacterium* spp. They share 99.99% of their nucleotide sequences and it has been decreted a maximum distance of ~2,500 SNPs between any pair of genomes belonging to the MTBC. No horizontal gene transfer or genomic recombination have been reported in these complex of bacteria, suggesting that these organisms are 'clones'. A pangenome analysis of these bacteria can help us to determine if MTBC has an open or a closed pengenome. 


**Selected genomes**


| Organism                | Host    | Lineage    | Country     | SRA accesion number |
|-------------------------|---------|------------|-------------|---------------------|
|*M. tuberculosis*  N0004 | Human   | L3         |USA          | ERR2704697          |
|*M. tuberculosis*  N0031 | Human   | L2         |USA          | ERR2704676          |
|*M. tuberculosis*  N0052 | Human   | L2         |USA          | ERR2704699          |
|*M. tuberculosis*  N0054 | Human   | L3         |USA          | ERR2704678          |
|*M. tuberculosis*  N0069 | Human   | L1         |China        | ERR2704679          |
|*M. tuberculosis*  N0072 | Human   | L1         |USA          | ERR2704680          |
|*M. tuberculosis*  N0091 | Human   | L6         |Gambia       | ERR2704681          |
|*M. tuberculosis*  N0136 | Human   | L4         |USA          | ERR2704700          |
|*M. tuberculosis*  N0145 | Human   | L2         |USA          | ERR2704702          |
|*M. tuberculosis*  N0155 | Human   | L2         |USA          | ERR2704684          |
|*M. tuberculosis*  N0157 | Human   | L1         |USA          | ERR2704704          |
|*M. tuberculosis*  N1176 | Human   | L5         |Ghana        | ERR2704686          |
|*M. tuberculosis*  N1201 | Human   | L6         |Ghana        | ERR2704687          |
|*M. tuberculosis*  N1202 | Human   | L6         |Ghana        | ERR2704688          |
|*M. tuberculosis*  N1216 | Human   | L4         |Ghana        | ERR2704689          |
|*M. tuberculosis*  N1268 | Human   | L5         |Sierra Leone | ERR2704690          |
|*M. tuberculosis*  N1272 | Human   | L5         |Germany      | ERR2704692          |
|*M. tuberculosis*  N1274 | Human   | L3         |Germany      | ERR2704693          |
|*M. tuberculosis*  N1283 | Human   | L4         |Germany      | ERR2704709          |
|*M. tuberculosis*  N3913 | Human   | L7         |Australia    | ERR2704711          |


{% include links.md %}

