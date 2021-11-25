---
title: "Anvi'o"
teaching: 15 min
exercises: 40 min
questions:
- "What is Anvi'o?"
- "How can I obtain a pangenome analysis in Anvi'o?"
objectives:
- "Establish a dataset of genomes to obtain their pangenome"
- "Perform a basic workflow to obtain a pangenome in Anvi'o"
- "Understand and interpret the pangenome results"
keypoints:
- "Anvi’o is an open-source, community-driven analysis and visualization platform for microbial ‘omics. "
---
## Anvi'o

Anvi’o is an open-source, community-driven analysis and visualization platform for microbial ‘omics.
It brings together many aspects of today's cutting-edge strategies including **genomics, metagenomics, metatranscriptomics, phylogenomics, microbial population genetis, pangenomics and metapangenomis** in an *integrated* and *easy-to-use* fashion thorugh extensive interactive visualization capabilities. 

![Figure 1. Anvi'o network representation](../fig/anvio-network.png)

## Installation of Anvi'o

* Install a stable release of anvi’o on a Mac, Linux, or Windows running computer (best option for end users).
* Use anvi’o from the active development branch (best option for developers).
* Run anvi’o through its Docker containers without any installation (best option for the lazy).

## Selection of genomes 

Selection of proper data to construct a pangenome analysis is a crucial task of the process. By definition, a pangenome represents the entire set of genes from all strains in a clade. Thus, let's define and understand what a clade is. 

A **clade** (*kládos*, 'branch'), also known as *monophyletic group*, is a group of organisms that share the same common ancestor. In a phylogenetic tree, clades can be distinguished by identified the common ancestor of branch subgroups. For instance, in Fig.2, the blue and red subgroups are considered clades but green subgruop is not a clade, it is a *parahyletic group*. Why is that? Becuse it excludes the blue clade which has descended from the same common ancestor. Instead, the green and blue subgroups together form a clade.

![Figure 2. Cladogram representation](../fig/Cladogram.png)

For this lesson purposes, we selected 20 representative *Mycobacterium tuberculosis* genomes. 

## Brief description of *Mycobacterium tuberculosis* dataset




{% include links.md %}
