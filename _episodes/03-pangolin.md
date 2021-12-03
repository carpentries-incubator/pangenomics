---
title: "PPanGGOLiN"
teaching: 20
exercises: 30
questions:
- "What is PPanGGOLiN?"
- "Which kind of pangenome analysis can be obtained with PPanGGOLiN?"
- "How to run a PPanGGOLiN analysis?"

objectives:
- "Understand the fundaments of the PPanGGOLiN tool."
- "Identify the principal differences between PPanGGOLiN and other pangenome tools"
- "Conduct a basic workflow with PPanGGOLiN"
keypoints:

---

PPanGGOLiN
===============================================
**Partitioned PanGenome Graph Of Linked Neighbors**


PPanGGOLiN is a software suite used to create and manipulate prokaryotic pangenomes from a set of either genomic DNA sequences or provided genome annotations. It is designed to scale up to tens of thousands of genomes. It has the specificity to partition the pangenome using a statistical approach rather than using fixed thresholds which gives it the ability to work with low-quality data such as Metagenomic Assembled Genomes (MAGs) or Single-cell Amplified Genomes (SAGs) thus taking advantage of large scale environmental studies and letting users study the pangenome of uncultivable species.

PPanGGOLiN builds pangenomes through a graphical model and a statistical method to partition gene families in persistent, shell and cloud genomes. It integrates both information on protein-coding genes and their genomic neighborhood to build a graph of gene families where each node is a gene family and each edge is a relation of genetic contiguity. The partitioning method promotes that two gene families that are consistent neighbors in the graph are more likely to belong to the same partition. It results in a Partitioned Pangenome Graph (PPG) made of persistent, shell and cloud nodes drawing genomes on rails like a subway map to help biologists navigate the great diversity of microbial life.

Moreover, the panRGP method (Bazin et al. 2020) included in PPanGGOLiN predicts, for each genome, Regions of Genome Plasticity (RGPs) that are clusters of genes made of shell and cloud genomes in the pangenome graph. Most of them arise from Horizontal gene transfer (HGT) and correspond to Genomic Islands (GIs). RGPs from different genomes are next grouped in spots of insertion based on their conserved flanking persistent genes.

{: .source}
{: .output}
{% include links.md %}
