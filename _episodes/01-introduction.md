---
title: "Introduction to Pangenomics"
teaching: 10
exercises: 10
questions:
- "How can I represent all the gene content of a species?" 
objectives:
- "Discuss the importance of the pangenome concept."
- "Define the different parts of the pangenome."
- "Define the types of pangenomes."
- "Know the data that we will use to build a pangenome."

keypoints:
- "A pangenome is the entire set of genes from all genomes in a group of interest (typically a species)."
- "The pangenome is composed of a core genome and an accesory genome."
- "The accessory genome is usually subdivided in shell genome and cloud genome."
- "An *open pangenome* occurs when the size of the pangenome increases considerably with every added genome."
- "A *closed pangenome* occurs when only a few gene families are incorporated to the pangenome when a new genome is added."
- "A pangenome can be obtained by comparing the complete genome sequences of all clade members."
---

## Unveiling the Genomic Complexity: The Pangenome Paradigm

The concept of **Pangenomics** originated from [Tettelin *et al.*](https://www.pnas.org/doi/10.1073/pnas.0506758102]), whose goal was to develop a vaccine against Group B *Streptococcus* (GBS, or *Streptococcus agalactiae*), a human pathogen causing neonatal infections. Previous to this, reverse vaccinology had been successfully applied to *Neisseria meningitidis* using a single genome. However, in the case of GBS, two complete gap-free sequences were available when the project started. These initial genomic studies revealed significant variability in gene content among closely related GBS isolates, challenging the notion that a single genome could represent the entire species. Consequently, the collaborative team decided to sequence six additional GBS genomes, representing the major disease-causing serotypes.   

Contrary to what was thought in the past, one genome is not enough to represent all of the genomic content of a species. This was observed 
by [Tettelin *et al.*](https://www.pnas.org/doi/10.1073/pnas.0506758102]),  who were working on developing a vaccine against the bacteria 
*Streptococcus agalactiae* (a human pathogen that causes neonatal infections) and realized that one genome was not enough to choose an appropriate 
gene to be the target of the vaccine.  

In 2005, these authors came up with the concept of **pangenome**, which derives from the Greek *pan*, meaning ‘whole’ or ‘everything’ and *genome*, to 
refer to the collection of all the genes present in a species. This concept can be applied to any group of interest, not only to the taxonomic 
level of species. And although it started from the study of bacteria, it has also been used in the research of eukaryotes, archaea, and viruses.  

The pangenome can be divided into two principal components: the **core genome**, which contains genes present in all the genomes compared, and the 
**dispensable** or **accessory genome**, that comprises the genes that are not shared by all the genomes. The accesory can be subdivided into the 
**shell genome**, the proportion of genes present in the majority of the genomes, and **cloud genome**, comprised of genes in the minority of the 
genomes. This is a somewhat ambiguous definition of the partitions, in the practice, the percentage of genomes that define each partition can be set 
differently by each software of pangenome analysis and by each researcher. And other terms such as
[persistent genome and soft-core genome](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007732) are also commonly 
used in the field.

One of the most groundbreaking ideas that pangenomics brought to the study of microbiology is that a species can be unbounded. When each newly 
included genome increases considerably the number of total genes, the species has an **open pangenome**. Whereas a **closed pangenome** results 
when only a few gene families are incorporated into the pangenome when a new member is added.

<a href="{{ page.root }}/fig/01-01-01.png">
   <img src="{{ page.root }}/fig/01-01-01.png" alt=" Venn diagram of a) a closed pangenome and b) an open pangenome, comparing the sizes of their core and accessory genomes. c) Graphic depicting the differences between closed and open pangenomes regarding their size, total genes in pangenome, and the number of sequenced genomes." />
  </a>

> ## Discussion: Open or closed?
>  The amount of gene transfer, the degree of interaction with other species in the environment, the number of niches inhabited and the lifestyle 
>  of the species can play a role in the size of the pangenome. 
>  
>  Between a human lung pathogen and a soil bacterium, which one do you think has a closed pangenome and which one has an open one?
>  
> > ## Solution
> > Specialized pathogens tend to have smaller genomes because they live in an environment with more or less the same challenges, 
> > so they have less variability in their gene repertoire, hence a closed pangenome.  
> > On the other hand, free-living species face a changing environment, so they will be more adaptable if the population has a more diverse 
> > set of tools to deal with the conditions. These tools are usually in the accessory genome, and the bigger the accessory genome is the more
> >  open the pangenome will be.  
> {: .solution}  
{: .discussion}

> ## Know more
> If you want to read more on pangenomics go to the book [The Pangenome](https://link.springer.com/book/10.1007/978-3-030-38281-0).
{: .callout}

## The *Streptococcus agalactiae* dataset

In this lesson, we will perform a standard pangenomics pipeline, which consist of the genomic annotation, clustering of gene families by identification
of orthologous sequences, and the description of the pangenome partitions and openness. For this we will work with the eight *Streptococcus agalactiae*
 strains included in the first pangenome made by *Tettelin et al., 2005*.

The genomes of strains 18RS21 and H36B are already in our `pan_workshop/data`, but the others will be downloaded and annotated in the following episodes.

**Metadata**


|Strain	| Host	| Serotype   |
|-------------------------|---------|------------|
|*S. agalactiae*  18RS21  | Human   | II       	|
|*S. agalactiae*  515 	| Human   | Ia       	|
|*S. agalactiae*  A909	| Human   | Ia       	|
|*S. agalactiae*  CJB111  | Human   | V       	|
|*S. agalactiae*  COH1	| Human   | III       	|
|*S. agalactiae*  H36B	| Human   | Ib       	|
|*S. agalactiae*  NEM316	| Human   | III     	|
|*S. agalactiae*  2603V/R 	| Human   | V      	|



> ## Prepare your genome database
> Make sure you have the `pan_workshop/` directory. If you do not have it, you can download it with the following instructions.
>
> ~~~
> $ cd ~ #Make sure you are in the home directory
> $ wget https://zenodo.org/record/7974915/files/pan_workshop.zip?download=1 #Download the `zip` file.
> $ unzip 'pan_workshop.zip?download=1' 
> $ rm 'pan_workshop.zip?download=1'
> ~~~
> {: .language-bash}
{: .checklist}

{% include links.md %}





