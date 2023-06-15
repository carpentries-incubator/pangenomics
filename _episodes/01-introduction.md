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

### A brief history of how the term pangenome was created

The concept of **Pangenomics** originated from [Tettelin *et al.*](https://www.pnas.org/doi/10.1073/pnas.0506758102]), whose goal was to develop a vaccine against Group B *Streptococcus* (GBS, or *Streptococcus agalactiae*), a human pathogen causing neonatal infections. Previous to this, reverse vaccinology had been successfully applied to *Neisseria meningitidis* using a single genome. However, in the case of GBS, two complete gap-free sequences were available when the project started. These initial genomic studies revealed significant variability in gene content among closely related GBS isolates, challenging the assumption that a single genome could represent an entire species. Consequently, the collaborative team decided to sequence six additional GBS genomes, representing the major disease-causing serotypes. The comparison of these genomes confirmed the presence of diverse regions, including differential pathogenicity islands, and revealed that the shared core set of genes accounted for only about 80% of an individual genome. The existence of broad genomic diversity prompted the question of how many genomes needed to be sequenced to identify all possible genes harbored by GBS as a whole. Motivated by the goal of identifying vaccine candidates, the collaborators engaged in active discussions, scientific drafts, and the development of a mathematical model to determine the optimal number of sequenced strains. And, this is how these pioneering authors introduced the revolutionary concept of the pangenome in 2005. 


The term "pangenome" is a fusion of the Greek words *pan*, which means 'whole', and *genome*, referring to the complete set of genes in an organism. By definition, a **pangenome** represents the entirety of genes present in a group of organisms belonging to the same species. Notably, the pangenome concept extends beyond bacteria and can be applied to any taxa of interest, including humans, animals, plants, fungi, archea, or viruses. 

> ## Analogy 
> Do you feel confused about what represents a pangenome? Look at this analogy!
>  
> > ## Solution
> > Imagine you're on a mission to open the finest pizza restaurant in New York City, aiming to offer your customers a wide variety of pizzas from around the world. To achieve this, you set out to gather all the pizza recipes ever created, including pepperoni, Hawaiian, vegetarian, and more. As you examine these recipes, you begin to notice that certain ingredients appear in multiple pizzas, while others are unique to specific recipes. This analogy helps us grasp the concept of pangenomics.

In this analogy, your collection of pizza recipes represents the **pangenome**, which encompasses the entire diversity of pizzas. Each pizza recipe corresponds to a **genome**, while each ingredient in the recipes represents a **gene**. The common ingredients shared among the recipes, such as flour, tomato sauce, and mozzarella cheese, can be thought of as **families of orthologous genes**. Within these common ingredients, there may be variations in brands or preparation styles, reflecting the **genetic diversity** within the orthologous gene families.

As you continue to add more recipes to your collection, you gain a better understanding of the vast diversity in pizza-making techniques. This enables you to fulfill your objective of offering your customers the most comprehensive and diverse pizza menu.
> {: .solution}  
{: .discussion}
 

### The components and classification of pangenomes

The pangenome consists of two primary components: **core genome** and **accessory genome**. The core genome comprises genes that are present in all the genomes being compared, while the accessory genome consists of genes that are not shared by all genomes. Within the accessory genome, we can further distinguish two partitions, the *shell genome*, which encompasses the genes found in the majority of genomes, and the *cloud genome*, which comprises genes present in only a minority of genomes. It is worth mentioning that, the specific definition and percentages of these partitions may vary across different pangenome analysis software and among researchers. Additional terms such as [*persistent genome* and *soft-core genome*](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007732) are also commonly used in the field.

The concept of pangenome encompasses two types: the **open pangenome** and the **closed pangenome**. An open pangenome refers to a scenario where the addition of new genomes to a species leads to a significant increase in the total number of genes. This indicates a high level of genomic diversity and the potential acquisition of new genetic traits with each newly aggregated genome. In contrast, a closed pangenome occurs when the incorporation of new genomes does not contribute significantly to the overall gene count. In a closed pangenome, the gene pool remains relatively stable and limited, indicating a lower degree of genomic diversity within the species.


<a href="{{ page.root }}/fig/01-01-01.png">
   <img src="{{ page.root }}/fig/01-01-01.png" alt=" Venn diagram of a) a closed pangenome and b) an open pangenome, comparing the sizes of their core and accessory genomes. c) Graphic depicting the differences between closed and open pangenomes regarding their size, total genes in pangenome, and the number of sequenced genomes." />
  </a>


To understand these concepts better, we can visualize the pangenome as a matrix representing the presence (1) or absence (0) of orthologous gene families. The columns represent the gene families, while the rows represent the genomes added to the pangenome. In an open pangenome, the number of columns increases significantly as new genomes are added. Conversely, in a closed pangenome, the number of columns remains relatively unchanged as the number of genomes increases. This suggests that the species maintains a consistent set of orthologous gene families. Consequently, the size of the core genome in a closed pangenome closely matches that of a single complete genome of the species. In contrast, the core size of an open pangenome is relatively smaller compared to the size of an individual genome.

In summary, the terms "open pangenome" and "closed pangenome" describe the dynamic nature of gene content in a species, with the former signifying an expanding gene pool and the latter representing a more stable gene repertoire  


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





