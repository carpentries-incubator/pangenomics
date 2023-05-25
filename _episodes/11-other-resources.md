---
title: "Other Resources"
teaching: 10
exercises: 10
questions:
- "What can I do after I have built a pangenome?"
- "What bioinformatic tools are available for downstream analysis of pangenomes?"
objectives:
- "Learn about the downstream analysis to describe a pangenome"
keypoints:
- "Downstream analysis of pangenomes could be focused on describing the core or the accessory genome of the organism studied."
- "Examples using the information obtained in the **CORE GENOME**:" 
-   "a) Selection of a conserved gene to design a molecular test for a diagnostic tool or a vaccine."
-   "b) Reconstruction of a species phylogenetic tree by using all the core genes."
- "Examples using the information obtained in the **ACCESSORY GENOME**:"
-   "a) Describe niche specific genes among the strains compared."
-   "b) Analysis of horizontal gene transfer or genetic recombination."
-   "c) Evolutionary studies of genes (duplication, gain-loss genes, etc.)."
---
## Downstream analysis of pangenomes

At this point, we already have learnt about how to build a pangenome. We have a list of the organisms we included in the analysis and we obtained a matrix of presence-absence of genes, which roughly represents our pangenome. 

Preferable, we can also have a metadata file that includes phenotypic features related to each one of our samples. For example, year of isolation, geographic origin of the sample, host, drug-resistance profile, serotype or lineage, etc. With all of this information, there exists a great number of downstream analysis we can conduct on a pangenome. In one hand, we can focus on the *core genome* which can be useful to reconstruct phylogenetic trees, design of molecular diagnostic tools or vaccines, among others. On the other hand, we can conduct comparative genomic analysis using the *accessory genome* to describe differences between subsets of samples in our dataset thath might be related with some specific traits   
Here we want to mention some of them:

### Genetic Variation Analysis

### Functional Annotation

### Pan-GWAS Analysis

### Clustering analysis of the accesory genome

### Benchmarking for Orthologs Predictions

In addition, you might be interested in evaluate the quality of your results regarding the methodology used to identify the ortholog genes. For that specific task, you can explore the website of [**Orthology Benchmarking**](https://orthology.benchmarkservice.org/), which shows you a pre-processed benchmarking analysis of several methodologies ranked for their efficiency to identify real orthologs (recal) and also gives you the opportunity to evaluate your own methods and results. 

Find more information about how to submit your own data in:

[**How to submit predictions to Orthology Benchmarking**](https://orthology.benchmarkservice.org/proxy/doc#submit)

**Original paper:**
Brigitte Boeckmann and others, Conceptual framework and pilot study to benchmark phylogenomic databases based on reference gene trees, Briefings in Bioinformatics, Volume 12, Issue 5, September 2011, Pages 423–435, https://doi.org/10.1093/bib/bbr034


### OrthoVenn3



## Carpentries Philosophy
A good lesson should be as complete and clear that becomes easy to teach by any instructor. 
Carpentries lessons are developed for the community, and now you are part of us. 
This lesson is being developed and we are sure that you can colaborate and help us improve it.
{% include links.md %}
