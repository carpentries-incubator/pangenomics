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

At this stage, we have acquired a comprehensive understanding of pangenome construction techniques. Our analysis has yielded a gene presence-absence matrix, which serves as a rudimentary representation of our pangenome. Ideally, we can complement this data by incorporating a metadata file encompassing various phenotypic features associated with each of our samples. These features may include the year of isolation, geographic origin, host information, drug-resistance profiles, levels of pathogenicity or virulence, serotypes, lineages, and more. Armed with this extensive metadata, a plethora of downstream analyses can be persued on the pangenome. 

Firstly, by focusing on the **core genome**, we can employ it to reconstruct phylogenetic species trees, devise molecular diagnostic tools, design vaccines, and explore other relevant avenues by using the highly conserved regions of the species studied. Simultaneously, utilizing the **accessory genome** enables us to conduct comparative genomic analyses, unrevealing discrepancies among subsets of samples within our dataset that may correlate with distinct characteristics or traits. Now, let us delve into a selection of those potential downstream analyses: 

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
