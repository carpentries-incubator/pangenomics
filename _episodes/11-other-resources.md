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

At this stage, we have acquired a comprehensive understanding of pangenome construction techniques. Our analysis has yielded a gene presence-absence matrix, which serves as a rudimentary representation of our pangenome. Ideally, we can complement this data by incorporating a metadata file encompassing various phenotypic features associated with each of our samples. These features may include the year of isolation, geographic origin, host information, drug-resistance profiles, levels of pathogenicity or virulence, serotypes, lineages, and more. Armed with this extensive metadata, numerous intriguing avenues for downstream analyses become available to us. By exploring the pangenome, we can uncover invaluable insights and elucidate various aspects of microbial genetics. 

Firstly, a focused investigation of the **core genome** offers us the opportunity to reconstruct robust phylogenetic trees. By examining the conserved genetic elements shared among organisms, we can establish evolutionary relationships, elucidate common ancestry, and gain a deeper understanding of the microbial diversity within our dataset. Furthermore, these phylogenetic trees serve as a foundation for designing effective molecular diagnostic tools and formulating targeted antimicrobial drugs or vaccines to combat a broad range of strains or lineages. 

Simultaneously, we can harness the power of the **accessory genome** to conduct comparative genomic analyses. This allows us to discern dissimilarities between subsets of samples within our dataset, providing crucial insights into specific traits or characteristics. By identifying genes or genetic variations unique to certain groups of organisms, we can unreveal the genetic basis of their phenotypic diversity, such as drug-resistance profiles, pathogenicity, or virulence levels. 

### Pan-GWAS Analysis

The Genome-Wide Association Studies (GWAS) represent a strategy to investigate the relationship between genetic variations and phenotypic traits in large populations. Typically, the GWAS process involves i) scanning the entire genomes to look for variations, such as single nucleotide polymorphisms (SNPs) and insertions/deletions (indels); and ii) identifying these variable regions within the genome that are statistically linked to a trait (*p.e.* hearth disease, fruit shape,  and more). Thus, GWAS is a fundamental technique for human and agricultural investigations. However, for prokaryote organisms this strategy is just recently began to be used. 

In Pangenomics, the GWAS strategy could also be applied 




### Clustering analysis of the accesory genome

### Benchmarking for Orthologs Predictions

You might be interested in evaluate the quality of your results concerning the methodology employed for indentifying orthologous genes. To facilitate this task, you can utilize the resources provided by the [**Orthology Benchmarking**](https://orthology.benchmarkservice.org/) website. This platform offers a pre-processed benchmarking analysis of various methodologies, ranking them based on their efficacy in identifying authentic orthologs (recal). Additionally, it provides the opportunity to assess your own methods and results within this framework. 

For detailed instruction on how to submit your own data for evaluation, please refer to the following information:

[**How to submit my ortholog predictions to Orthology Benchmarking**](https://orthology.benchmarkservice.org/proxy/doc#submit)

**Original paper:**
Brigitte Boeckmann and others, Conceptual framework and pilot study to benchmark phylogenomic databases based on reference gene trees, Briefings in Bioinformatics, Volume 12, Issue 5, September 2011, Pages 423–435, https://doi.org/10.1093/bib/bbr034


### OrthoVenn3


In conclusion, the vast array of downstream analyses made possible by our pangenome dataset and accompanying metadata file holds immense promise. By exploring the core and accessory genomes, along with phenotypic correlations, we can gain valuable insights into evolutionary relationships, understand the genetic basis of specific traits, and contribute to a better knowledge about the adaptive processes bacteria have achieved in their respective lifestyle. 


## Carpentries Philosophy
A good lesson should be as complete and clear that becomes easy to teach by any instructor. 
Carpentries lessons are developed for the community, and now you are part of us. 
This lesson is being developed and we are sure that you can colaborate and help us improve it.
{% include links.md %}
