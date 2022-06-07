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
- "PPanGGOLiN is a sotfware to create and manipulate prokaryotic pangenomes."
- "PPanGGOLiN integrates protein-coding genes and their genomic neighborhood to build a graph."
- "panRGP method predicts Regions of Genomic Plasticity which are grouped into insertion sites based on their conserved persistent flanking genes."
- "PPanGGOLiN is designed to scale up to tens of thousand of geneomes, including whole genomes, metagenomes and single-cell annotated genomes."
---


![Figure 1. PPanGGOLiN logo](../fig/logo_ppanggolin.png)

PPanGGOLiN
===============================================
**Partitioned PanGenome Graph Of Linked Neighbors**

PPanGGOLiN is a software to create and manipulate prokaryotic pangenomes. It partitionates a pagenome into persistent-, shell- and, cloud-gene families through a graphical model and a statistical approach rather than using fixed thresholds. Unlike other methods, PPanGGOLiN integrates both information about protein-coding genes and their genomic neighborhood to build a graph of gene families where each node is a gene family and the edges represent a relation of genetic contiguity. Therefore, two gene families that are consistent neighbors in the graph are more likely to belong to the same partition, yielding a partitioned pangenome graph (PPG) made up of persistent, shell, and cloud nodes. The resulting plot looks like a subway map, where the rails represent the genomes.

PPanGGOLiN analysis can start from genomic DNA sequences (.fasta) or annotated genomes (.gbk) of whole genomes, Metagenomic Assembled Genomes (MAG), and Single-cell Amplified Genomes (SAG), useful for large-scale environmental studies, including the non-cultivable species pangenome.  It is designed to scale up to tens of thousands of genomes. In addition, PPanGGOLiN includes the panRGP method (Bazin et al. 2020) that predicts Regions of Genomic Plasticity (RGP) for each genome. RGPs are groups of genes made of shell and cloud genomes in the pangenome chart, most of which arise from horizontal gene transfer and correspond to genomic islands. RGPs from different genomes are then grouped into insertion sites based on their conserved persistent flanking genes.

Step by step pangenome analysis with PPanGGOLiN
===============================================

Before start using PPanGGOLiN, activate the Pangenomics environment 

~~~
conda activate Pangenomics
~~~
{: .source}

~~~
(Pangenomics) betterlab@betterlabub:~$
~~~
{: .output}


Step 1
===============================================
**Identify and explore the genome files (.gbk)**

~~~
cd ~/GenomeMining/datos/gbk
ls *.gbk
~~~
{: .source}

~~~
Streptococcus_agalactiae_18RS21.gbk  Streptococcus_agalactiae_CJB111.gbk
Streptococcus_agalactiae_515.gbk     Streptococcus_agalactiae_COH1.gbk
Streptococcus_agalactiae_A909.gbk    Streptococcus_agalactiae_H36B.gbk
~~~
{: .output}


Step 2
===============================================
**Obtain a tsv-separated file with the genomes information**

Each line of this file represent one organism, first column contains a unique organism name and the second column contains the path to the associate gbk file.

~~~
ls *.gbk | cut -d'.' -f1|while read line; do echo $line$'\t/home/betterlab/GenomeMining/datos/gbk/'$line.gbk >> organisms.gbk.list; done
head organism.gbk.list
~~~
{: .source}

~~~
Streptococcus_agalactiae_18RS21	~/GenomeMining/datos/gbk/Streptococcus_agalactiae_18RS21.gbk
Streptococcus_agalactiae_515	~/GenomeMining/datos/gbk/Streptococcus_agalactiae_515.gbk
Streptococcus_agalactiae_A909	~/GenomeMining/datos/gbk/Streptococcus_agalactiae_A909.gbk
Streptococcus_agalactiae_CJB111	~/GenomeMining/datos/gbk/Streptococcus_agalactiae_CJB111.gbk
Streptococcus_agalactiae_COH1	~/GenomeMining/datos/gbk/Streptococcus_agalactiae_COH1.gbk
Streptococcus_agalactiae_H36B	~/GenomeMining/datos/gbk/Streptococcus_agalactiae_H36B.gbk
~~~
{: .output}


Step 3
===============================================
**Create a work directory for PPanGGOLiN analysis**

~~~
cd 
cd ~/Pangenomics
mkdir ppanggolin
ls -la
~~~
{: .source}

~~~
drwxrwxr-x  3 betterlab betterlab 4096 Jun  6 15:10 .
drwxrwxr-x 16 betterlab betterlab 4096 Jun  6 15:09 ..
drwxrwxr-x  2 betterlab betterlab 4096 Jun  6 15:10 ppanggolin
~~~
{: .output}


Step 4
===============================================
**Copy the organisms.gbk.list file into the work directory**

~~~
cd ppanggolin
cp ~/GenomeMining/datos/gbk/organisms.gbk.list .
ls
~~~
{: .source}

~~~
organisms.gbk.list
~~~
{: .output}


Step 5
===============================================
**Genome annotation**

Using the organisms list, annotation of genomes is made with the 'annotate' module of PPanGGOLiN

~~~
ppanggolin annotate --anno organisms.gbk.list --output pangenome
~~~
{: .source}

~~~
2022-06-06 15:55:00 main.py:l180 INFO   Command: /home/betterlab/.conda/envs/Pangenomics/bin/ppanggolin annotate --anno organisms.gbk.list --output pangenome
2022-06-06 15:55:00 main.py:l181 INFO   PPanGGOLiN version: 1.1.136
2022-06-06 15:55:00 annotate.py:l338 INFO       Reading organisms.gbk.list the list of organism files ...
100%|███████████████████████████████████████████████████████████████████████████████████████████| 6/6 [00:00<00:00,  6.53file/s]
2022-06-06 15:55:01 writeBinaries.py:l481 INFO  Writing genome annotations...
100%|████████████████████████████████████████████████████████████████████████████████████████| 6/6 [00:00<00:00, 133.20genome/s]
2022-06-06 15:55:02 writeBinaries.py:l494 INFO  writing the protein coding gene dna sequences
100%|███████████████████████████████████████████████████████████████████████████████| 13633/13633 [00:00<00:00, 129708.48gene/s]
2022-06-06 15:55:02 writeBinaries.py:l530 INFO  Done writing the pangenome. It is in file : pangenome/pangenome.h5
~~~
{: .output}

Now a new directory was created
~~~
ls
~~~
{: .source}

~~~
organisms.gbk.list  pangenome
~~~
{: .output}

Move into the pangenome/ directory and explore it. 
~~~
cd pangenome/
ls
~~~
{: .source}

~~~
pangenome.h5
~~~
{: .output}

The pangenome.h5 file will be used as input and output for all subsequent analysis 

Step 6
===============================================
**Gene clustering**

~~~
ppanggolin cluster --pangenome pangenome.h5 --cpu 8
~~~
{: .source}

~~~
2022-06-06 16:01:06 main.py:l180 INFO   Command: /home/betterlab/.conda/envs/Pangenomics/bin/ppanggolin cluster --pangenome pangenome.h5 --cpu 8
2022-06-06 16:01:06 main.py:l181 INFO   PPanGGOLiN version: 1.1.136
2022-06-06 16:01:06 readBinaries.py:l37 INFO    Getting the current pangenome's status
2022-06-06 16:01:06 readBinaries.py:l78 INFO    Extracting and writing CDS sequences from a .h5 pangenome file to a fasta file...
100%|███████████████████████████████████████████████████████████████████████████████| 13633/13633 [00:00<00:00, 113867.25gene/s]
2022-06-06 16:01:06 cluster.py:l187 INFO        Clustering all of the genes sequences...
2022-06-06 16:01:06 cluster.py:l45 INFO Creating sequence database...
2022-06-06 16:01:06 cluster.py:l54 INFO Clustering sequences...
2022-06-06 16:01:07 cluster.py:l56 INFO Extracting cluster representatives...
2022-06-06 16:01:07 cluster.py:l68 INFO Writing gene to family informations
2022-06-06 16:01:07 cluster.py:l195 INFO        Associating fragments to their original gene family...
2022-06-06 16:01:07 cluster.py:l30 INFO Aligning cluster representatives...
2022-06-06 16:01:09 cluster.py:l35 INFO Extracting alignments...
2022-06-06 16:01:09 cluster.py:l97 INFO Starting with 4565 families
2022-06-06 16:01:09 cluster.py:l126 INFO        Ending with 2894 gene families
2022-06-06 16:01:09 cluster.py:l148 INFO        Adding protein sequences to the gene families
2022-06-06 16:01:09 cluster.py:l130 INFO        Adding 13633 genes to the gene families
100%|███████████████████████████████████████████████████████████████████████████████| 13633/13633 [00:00<00:00, 611201.39gene/s]
2022-06-06 16:01:09 cluster.py:l286 INFO        Done with the clustering
2022-06-06 16:01:09 writeBinaries.py:l499 INFO  Writing gene families and gene associations...
100%|██████████████████████████████████████████████████████████████████████████| 2894/2894 [00:00<00:00, 670180.86gene family/s]
2022-06-06 16:01:09 writeBinaries.py:l501 INFO  Writing gene families information...
100%|██████████████████████████████████████████████████████████████████████████| 2894/2894 [00:00<00:00, 375787.62gene family/s]
2022-06-06 16:01:09 writeBinaries.py:l421 INFO  Updating annotations with fragment information
100%|███████████████████████████████████████████████████████████████████████████████| 14288/14288 [00:00<00:00, 459857.85gene/s]
2022-06-06 16:01:09 writeBinaries.py:l530 INFO  Done writing the pangenome. It is in file : pangenome.h5
~~~
{: .output}

The results are saved in the pangenome.h5 file given as input

Step 7
===============================================
**Build the pangenome graph**

~~~
ppanggolin graph --pangenome pangenome.h5 --cpu 8 
~~~
{: .source}

~~~
2022-06-06 16:01:49 main.py:l180 INFO   Command: /home/betterlab/.conda/envs/Pangenomics/bin/ppanggolin graph --pangenome pangenome.h5 --cpu 8
2022-06-06 16:01:49 main.py:l181 INFO   PPanGGOLiN version: 1.1.136
2022-06-06 16:01:49 readBinaries.py:l37 INFO    Getting the current pangenome's status
2022-06-06 16:01:49 readBinaries.py:l294 INFO   Reading pangenome annotations...
100%|███████████████████████████████████████████████████████████████████████████████| 14288/14288 [00:00<00:00, 350954.07gene/s]
100%|███████████████████████████████████████████████████████████████████████████████████████| 6/6 [00:00<00:00, 35.74organism/s]
2022-06-06 16:01:49 readBinaries.py:l307 INFO   Reading pangenome gene families...
100%|███████████████████████████████████████████████████████████████████████████████| 13633/13633 [00:00<00:00, 311021.25gene/s]
100%|██████████████████████████████████████████████████████████████████████████| 2894/2894 [00:00<00:00, 199554.73gene family/s]
2022-06-06 16:01:49 makeGraph.py:l56 INFO       Computing the neighbors graph...
Processing Streptococcus_agalactiae_H36B: 100%|████████████████████████████████████████████| 6/6 [00:00<00:00, 316.89organism/s]
2022-06-06 16:01:49 makeGraph.py:l74 INFO       Done making the neighbors graph.
2022-06-06 16:01:49 writeBinaries.py:l508 INFO  Writing the edges...
100%|█████████████████████████████████████████████████████████████████████████████████| 3188/3188 [00:00<00:00, 719746.00edge/s]
2022-06-06 16:01:49 writeBinaries.py:l530 INFO  Done writing the pangenome. It is in file : pangenome.h5
~~~
{: .output}

The results are saved in the pangenome.h5 file given as input

Step 8
===============================================
**Pangenome partition**

This is the step that will assign gene families to the 'persistent', 'shell', or 'cloud' partitions.
The one parameter that might be of importance is the '-K', or '--nb_of_partitions' parameter. This will define the number of classes used to partition the pangenome. This may be of use if you expect to have well-defined subpopulations in your pangenome and you know exactly how many. If not, that number if detected automatically through an ICL criterion. The idea is that the most present partition will be 'persistent', the least present will be 'cloud', and all the others will be 'shell'. The number of partitions corresponding to the shell will be the number of expected subpopulations in your pangenome. (So if you expect 5 subpopulations, you could use -K 7).

In most cases, you should let the statistical criterion used by PPanGGOLiN find the optimal number of partitions for you.

~~~
ppanggolin partition --pangenome pangenome.h5 --cpu 8
~~~
{: .source}

~~~
2022-06-07 11:28:56 main.py:l180 INFO   Command: /home/betterlab/.conda/envs/Pangenomics/bin/ppanggolin partition --pangenome pangenome.h5 --cpu 8
2022-06-07 11:28:56 main.py:l181 INFO   PPanGGOLiN version: 1.1.136
2022-06-07 11:28:56 readBinaries.py:l37 INFO    Getting the current pangenome's status
2022-06-07 11:28:56 readBinaries.py:l294 INFO   Reading pangenome annotations...
100%|███████████████████████████████████████████████████████████████████████| 14288/14288 [00:00<00:00, 340902.17gene/s]
100%|███████████████████████████████████████████████████████████████████████████████| 6/6 [00:00<00:00, 35.44organism/s]
2022-06-07 11:28:56 readBinaries.py:l307 INFO   Reading pangenome gene families...
100%|███████████████████████████████████████████████████████████████████████| 13633/13633 [00:00<00:00, 305842.61gene/s]
100%|██████████████████████████████████████████████████████████████████| 2894/2894 [00:00<00:00, 200557.07gene family/s]
2022-06-07 11:28:56 readBinaries.py:l314 INFO   Reading the neighbors graph edges...
100%|███████████████████████████████████████████████████████████| 11958/11958 [00:00<00:00, 289813.92contig adjacency/s]
2022-06-07 11:28:56 partition.py:l343 WARNING   The number of selected organisms is too low (6 organisms used) to robustly partition the graph
2022-06-07 11:28:56 partition.py:l356 INFO      Estimating the optimal number of partitions...
100%|███████████████████████████████████████████████████████| 19/19 [00:00<00:00, 59.88Number of number of partitions/s]
2022-06-07 11:28:57 partition.py:l358 INFO      The number of partitions has been evaluated at 3
2022-06-07 11:28:57 partition.py:l376 INFO      Partitioning...
2022-06-07 11:28:57 partition.py:l436 INFO      Partitionned 6 genomes in 0.07 seconds.
2022-06-07 11:28:57 writeBinaries.py:l408 INFO  Updating gene families with partition information
100%|██████████████████████████████████████████████████████████████████| 2894/2894 [00:00<00:00, 195908.84gene family/s]
2022-06-07 11:28:57 writeBinaries.py:l530 INFO  Done writing the pangenome. It is in file : pangenome.h5
~~~
{: .output}

All the results will be added to the given 'pangenome.h5' input file.

Step 9
===============================================
**Predict the regions of genome plasticity with RGP module**

~~~
ppanggolin rgp --pangenome pangenome.h5 --cpu 8
~~~
{: .source}

~~~
2021-12-06 11:28:59 main.py:l180 INFO	Command: /home/betterlab/.conda/envs/Pangenomics/bin/ppanggolin rgp --pangenome pangenome.h5 --cpu 8
2021-12-06 11:28:59 main.py:l181 INFO	PPanGGOLiN version: 1.1.136
2021-12-06 11:28:59 readBinaries.py:l37 INFO	Getting the current pangenome's status
2021-12-06 11:28:59 readBinaries.py:l294 INFO	Reading pangenome annotations...
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 86820/86820 [00:00<00:00, 389209.92gene/s]
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 20/20 [00:00<00:00, 22.12organism/s]
2021-12-06 11:29:00 readBinaries.py:l307 INFO	Reading pangenome gene families...
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 85940/85940 [00:00<00:00, 258946.34gene/s]
100%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 4945/4945 [00:00<00:00, 177141.87gene family/s]
2021-12-06 11:29:01 genomicIsland.py:l197 INFO	Detecting multigenic families...
2021-12-06 11:29:01 pangenome.py:l311 INFO	20 gene families are defined as being multigenic. (duplicated in more than 0.05 of the genomes)
2021-12-06 11:29:01 genomicIsland.py:l199 INFO	Compute Regions of Genomic Plasticity ...
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 20/20 [00:00<00:00, 183.34genomes/s]
2021-12-06 11:29:01 genomicIsland.py:l204 INFO	Predicted 130 RGP
2021-12-06 11:29:01 writeBinaries.py:l517 INFO	Writing Regions of Genomic Plasticity...
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 130/130 [00:00<00:00, 430354.79region/s]
2021-12-06 11:29:01 writeBinaries.py:l530 INFO	Done writing the pangenome. It is in file : pangenome.h5
~~~
{: .output}

You also can obtain a list of the plastic regions (RGPs) for each genome by using the module write 

~~~
ppanggolin write -p pangenome.h5 --regions --output rgp
~~~
{: .source}

~~~
2021-12-06 11:30:05 main.py:l180 INFO	Command: /home/betterlab/.conda/envs/Pangenomics/bin/ppanggolin write -p pangenome.h5 --regions --output rgp
2021-12-06 11:30:05 main.py:l181 INFO	PPanGGOLiN version: 1.1.136
2021-12-06 11:30:05 readBinaries.py:l37 INFO	Getting the current pangenome's status
2021-12-06 11:30:05 readBinaries.py:l294 INFO	Reading pangenome annotations...
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 86820/86820 [00:00<00:00, 382028.80gene/s]
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 20/20 [00:00<00:00, 21.61organism/s]
2021-12-06 11:30:06 readBinaries.py:l307 INFO	Reading pangenome gene families...
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 85940/85940 [00:00<00:00, 254608.91gene/s]
100%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 4945/4945 [00:00<00:00, 175441.19gene family/s]
2021-12-06 11:30:07 readBinaries.py:l320 INFO	Reading the RGP...
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 1237/1237 [00:00<00:00, 436877.24gene/s]
~~~
{: .output}

Explore the rgp results

~~~
cd rgp/
ls
~~~
{: .source}

~~~
plastic_regions.tsv
~~~
{: .output}

~~~
head plastic_regions.tsv
~~~
{: .source}

~~~
region	organism	contig	start	stop	genes	contigBorder	wholeContig
MtbN0004_contig_15_RGP_0	Mtb_N0004_L3	MtbN0004_contig_15	83651	93061	17	True	False
MtbN0004_contig_16_RGP_0	Mtb_N0004_L3	MtbN0004_contig_16	33871	50925	23	False	False
MtbN0004_contig_3_RGP_0         Mtb_N0004_L3	MtbN0004_contig_3	178698	182189	5	True	False
MtbN0004_contig_58_RGP_0	Mtb_N0004_L3	MtbN0004_contig_58	8243	18251	16	False	False
MtbN0031_contig_10_RGP_0	Mtb_N0031_L2	MtbN0031_contig_10	79495	97528	24	False	False
MtbN0031_contig_14_RGP_0	Mtb_N0031_L2	MtbN0031_contig_14	81684	86828	6	False	False
MtbN0031_contig_29_RGP_0	Mtb_N0031_L2	MtbN0031_contig_29	23141	28818	8	False	False
MtbN0031_contig_55_RGP_0	Mtb_N0031_L2	MtbN0031_contig_55	1457	11465	16	False	False
MtbN0052_contig_17_RGP_0	Mtb_N0052_L2	MtbN0052_contig_17	43797	61393	23	False	False
~~~
{: .output}

Return to the working directory
~~~
cd ..
~~~
{: .source}

Step 8
===============================================
**Compute the spots of insertion**

~~~
ppanggolin spot --pangenome pangenome.h5 --cpu 8
~~~
{: .source}

~~~
2021-12-06 11:46:59 main.py:l180 INFO	Command: /home/betterlab/.conda/envs/Pangenomics/bin/ppanggolin spot --pangenome pangenome.h5 --cpu 8
2021-12-06 11:46:59 main.py:l181 INFO	PPanGGOLiN version: 1.1.136
2021-12-06 11:46:59 readBinaries.py:l37 INFO	Getting the current pangenome's status
2021-12-06 11:46:59 readBinaries.py:l294 INFO	Reading pangenome annotations...
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 86820/86820 [00:00<00:00, 384355.02gene/s]
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 20/20 [00:00<00:00, 21.97organism/s]
2021-12-06 11:47:00 readBinaries.py:l307 INFO	Reading pangenome gene families...
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 85940/85940 [00:00<00:00, 255481.24gene/s]
100%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 4945/4945 [00:00<00:00, 177932.10gene family/s]
2021-12-06 11:47:01 readBinaries.py:l320 INFO	Reading the RGP...
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 1237/1237 [00:00<00:00, 438501.86gene/s]
2021-12-06 11:47:01 spot.py:l129 INFO	Detecting multigenic families...
2021-12-06 11:47:01 pangenome.py:l311 INFO	20 gene families are defined as being multigenic. (duplicated in more than 0.05 of the genomes)
2021-12-06 11:47:01 spot.py:l132 INFO	Detecting hotspots in the pangenome...
2021-12-06 11:47:01 spot.py:l82 INFO	66 RGPs were not used as they are on a contig border (or have less than 3 persistent gene families until the contig border)
2021-12-06 11:47:01 spot.py:l83 INFO	64 RGPs are being used to predict spots of insertion
2021-12-06 11:47:01 spot.py:l85 INFO	14 number of different pairs of flanking gene families
2021-12-06 11:47:01 spot.py:l140 INFO	11 spots were detected
2021-12-06 11:47:01 writeBinaries.py:l522 INFO	Writing Spots of Insertion...
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 11/11 [00:00<00:00, 228402.69spot/s]
2021-12-06 11:47:01 writeBinaries.py:l530 INFO	Done writing the pangenome. It is in file : pangenome.h5
~~~
{: .output}


**THIS VERSION DO NOT ALLOW 'MODULE' NOR 'CONTEXT' ANALYSIS**


Step 9
===============================================
**Compute the pangenome results**

PPanGGOLiN provides multiple outputs to describe a pangenome. In this section the different outputs will be described

**1. DRAW**

**1.1 U-shaped plot**

A U-shaped plot is a figure presenting the number of families (y axis) per number of organisms (x axis). It is a .html file that can be opened with any browser and with which you can interact, zoom, move around, mouseover to see numbers in more detail, and you can save what you are seeing as a .png image file.

~~~
ppanggolin draw --pangenome pangenome.h5 --ucurve --output draw_ucurve
~~~
{: .source}

~~~
2021-12-06 12:10:30 main.py:l180 INFO	Command: /home/betterlab/.conda/envs/Pangenomics/bin/ppanggolin draw --pangenome pangenome.h5 --ucurve --output draw_ucurve
2021-12-06 12:10:30 main.py:l181 INFO	PPanGGOLiN version: 1.1.136
2021-12-06 12:10:30 readBinaries.py:l37 INFO	Getting the current pangenome's status
2021-12-06 12:10:30 readBinaries.py:l294 INFO	Reading pangenome annotations...
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 86820/86820 [00:00<00:00, 383434.30gene/s]
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 20/20 [00:00<00:00, 21.97organism/s]
2021-12-06 12:10:31 readBinaries.py:l307 INFO	Reading pangenome gene families...
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 85940/85940 [00:00<00:00, 255888.23gene/s]
100%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 4945/4945 [00:00<00:00, 178552.47gene family/s]
2021-12-06 12:10:31 readBinaries.py:l314 INFO	Reading the neighbors graph edges...
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 83033/83033 [00:00<00:00, 192821.36contig adjacency/s]
2021-12-06 12:10:32 ucurve.py:l13 INFO	Drawing the U-shaped curve...
2021-12-06 12:10:32 ucurve.py:l60 INFO	Done drawing the U-shaped curve : 'draw_ucurve/Ushaped_plot.html'
~~~
{: .output}

~~~
cd draw_ucurve/
ls
~~~
{: .source}

~~~
Ushaped_plot.html
~~~
{: .output}

**1.2 Tile plot**

A tile plot is a heatmap representing the gene families (y axis) in the organisms (x axis) making up your pangenome. The tiles on the graph will be colored if the gene family is present in an organism and uncolored if absent. The gene families are ordered by partition, and the genomes are ordered by a hierarchical clustering based on their shared gene families (basically two genomes that are close together in terms of gene family composition will be close together on the figure).

This plot is quite helpful to observe potential structures in your pangenome, and can also help you to identify eventual outliers. You can interact with it, and mousing over a tile in the plot will indicate to you which is the gene identifier(s), the gene family and the organism that corresponds to the tile.

If you build your pangenome using the 'workflow' subcommand and you have more than 500 organisms, only the 'shell' and the 'persistent' partitions will be drawn, leaving out the 'cloud' as the figure tends to be too heavy for a browser to open it otherwise.

~~~
ppanggolin draw --pangenome pangenome.h5 --tile_plot --output draw_tile
~~~
{: .source}

~~~
2021-12-06 12:09:32 main.py:l180 INFO	Command: /home/betterlab/.conda/envs/Pangenomics/bin/ppanggolin draw --pangenome pangenome.h5 --tile_plot --output draw_tile
2021-12-06 12:09:32 main.py:l181 INFO	PPanGGOLiN version: 1.1.136
2021-12-06 12:09:32 readBinaries.py:l37 INFO	Getting the current pangenome's status
2021-12-06 12:09:32 readBinaries.py:l294 INFO	Reading pangenome annotations...
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 86820/86820 [00:00<00:00, 387154.17gene/s]
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 20/20 [00:00<00:00, 21.51organism/s]
2021-12-06 12:09:33 readBinaries.py:l307 INFO	Reading pangenome gene families...
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 85940/85940 [00:00<00:00, 258514.00gene/s]
100%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 4945/4945 [00:00<00:00, 175470.88gene family/s]
2021-12-06 12:09:33 readBinaries.py:l314 INFO	Reading the neighbors graph edges...
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 83033/83033 [00:00<00:00, 196168.75contig adjacency/s]
2021-12-06 12:09:34 tile_plot.py:l26 INFO	Drawing the tile plot...
2021-12-06 12:09:34 tile_plot.py:l42 INFO	start with matrice
2021-12-06 12:09:35 tile_plot.py:l57 INFO	done with making the dendrogram to order the organisms on the plot
2021-12-06 12:09:35 tile_plot.py:l92 INFO	Getting the gene name(s) and the number for each tile of the plot ...
2021-12-06 12:09:35 tile_plot.py:l101 INFO	Done extracting names and numbers. Making the heatmap ...
2021-12-06 12:09:35 tile_plot.py:l157 INFO	Drawing the figure itself...
2021-12-06 12:09:37 tile_plot.py:l159 INFO	Done with the tile plot : 'draw_tile/tile_plot.html' 
~~~
{: .output}

If you do not want the 'cloud' gene families as it is a lot of data and can be hard to open with a browser sometimes, you can use the following option:

~~~
ppanggolin draw --pangenome pangenome.h5 --tile_plot --nocloud --output draw_tile_nocloud
~~~
{: .source}

~~~
2021-12-06 12:14:19 main.py:l180 INFO	Command: /home/betterlab/.conda/envs/Pangenomics/bin/ppanggolin draw --pangenome pangenome.h5 --tile_plot --nocloud --output draw_tile_nocloud
2021-12-06 12:14:19 main.py:l181 INFO	PPanGGOLiN version: 1.1.136
2021-12-06 12:14:19 readBinaries.py:l37 INFO	Getting the current pangenome's status
2021-12-06 12:14:19 readBinaries.py:l294 INFO	Reading pangenome annotations...
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 86820/86820 [00:00<00:00, 380327.86gene/s]
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 20/20 [00:00<00:00, 21.65organism/s]
2021-12-06 12:14:21 readBinaries.py:l307 INFO	Reading pangenome gene families...
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 85940/85940 [00:00<00:00, 255772.75gene/s]
100%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 4945/4945 [00:00<00:00, 178563.23gene family/s]
2021-12-06 12:14:21 readBinaries.py:l314 INFO	Reading the neighbors graph edges...
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 83033/83033 [00:00<00:00, 191715.85contig adjacency/s]
2021-12-06 12:14:21 tile_plot.py:l26 INFO	Drawing the tile plot...
2021-12-06 12:14:21 tile_plot.py:l42 INFO	start with matrice
2021-12-06 12:14:21 tile_plot.py:l57 INFO	done with making the dendrogram to order the organisms on the plot
2021-12-06 12:14:21 tile_plot.py:l92 INFO	Getting the gene name(s) and the number for each tile of the plot ...
2021-12-06 12:14:21 tile_plot.py:l101 INFO	Done extracting names and numbers. Making the heatmap ...
2021-12-06 12:14:22 tile_plot.py:l157 INFO	Drawing the figure itself...
2021-12-06 12:14:24 tile_plot.py:l159 INFO	Done with the tile plot : 'draw_tile_nocloud/tile_plot.html' 
~~~
{: .output}

**1.3 Spots plot**

**THIS VERSION DO NOT ALLOW 'MODULE' NOR 'CONTEXT' ANALYSIS**



{% include links.md %}
