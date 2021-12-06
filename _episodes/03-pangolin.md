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

conda activate Pangenomics


Step 1
===============================================
**Create a work directory for PPanGGOLiN analysis**

cd Pangenomics/
mkdir PPanGGOLiN
ls -la

drwxrwxr-x  6 betterlab betterlab 4096 Dec  6 08:31 .
drwxr-xr-x 42 betterlab betterlab 4096 Dec  3 18:37 ..
drwxrwxr-x  3 betterlab betterlab 4096 Dec  3 16:12 Anvio
drwxrwxr-x  2 betterlab betterlab 4096 Dec  6 08:32 PPanGGOLiN
drwxrwxr-x  2 betterlab betterlab 4096 Dec  6 09:08 gbk


Step 2
===============================================
**Copy the genome files (.gbk) into the PPanGGOLiN directory**

cd gbk
cp *.gbk ../PPanGGOLiN/.
cd ../PPanGGOLiN/
ls

Mtb_N0004_L3.gbk  Mtb_N0069_L1.gbk  Mtb_N0145_L2.gbk  Mtb_N1201_L6.gbk  Mtb_N1272_L5.gbk
Mtb_N0031_L2.gbk  Mtb_N0072_L1.gbk  Mtb_N0155_L2.gbk  Mtb_N1202_L6.gbk  Mtb_N1274_L3.gbk
Mtb_N0052_L2.gbk  Mtb_N0091_L6.gbk  Mtb_N0157_L1.gbk  Mtb_N1216_L4.gbk  Mtb_N1283_L4.gbk
Mtb_N0054_L3.gbk  Mtb_N0136_L4.gbk  Mtb_N1176_L5.gbk  Mtb_N1268_L5.gbk  Mtb_N3913_L7.gbk

Step 3
===============================================
**Obtain a tsv-separated file with the genomes information**

Each line of this file represent one organism, first column contains a unique organism name and the second column contains the associate gbk file.
First, we move into the gbk/ directory which contains all the genome files. 

cd ~/Pangenomics/gbk
ls *.gbk | cut -d'.' -f1|while read line; do echo $line$'\t'$line.gbk >> organisms.gbk.list; done
head organism.gbk.list

Mtb_N0004_L3	Mtb_N0004_L3.gbk
Mtb_N0031_L2	Mtb_N0031_L2.gbk
Mtb_N0052_L2	Mtb_N0052_L2.gbk
Mtb_N0054_L3	Mtb_N0054_L3.gbk
Mtb_N0069_L1	Mtb_N0069_L1.gbk
Mtb_N0072_L1	Mtb_N0072_L1.gbk
Mtb_N0091_L6	Mtb_N0091_L6.gbk
Mtb_N0136_L4	Mtb_N0136_L4.gbk
Mtb_N0145_L2	Mtb_N0145_L2.gbk
Mtb_N0155_L2	Mtb_N0155_L2.gbk

Then, we copy this list into the PPanGGOLiN/ work directory 

cp organisms.gbk.list ../PPanGGOLiN/.

Move into the work directory and verify your organisms list

cd ../PPanGGOLiN
ls

organisms.gbk.list

Step 3
===============================================
**Genome annotation**

Using the organisms list, annotation of genomes is made with the 'annotate' module of PPanGGOLiN

ppanggolin annotate --anno organisms.gbk.list --output pangenome

2021-12-06 10:36:41 main.py:l180 INFO	Command: /home/betterlab/.conda/envs/Pangenomics/bin/ppanggolin annotate --anno organisms.gbk.list --output pangenome
2021-12-06 10:36:41 main.py:l181 INFO	PPanGGOLiN version: 1.1.136
2021-12-06 10:36:41 annotate.py:l338 INFO	Reading organisms.gbk.list the list of organism files ...
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 20/20 [00:04<00:00,  4.34file/s]
2021-12-06 10:36:46 writeBinaries.py:l481 INFO	Writing genome annotations...
100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 20/20 [00:00<00:00, 104.24genome/s]
2021-12-06 10:36:46 writeBinaries.py:l494 INFO	writing the protein coding gene dna sequences
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 85940/85940 [00:00<00:00, 100732.00gene/s]
2021-12-06 10:36:47 writeBinaries.py:l530 INFO	Done writing the pangenome. It is in file : pangenome/pangenome.h5
(Pangenomics) betterlab@betterlabub:~/Pangenomics/PPanGGOLiN$ 



ls
Mtb_N0004_L3.gbk  Mtb_N0069_L1.gbk  Mtb_N0145_L2.gbk  Mtb_N1201_L6.gbk  Mtb_N1272_L5.gbk  pangenome
Mtb_N0031_L2.gbk  Mtb_N0072_L1.gbk  Mtb_N0155_L2.gbk  Mtb_N1202_L6.gbk  Mtb_N1274_L3.gbk  organisms.gbk.list
Mtb_N0052_L2.gbk  Mtb_N0091_L6.gbk  Mtb_N0157_L1.gbk  Mtb_N1216_L4.gbk  Mtb_N1283_L4.gbk
Mtb_N0054_L3.gbk  Mtb_N0136_L4.gbk  Mtb_N1176_L5.gbk  Mtb_N1268_L5.gbk  Mtb_N3913_L7.gbk

cd pangenome/
ls

pangenome.h5


Step 4
===============================================
**Gene clustering**

ppanggolin cluster --pangenome pangenome.h5 --cpu 8

2021-12-06 10:39:37 main.py:l180 INFO	Command: /home/betterlab/.conda/envs/Pangenomics/bin/ppanggolin cluster --pangenome pangenome.h5 --cpu 8
2021-12-06 10:39:37 main.py:l181 INFO	PPanGGOLiN version: 1.1.136
2021-12-06 10:39:37 readBinaries.py:l37 INFO	Getting the current pangenome's status
2021-12-06 10:39:37 readBinaries.py:l78 INFO	Extracting and writing CDS sequences from a .h5 pangenome file to a fasta file...
100%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 85940/85940 [00:01<00:00, 78048.92gene/s]
2021-12-06 10:39:38 cluster.py:l187 INFO	Clustering all of the genes sequences...
2021-12-06 10:39:38 cluster.py:l45 INFO	Creating sequence database...
2021-12-06 10:39:39 cluster.py:l54 INFO	Clustering sequences...
2021-12-06 10:39:41 cluster.py:l56 INFO	Extracting cluster representatives...
2021-12-06 10:39:41 cluster.py:l68 INFO	Writing gene to family informations
2021-12-06 10:39:41 cluster.py:l195 INFO	Associating fragments to their original gene family...
2021-12-06 10:39:41 cluster.py:l30 INFO	Aligning cluster representatives...
2021-12-06 10:39:45 cluster.py:l35 INFO	Extracting alignments...
2021-12-06 10:39:45 cluster.py:l97 INFO	Starting with 6184 families
2021-12-06 10:39:45 cluster.py:l126 INFO	Ending with 4945 gene families
2021-12-06 10:39:45 cluster.py:l148 INFO	Adding protein sequences to the gene families
2021-12-06 10:39:45 cluster.py:l130 INFO	Adding 85940 genes to the gene families
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 85940/85940 [00:00<00:00, 429762.91gene/s]
2021-12-06 10:39:45 cluster.py:l286 INFO	Done with the clustering
2021-12-06 10:39:45 writeBinaries.py:l499 INFO	Writing gene families and gene associations...
100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 4945/4945 [00:00<00:00, 207957.34gene family/s]
2021-12-06 10:39:45 writeBinaries.py:l501 INFO	Writing gene families information...
100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 4945/4945 [00:00<00:00, 263469.34gene family/s]
2021-12-06 10:39:45 writeBinaries.py:l421 INFO	Updating annotations with fragment information
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 86820/86820 [00:00<00:00, 333993.84gene/s]
2021-12-06 10:39:45 writeBinaries.py:l530 INFO	Done writing the pangenome. It is in file : pangenome.h5

The results are saved in the pangenome.h5 file given as input

Step 5
===============================================
**Build the pangenome graph**

ppanggolin graph --pangenome pangenome.h5 --cpu 8 

2021-12-06 10:46:22 main.py:l180 INFO	Command: /home/betterlab/.conda/envs/Pangenomics/bin/ppanggolin graph --pangenome pangenome.h5 --cpu 8
2021-12-06 10:46:22 main.py:l181 INFO	PPanGGOLiN version: 1.1.136
2021-12-06 10:46:22 readBinaries.py:l37 INFO	Getting the current pangenome's status
2021-12-06 10:46:22 readBinaries.py:l294 INFO	Reading pangenome annotations...
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 86820/86820 [00:00<00:00, 382162.71gene/s]
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 20/20 [00:00<00:00, 21.38organism/s]
2021-12-06 10:46:23 readBinaries.py:l307 INFO	Reading pangenome gene families...
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 85940/85940 [00:00<00:00, 254976.13gene/s]
100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 4945/4945 [00:00<00:00, 176288.18gene family/s]
2021-12-06 10:46:23 makeGraph.py:l56 INFO	Computing the neighbors graph...
Processing Mtb_N3913_L7: 100%|██████████████████████████████████████████████████████████████████████████████████████████████████████████| 20/20 [00:00<00:00, 70.11organism/s]
2021-12-06 10:46:24 makeGraph.py:l74 INFO	Done making the neighbors graph.
2021-12-06 10:46:24 writeBinaries.py:l508 INFO	Writing the edges...
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 5686/5686 [00:00<00:00, 130941.02edge/s]
2021-12-06 10:46:24 writeBinaries.py:l530 INFO	Done writing the pangenome. It is in file : pangenome.h5

The results are saved in the pangenome.h5 file given as input

Step 6
===============================================
**Pangenome partition**

This is the step that will assign gene families to the 'persistent', 'shell', or 'cloud' partitions.
The one parameter that might be of importance is the '-K', or '--nb_of_partitions' parameter. This will define the number of classes used to partition the pangenome. This may be of use if you expect to have well-defined subpopulations in your pangenome and you know exactly how many. If not, that number if detected automatically through an ICL criterion. The idea is that the most present partition will be 'persistent', the least present will be 'cloud', and all the others will be 'shell'. The number of partitions corresponding to the shell will be the number of expected subpopulations in your pangenome. (So if you expect 5 subpopulations, you could use -K 7).

In most cases, you should let the statistical criterion used by PPanGGOLiN find the optimal number of partitions for you.

ppanggolin partition --pangenome pangenome.h5 --cpu 8

2021-12-06 11:09:56 main.py:l180 INFO	Command: /home/betterlab/.conda/envs/Pangenomics/bin/ppanggolin partition --pangenome pangenome.h5 --cpu 8 --force
2021-12-06 11:09:56 main.py:l181 INFO	PPanGGOLiN version: 1.1.136
2021-12-06 11:09:56 readBinaries.py:l37 INFO	Getting the current pangenome's status
2021-12-06 11:09:56 readBinaries.py:l294 INFO	Reading pangenome annotations...
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 86820/86820 [00:00<00:00, 387372.86gene/s]
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 20/20 [00:00<00:00, 21.51organism/s]
2021-12-06 11:09:57 readBinaries.py:l307 INFO	Reading pangenome gene families...
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 85940/85940 [00:00<00:00, 256223.45gene/s]
100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 4945/4945 [00:00<00:00, 177297.84gene family/s]
2021-12-06 11:09:57 readBinaries.py:l314 INFO	Reading the neighbors graph edges...
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 83033/83033 [00:00<00:00, 192353.08contig adjacency/s]
2021-12-06 11:09:58 partition.py:l356 INFO	Estimating the optimal number of partitions...
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████| 19/19 [00:01<00:00, 18.50Number of number of partitions/s]
2021-12-06 11:09:59 partition.py:l358 INFO	The number of partitions has been evaluated at 3
2021-12-06 11:09:59 partition.py:l376 INFO	Partitioning...
2021-12-06 11:09:59 partition.py:l436 INFO	Partitionned 20 genomes in 0.35 seconds.
2021-12-06 11:09:59 writeBinaries.py:l408 INFO	Updating gene families with partition information
100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 4945/4945 [00:00<00:00, 162664.27gene family/s]
2021-12-06 11:09:59 writeBinaries.py:l530 INFO	Done writing the pangenome. It is in file : pangenome.h5

All the results will be added to the given 'pangenome.h5' input file.

Step 7
===============================================
**Predict the regions of genome plasticity with RGP module**

ppanggolin rgp --pangenome pangenome.h5 --cpu 8

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


You also can obtain a list of the plastic regions (RGPs) for each genome by using the module write 

ppanggolin write -p pangenome.h5 --regions --output rpg

2021-12-06 11:30:05 main.py:l180 INFO	Command: /home/betterlab/.conda/envs/Pangenomics/bin/ppanggolin write -p pangenome.h5 --regions --output rpg
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

Explore the rgp results

cd rgp/
ls

plastic_regions.tsv

head plastic_regions.tsv

region	                  organism	    contig	            start	 stop	    genes	 contigBorder	 wholeContig
MtbN0004_contig_15_RGP_0	Mtb_N0004_L3	MtbN0004_contig_15	83651	  93061	  17	   True	         False
MtbN0004_contig_16_RGP_0	Mtb_N0004_L3	MtbN0004_contig_16	33871	  50925	  23	   False	       False
MtbN0004_contig_3_RGP_0	  Mtb_N0004_L3	MtbN0004_contig_3	  178698	182189	5	     True	         False
MtbN0004_contig_58_RGP_0	Mtb_N0004_L3	MtbN0004_contig_58	8243	  18251	  16	   False	       False
MtbN0031_contig_10_RGP_0	Mtb_N0031_L2	MtbN0031_contig_10	79495 	97528	  24	   False	       False
MtbN0031_contig_14_RGP_0	Mtb_N0031_L2	MtbN0031_contig_14	81684	  86828	  6	     False	       False
MtbN0031_contig_29_RGP_0	Mtb_N0031_L2	MtbN0031_contig_29	23141	  28818	  8 	   False	       False
MtbN0031_contig_55_RGP_0	Mtb_N0031_L2	MtbN0031_contig_55	1457	  11465	  16	   False	       False
MtbN0052_contig_17_RGP_0	Mtb_N0052_L2	MtbN0052_contig_17	43797	  61393	  23	   False	       False

Return to the working directory
cd ..

Step 8
===============================================
**Compute the spots of insertion**

ppanggolin spot --pangenome pangenome.h5 --cpu 8

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

**THIS VERSION DO NOT ALLOW 'MODULE' NOR 'CONTEXT' ANALYSIS**

Step 9
===============================================
**Compute the pangenome results**

PPanGGOLiN provides multiple outputs to describe a pangenome. In this section the different outputs will be described

**I. DRAW**
    a) *U


{: .source}
{: .output}
{% include links.md %}
