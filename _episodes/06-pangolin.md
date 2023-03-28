---
title: "Neighboring Analysis of Gene Families in a Pangenome"
teaching: 20
exercises: 30
questions:
- "How can you predict the regions of genome plasticity?"
- "How can you identify the spots of insertions?"
- "How can you visualize the relationship between gene families?"

objectives:
- "Understand the fundamentals of the PPanGGOLiN tool."
- "Identify the principal differences between PPanGGOLiN and other pangenome tools."
- "Conduct a basic workflow with PPanGGOLiN."
- "Interpret the main results of PPanGGOLiN."
keypoints:
- "PPanGGOLiN is a sotfware to create and manipulate prokaryotic pangenomes."
- "PPanGGOLiN integrates protein-coding genes and their genomic neighborhood to build a graph."
- "panRGP method predicts Regions of Genomic Plasticity which are grouped into insertion sites based on their conserved persistent flanking genes."
- "PPanGGOLiN is designed to scale up to tens of thousands of genomes, including whole genomes, metagenomes and single-cell annotated genomes."
---

> ## Requirements:
> Install [gephi](https://gephi.org/) in your local machine to visualize the graphs.
{: .prereq}

<a href="../fig/01-04-01.png">
  <img src="../fig/01-04-01.png" alt="PPanGGOLiN's logo." />
</a>

# PPanGGOLiN

**Partitioned PanGenome Graph Of Linked Neighbors**

PPanGGOLiN is a software to create and manipulate prokaryotic pangenomes. It partitions a pangenome into persistent-, shell- and cloud-gene families through a graphical model and a statistical approach rather than using fixed thresholds. Unlike other methods, PPanGGOLiN integrates information about protein-coding genes and their genomic neighborhood to build a graph of gene families. Each node in the graph is a gene family and the edges represent a relation of genetic contiguity. Therefore, two gene families that are consistent neighbors in the graph are more likely to belong to the same partition, yielding a partitioned pangenome graph (PPG) made up of persistent, shell, and cloud nodes. The resulting plot looks like a subway map, where the rails represent the genomes. The following table shows how the classes are defined.

|    	Classes   		 |                           	Definition                         		 |
|:---------------------:    |:---------------------------------------------------------------------:    |
| **Persistent genome**     |      	For gene families present in almost all genomes.    		 |
|	**Shell genome**  	 | For gene families present at intermediate frequencies in the genomes.     |
|	**Cloud genome**  	 |   	For gene familes present at low frequency in the species.  		 |

## Input files

PPanGGOLiN analysis can start from genomic DNA sequences ([.fasta](https://raw.githubusercontent.com/paumayell/pangenomics/gh-pages/files/ExampleFASTA.fasta)) or annotated genomes ([.gbk](https://raw.githubusercontent.com/paumayell/pangenomics/gh-pages/files/ExampleGBK.gbk)) of whole genomes, Metagenomic Assembled Genomes (MAG), and Single-cell Amplified Genomes (SAG), useful for large-scale environmental studies, including the non-cultivable species pangenome.  It is designed to scale up to tens of thousands of genomes.

In addition, PPanGGOLiN includes the panRGP method (Bazin et al. 2020) that predicts Regions of Genomic Plasticity (RGP) for each genome. RGPs are groups of genes made of shell and cloud genomes in the pangenome chart, most of which arise from horizontal gene transfer and correspond to genomic islands. RGPs from different genomes are then grouped into insertion sites based on their conserved persistent flanking genes.

## Outputs

PPanGGOLiN provides multiple outputs to describe a pangenome. In most cases it will provide with a HDF-5 file named `pangenome.h5`. This file stores all the information about your pangenome and the analyses that were run. You can extract information from this file to get a graphical representation of your data.


> ## Exercise 1: Partitions.
>  As we see before, usually the pangenome is divided into core, dispensable and accesory genome. What is the term that PPanGGOLiN authors used for core? Which are the other pangenome partitions made by PPanGGOLiN?
>   
> a) Persistent, shell and cloud-gene.
>
> b) Softcore, shell and cloud-gene.
>
> c) Extended core, soft core and shell.
>
> d) Hard core, extended core and shell.
> > ## Solution
> >a. As it was said before, PPanGGOLiN partitions a pangenome into persistent-, shell- and cloud-gene families.
> {: .solution}
{: .challenge}

## Step by step pangenome analysis with PPanGGOLiN

Before starting using PPanGGOLiN, activate the Pangenomics environment.

~~~
$ conda activate Pangenomics_Global
~~~
{: .language-bash}

~~~
(Pangenomics_Global) ~$
~~~
{: .output}


### Step 1: Create a working directory

~~~
$ mkdir -p  ~/pan_workshop/results/pangenome/ppanggolin
~~~
{: .language-bash}

### Step 2: Identify and explore the genome files

~~~
$ ls ~/pan_workshop/results/annotated/*.gbk
~~~
{: .language-bash}

~~~
Streptococcus_agalactiae_18RS21_prokka.gbk  Streptococcus_agalactiae_CJB111_prokka.gbk
Streptococcus_agalactiae_2603V_prokka.gbk   Streptococcus_agalactiae_COH1_prokka.gbk
Streptococcus_agalactiae_515_prokka.gbk     Streptococcus_agalactiae_H36B_prokka.gbk
Streptococcus_agalactiae_A909_prokka.gbk    Streptococcus_agalactiae_NEM316_prokka.gbk
~~~
{: .output}

Create a Symbolic link with the file *.gbk*  (remember the previous episode)
~~~
$ cd ~/pan_workshop/results/pangenome/ppangolin
$ find ~/pan_workshop/results/annotated/. -name "*.gbk*" -exec ln -s {} . ';'
~~~
{: .language-bash}

### Step 3: Obtain the genome list

Each line of this file represents one organism, the first column contains a unique organism name and the second column contains the path to the associate `.gbk` file.

~~~
$ ls Streptococcus_agalactiae* | cut -d'.' -f1|while read line; do echo $line$'\t'$line.gbk >> organisms.gbk.list; done
~~~
{: .language-bash}

Move to the working directory.
~~~
$ cd
$ cd ~/pan_workshop/results/pangenome/ppangolin
$ ls
$ head organisms.gbk.list
~~~
{: .language-bash}

~~~
Streptococcus_agalactiae_18RS21_prokka	Streptococcus_agalactiae_18RS21_prokka.gbk
Streptococcus_agalactiae_2603V_prokka	Streptococcus_agalactiae_2603V_prokka.gbk
Streptococcus_agalactiae_515_prokka	Streptococcus_agalactiae_515_prokka.gbk
Streptococcus_agalactiae_A909_prokka	Streptococcus_agalactiae_A909_prokka.gbk
Streptococcus_agalactiae_CJB111_prokka	Streptococcus_agalactiae_CJB111_prokka.gbk
Streptococcus_agalactiae_COH1_prokka	Streptococcus_agalactiae_COH1_prokka.gbk
Streptococcus_agalactiae_H36B_prokka	Streptococcus_agalactiae_H36B_prokka.gbk
Streptococcus_agalactiae_NEM316_prokka	Streptococcus_agalactiae_NEM316_prokka.gbk
~~~
{: .output}

### Step 4: Genome annotation

Using the organisms list, the annotation of genomes is made with the `annotate` module of PPanGGOLiN.

~~~
$ ppanggolin annotate --anno organisms.gbk.list --output pangenome
~~~
{: .language-bash}

~~~
2023-03-17 13:37:02 main.py:l181 INFO   PPanGGOLiN version: 1.1.136
2023-03-17 13:37:02 annotate.py:l338 INFO       Reading organisms.gbk.list the list of organism files ...
100%|███████████████████████████████████████████| 8/8 [00:01<00:00,  7.31file/s]
2023-03-17 13:37:03 writeBinaries.py:l481 INFO  Writing genome annotations...
100%|████████████████████████████████████████| 8/8 [00:00<00:00, 279.31genome/s]
2023-03-17 13:37:03 writeBinaries.py:l494 INFO  writing the protein coding gene dna sequences
100%|███████████████████████████████| 16439/16439 [00:00<00:00, 142550.61gene/s]
2023-03-17 13:37:03 writeBinaries.py:l530 INFO  Done writing the pangenome. It is in file : pangenome/pangenome.h5

~~~
{: .output}

Now a new directory was created.
~~~
$ ls
~~~
{: .language-bash}

~~~
organisms.gbk.list  pangenome
~~~
{: .output}

Move into the `pangenome/` directory and explore it.
~~~
$ cd pangenome/
$ ls -lah pangenome.h5
~~~
{: .language-bash}

~~~
-rw-r--r-- 1 user rstudio-users 8.9M pangenome.h5
~~~
{: .output}

The `pangenome.h5` file will be used as input and output for all subsequent analysis.

### Step 5: Gene clustering

Ppanggolin use by default MMseqs2 to run clustering algoritm in all proteins. You can provide your own gene families or cluster adding the flag `--clusters` and providing the tsv file with your families but previously in step 4 you need to provide the annotations. 

~~~
$ ppanggolin cluster --pangenome pangenome.h5 --cpu 8
~~~
{: .language-bash}

~~~
2023-03-17 13:42:39 cluster.py:l130 INFO        Adding 16439 genes to the gene families
100%|███████████████████████████████████████████████████████████████████████████████| 16439/16439 [00:00<00:00, 619888.19gene/s]
2023-03-17 13:42:39 cluster.py:l286 INFO        Done with the clustering
2023-03-17 13:42:39 writeBinaries.py:l499 INFO  Writing gene families and gene associations...
100%|██████████████████████████████████████████████████████████████████████████| 2867/2867 [00:00<00:00, 645364.12gene family/s]
2023-03-17 13:42:39 writeBinaries.py:l501 INFO  Writing gene families information...
100%|██████████████████████████████████████████████████████████████████████████| 2867/2867 [00:00<00:00, 388519.58gene family/s]
2023-03-17 13:42:39 writeBinaries.py:l421 INFO  Updating annotations with fragment information
100%|███████████████████████████████████████████████████████████████████████████████| 17542/17542 [00:00<00:00, 473462.08gene/s]
2023-03-17 13:42:39 writeBinaries.py:l530 INFO  Done writing the pangenome. It is in file : pangenome.h5
~~~
{: .output}

You can change the cluster parameters adding the flags `--coverage`(default 0.8) `--identity`(default 0.8) to the previous command. You can also provide your 

The results are saved in the `pangenome.h5` file given as input. We can notice that the size of the file has increased.

~~~
$ ls -lah pangenome.h5
~~~
{: .language-bash}

~~~
-rw-r--r-- 1 user rstudio-users 9.6M pangenome.h5
~~~
{: .output}

### Step 6: Build the pangenome graph

In order to obtain the partitios of the pangenome, you need to construct the pangenome graph. You can specify if you want to eliminate genes families that are too duplicated with the option `remove-high-copy-number`.
~~~
$ ppanggolin graph --pangenome pangenome.h5 --cpu 8
~~~
{: .language-bash}

~~~
2023-03-27 12:50:01 readBinaries.py:l307 INFO   Reading pangenome gene families...
100%|████████████████████████████████████████████████████████████████████████████| 16439/16439 [00:00<00:00, 288907.82gene/s]
100%|███████████████████████████████████████████████████████████████████████| 2867/2867 [00:00<00:00, 210364.56gene family/s]
2023-03-27 12:50:01 makeGraph.py:l56 INFO       Computing the neighbors graph...
Processing Streptococcus_agalactiae_NEM316_prokka: 100%|████████████████████████████████| 8/8 [00:00<00:00, 138.42organism/s]
2023-03-27 12:50:01 makeGraph.py:l74 INFO       Done making the neighbors graph.
2023-03-27 12:50:01 writeBinaries.py:l508 INFO  Writing the edges...
100%|██████████████████████████████████████████████████████████████████████████████| 3222/3222 [00:00<00:00, 494730.10edge/s]
2023-03-27 12:50:01 writeBinaries.py:l530 INFO  Done writing the pangenome. It is in file : pangenome.h5
~~~
{: .output}

The results are saved in the `pangenome.h5` file given as input.


### Step 7: Pangenome partition

This is the step that will assign gene families to the 'persistent', 'shell', or 'cloud' partitions.
The one parameter that might be of importance is the `-K`, or `--nb_of_partitions` parameter. This will define the number of classes used to partition the pangenome. This may be of use if you expect to have well-defined subpopulations in your pangenome and you know exactly how many. If not, that number is detected automatically through an Integrated Completed Likelihood (ICL) criterion. The idea is that the most present partition will be 'persistent', the least present will be 'cloud', and all the others will be 'shell'. The number of partitions corresponding to the shell will be the number of expected subpopulations in your pangenome. (So if you expect 5 subpopulations, you could use `-K 7`).

In most cases, you should let the statistical criterion used by PPanGGOLiN find the optimal number of partitions for you.

~~~
$ ppanggolin partition --pangenome pangenome.h5 --cpu 8
~~~
{: .language-bash}

~~~
2023-03-27 12:54:01 readBinaries.py:l307 INFO   Reading pangenome gene families...
100%|████████████████████████████████████████████████████████████████████████████| 16439/16439 [00:00<00:00, 313436.90gene/s]
100%|███████████████████████████████████████████████████████████████████████| 2867/2867 [00:00<00:00, 232032.22gene family/s]
2023-03-27 12:54:01 readBinaries.py:l314 INFO   Reading the neighbors graph edges...
100%|████████████████████████████████████████████████████████████████| 15609/15609 [00:00<00:00, 179387.96contig adjacency/s]
2023-03-27 12:54:01 partition.py:l343 WARNING   The number of selected organisms is too low (8 organisms used) to robustly partition the graph
2023-03-27 12:54:01 partition.py:l356 INFO      Estimating the optimal number of partitions...
100%|████████████████████████████████████████████████████████████| 19/19 [00:00<00:00, 40.75Number of number of partitions/s]
2023-03-27 12:54:01 partition.py:l358 INFO      The number of partitions has been evaluated at 3
2023-03-27 12:54:01 partition.py:l376 INFO      Partitioning...
2023-03-27 12:54:01 partition.py:l436 INFO      Partitionned 8 genomes in 0.08 seconds.
2023-03-27 12:54:01 writeBinaries.py:l408 INFO  Updating gene families with partition information
100%|███████████████████████████████████████████████████████████████████████| 2867/2867 [00:00<00:00, 227639.75gene family/s]
2023-03-27 12:54:01 writeBinaries.py:l530 INFO  Done writing the pangenome. It is in file : pangenome.h5
~~~
{: .output}

All the results will be added to the given `pangenome.h5` input file.


### Step 8: Plasticity regions

Region of Genome Plasticity (RGP) correspond to genomic islands, plasmid and regions that have been lost in multiples strains. You can do this analysis directly from your fasta files using the commando `ppanggolin panrgp`. To predict the RGP after we perfom the partition we use the following command.

~~~
$ ppanggolin rgp --pangenome pangenome.h5 --cpu 8
~~~
{: .language-bash}

~~~
2023-03-27 12:59:42 genomicIsland.py:l197 INFO  Detecting multigenic families...
2023-03-27 12:59:42 pangenome.py:l311 INFO      45 gene families are defined as being multigenic. (duplicated in more than 0.05 of the genomes)
2023-03-27 12:59:42 genomicIsland.py:l199 INFO  Compute Regions of Genomic Plasticity ...
100%|████████████████████████████████████████████████████████████████████████████████████| 8/8 [00:00<00:00, 292.24genomes/s]
2023-03-27 12:59:42 genomicIsland.py:l204 INFO  Predicted 105 RGP
2023-03-27 12:59:42 writeBinaries.py:l517 INFO  Writing Regions of Genomic Plasticity...
100%|██████████████████████████████████████████████████████████████████████████████| 105/105 [00:00<00:00, 251514.52region/s]
2023-03-27 12:59:42 writeBinaries.py:l530 INFO  Done writing the pangenome. It is in file : pangenome.h5
~~~
{: .output}

To obtain a list of the plastic regions (RGPs) for each genome you can use the write module.

~~~
$ ppanggolin write -p pangenome.h5 --regions --output rgp
~~~
{: .language-bash}

~~~
100%|████████████████████████████████████████████████████████████████████████████| 17542/17542 [00:00<00:00, 395866.18gene/s]
100%|████████████████████████████████████████████████████████████████████████████████████| 8/8 [00:00<00:00, 49.21organism/s]
2023-03-27 13:00:52 readBinaries.py:l307 INFO   Reading pangenome gene families...
100%|████████████████████████████████████████████████████████████████████████████| 16439/16439 [00:00<00:00, 309036.96gene/s]
100%|███████████████████████████████████████████████████████████████████████| 2867/2867 [00:00<00:00, 230936.02gene family/s]
2023-03-27 13:00:52 readBinaries.py:l320 INFO   Reading the RGP...
100%|██████████████████████████████████████████████████████████████████████████████| 1848/1848 [00:00<00:00, 470532.01gene/s]
(Pangenomics_Global) alumno4@betterlabub:~
~~~
{: .output}

Explore the rgp results.

~~~
$ cd rgp/
$ ls
~~~
{: .language-bash}

~~~
plastic_regions.tsv
~~~
{: .output}

~~~
$ head plastic_regions.tsv
~~~
{: .language-bash}

~~~
region                  organism                                contig          start   stop    genes   contigBorder    wholeContig
AAJO01000011.1_RGP_0    Streptococcus_agalactiae_18RS21_prokka  AAJO01000011.1  6863    27451   20      True            False
AAJO01000013.1_RGP_0    Streptococcus_agalactiae_18RS21_prokka  AAJO01000013.1  564     25430   36      True            True
AAJO01000034.1_RGP_0    Streptococcus_agalactiae_18RS21_prokka  AAJO01000034.1  95      5670    6       True            False
AAJO01000044.1_RGP_0    Streptococcus_agalactiae_18RS21_prokka  AAJO01000044.1  14      13435   16      True            True
AAJO01000046.1_RGP_0    Streptococcus_agalactiae_18RS21_prokka  AAJO01000046.1  156     13006   13      True            True
AAJO01000061.1_RGP_0    Streptococcus_agalactiae_18RS21_prokka  AAJO01000061.1  84      10318   9       True            True
AAJO01000073.1_RGP_0    Streptococcus_agalactiae_18RS21_prokka  AAJO01000073.1  91      5837    6       True            False
AAJO01000077.1_RGP_0    Streptococcus_agalactiae_18RS21_prokka  AAJO01000077.1  1440    7746    7       True            False
AAJO01000081.1_RGP_0    Streptococcus_agalactiae_18RS21_prokka  AAJO01000081.1  4585    8617    6       True            False
~~~
{: .output}

Return to the `pangenome/` directory.
~~~
$ cd ..
~~~
{: .language-bash}

### Step 9: Spots of insertion

The RGPs that was found in the same area in different genomes can be gather into spots of insertions. Those sports are groups of RGPs that have similar bordering persistent genes but not necessarialy the same gene content. This analysis allows us to study the dynamic of gene turnover of large regions in bacterial genomes. 

~~~
$ ppanggolin spot --pangenome pangenome.h5 --cpu 8
~~~
{: .language-bash}

~~~
2023-03-27 13:11:20 spot.py:l132 INFO   Detecting hotspots in the pangenome...
2023-03-27 13:11:20 spot.py:l82 INFO    39 RGPs were not used as they are on a contig border (or have less than 3 persistent gene families until the contig border)
2023-03-27 13:11:20 spot.py:l83 INFO    66 RGPs are being used to predict spots of insertion
2023-03-27 13:11:20 spot.py:l85 INFO    40 number of different pairs of flanking gene families
2023-03-27 13:11:20 spot.py:l140 INFO   35 spots were detected
2023-03-27 13:11:20 writeBinaries.py:l522 INFO  Writing Spots of Insertion...
100%|██████████████████████████████████████████████████████████████████████████████████| 35/35 [00:00<00:00, 719610.98spot/s]
2023-03-27 13:11:20 writeBinaries.py:l530 INFO  Done writing the pangenome. It is in file : pangenome.h5
~~~
{: .output}

You can obtain a list of the spots for each genome by using the module write.

~~~
$ ppanggolin write -p pangenome.h5 --spots --output spots
~~~
{: .language-bash}
~~~
2023-03-27 13:11:53 readBinaries.py:l320 INFO   Reading the RGP...
100%|██████████████████████████████████████████████████████████████████████████████| 1848/1848 [00:00<00:00, 496959.27gene/s]
2023-03-27 13:11:53 readBinaries.py:l326 INFO   Reading the spots...
100%|████████████████████████████████████████████████████████████████████████████████| 66/66 [00:00<00:00, 241135.94region/s]
2023-03-27 13:11:53 writeFlat.py:l504 INFO      Done writing spots in : 'spots/summarize_spots.tsv'
~~~
{: .output}

Explore the spots results.

~~~
$ cd spots/
$ ls
~~~
{: .language-bash}

~~~
spots.tsv  summarize_spots.tsv
~~~
{: .output}

~~~
$ head spots.tsv
~~~
{: .language-bash}

~~~
spot_id rgp_id
spot_3  NC_004116.1_RGP_0
spot_3  NZ_CP051004.1_RGP_0
spot_28 NC_004368.1_RGP_3
spot_6  NC_004116.1_RGP_6
spot_15 NC_007432.1_RGP_11
spot_8  NC_004116.1_RGP_7
spot_8  NC_007432.1_RGP_3
spot_8  NC_004368.1_RGP_6
spot_8  NZ_HG939456.1_RGP_7
~~~
{: .output}

~~~
$ head summarize_spots.tsv
~~~
{: .language-bash}

~~~
spot    nb_rgp  nb_families     nb_unique_family_sets   mean_nb_genes   stdev_nb_genes  max_nb_genes    min_nb_genes
spot_2  7       6       1       6       0.0     6       6
spot_1  6       40      5       13.833  4.262   20      8
spot_5  6       15      4       4.5     0.837   6       4
spot_7  5       93      5       36.6    17.315  54      16
spot_8  4       47      3       14.75   8.617   26      8
spot_22 3       40      3       22.667  2.082   25      21
spot_3  2       107     2       66      1.414   67      65
spot_20 2       19      2       18      0.0     18      18
spot_13 2       10      1       10.5    0.707   11      10
~~~
{: .output}

Return to the `pangenome/` directory.
~~~
$ cd ..
~~~
{: .language-bash}


> ## Discussion
> What is the difference between RGP regions and spots of insertion?
>
> How can you use this information?
> > ## Solution
> > The RGPs are genomic islands, plasmid and regions that have been lost in multiple strains and the spots of insertions are groups of RGPs. 
> > Those analysis are usefult to study the dinamic of gene turnover of large regions in bacterial genomes. Then, spots of the same pangenome can be compared and if we compare the different metrics together we can stablished the dynamic.
> {: .solution}
{: .discussion}

### Step 10: Draw pangenome plots

PPanGGOLiN provides multiple outputs to describe a pangenome. In this section the different outputs will be described.

### U-shaped plot

The U-shaped plot represent the number of gene families (y axis) per number of organisms (x axis). It is a `.html` file that can be opened with any browser and with which you can interact, zoom, move around, mouseover to see numbers in more detail, and you can save what you are seeing as a `.png` image file.

~~~
$ ppanggolin draw --pangenome pangenome.h5 --ucurve --output draw_ucurve
~~~
{: .language-bash}

~~~
100%|███████████████████████████████| 16439/16439 [00:00<00:00, 312666.54gene/s]
100%|██████████████████████████| 2867/2867 [00:00<00:00, 235071.25gene family/s]
2023-03-27 19:58:49 readBinaries.py:l314 INFO   Reading the neighbors graph edges...
100%|███████████████████| 15609/15609 [00:00<00:00, 182167.72contig adjacency/s]
2023-03-27 19:58:49 ucurve.py:l13 INFO  Drawing the U-shaped curve...
2023-03-27 19:58:49 ucurve.py:l60 INFO  Done drawing the U-shaped curve : 'draw_ucurve/Ushaped_plot.html'
~~~
{: .output}

~~~
$ cd draw_ucurve/
$ ls
~~~
{: .language-bash}

~~~
Ushaped_plot.html
~~~
{: .output}

Return to the working directory.
~~~
$ cd ..
~~~
{: .language-bash}


#### Visualize the result

Open a new terminal locally. Then move to the desired directory where the images will be downloaded.
~~~
$ cd .\Desktop\Workshop\
~~~
{: .language-bash}

Copy the image to your directory using `scp` and write the password of the server.
~~~
$ scp user@bioinformatica:~/pan_workshop/results/pangenome/ppangolin/pangenome/draw_ucurve/Ushaped_plot.html .
~~~
{: .language-bash}

~~~
Ushaped_plot.html                             100% 3405KB   2.6MB/s   00:01
~~~
{: .output}

You can open the html file locally.

<a href="../fig/01-06-02.png">
  <img src="../fig/01-06-02.png" width="960" height="438" alt="Bar graph depicting the gene family frequency distribution, represented by a U-shaped plot.
                                                           	The number of organisms is plotted in the x axis and the number of gene families in the y axis." />
</a>


### Tile plot

A tile plot is a heatmap representing the gene families (y axis) in the organisms (x axis) making up your pangenome. The tiles on the graph will be colored if the gene family is present in an organism and uncolored if absent. The gene families are ordered by partition, and the genomes are ordered by a hierarchical clustering based on their shared gene families (basically two genomes that are close together in terms of gene family composition will be close together in the figure).

This plot is quite helpful to observe potential structures in your pangenome, and can also help you to identify eventual outliers. You can interact with it, and mousing over a tile in the plot indicates which is the gene identifier(s), the gene family and the organism that corresponds to the tile.

If you build your pangenome using the 'workflow' subcommand and you have more than 500 organisms, only the 'shell' and the 'persistent' partitions will be drawn, leaving out the 'cloud' as the figure tends to be too heavy for a browser to open it otherwise.

~~~
$ ppanggolin draw --pangenome pangenome.h5 --tile_plot --output draw_tile
~~~
{: .language-bash}

~~~
100%|█████████████████████████████████████████| 15609/15609 [00:00<00:00, 180972.77contig adjacency/s]
2023-03-27 20:05:53 tile_plot.py:l26 INFO       Drawing the tile plot...
2023-03-27 20:05:53 tile_plot.py:l42 INFO       start with matrice
2023-03-27 20:05:54 tile_plot.py:l57 INFO       done with making the dendrogram to order the organisms on the plot
2023-03-27 20:05:54 tile_plot.py:l92 INFO       Getting the gene name(s) and the number for each tile of the plot ...
2023-03-27 20:05:54 tile_plot.py:l101 INFO      Done extracting names and numbers. Making the heatmap ...
2023-03-27 20:05:54 tile_plot.py:l157 INFO      Drawing the figure itself...
2023-03-27 20:05:54 tile_plot.py:l159 INFO      Done with the tile plot : 'draw_tile/tile_plot.html'
~~~
{: .output}

You can download the plot and explore it in you computer.

<a href="../fig/01-06-03.png">
  <img src="../fig/01-06-03.png" width="956.5" height="453.5" alt="Tile plot displaying the gene families present within six strains of Streptococcus agalactiae, including the cloud gene families" />
</a>

If you do not want the 'cloud' gene families to be displayed, as it is a lot of data and can be hard to open with a browser, you can use the following option:

~~~
$ ppanggolin draw --pangenome pangenome.h5 --tile_plot --nocloud --output draw_tile_nocloud
~~~
{: .language-bash}

~~~
2023-03-27 20:06:23 tile_plot.py:l57 INFO       done with making the dendrogram to order the organisms on the plot
2023-03-27 20:06:23 tile_plot.py:l92 INFO       Getting the gene name(s) and the number for each tile of the plot ...
2023-03-27 20:06:23 tile_plot.py:l101 INFO      Done extracting names and numbers. Making the heatmap ...
2023-03-27 20:06:23 tile_plot.py:l157 INFO      Drawing the figure itself...
2023-03-27 20:06:24 tile_plot.py:l159 INFO      Done with the tile plot : 'draw_tile_nocloud/tile_plot.html'
~~~
{: .output}

<a href="../fig/01-06-04.png">
  <img src="../fig/01-06-04.png" width="956.5" height="434.5" alt="Tile plot displaying the gene families present within six strains of Streptococcus agalactiae, in the absence of cloud gene families" />
</a>

> ## Exercise 2: Basic commands.
>   Choose the indispensable commands to create a U-shaped plot.
>
> Commands:
> 1. cluster: Cluster proteins in protein families.
> 2. partition: Partition the pangenome graph.
> 3. rgp: Predicts Regions of Genomic Plasticity in the genomes of your pangenome.
> 4. annotate: Annotate genomes.
> 5. graph: Create the pangenome graph.
> 6. spot: Predicts spots in your pangenome.
> 7. draw: Draw figures representing the pangenome through different aspects.
>
> a) 1, 2, 3, 4, 5.
>
> b) 4, 1, 5, 7, 6.
>
> c) 4, 1, 5, 2, 7.
>
> d) 4, 2, 1, 6, 3.
> > ## Solution
> >c. The first step is always to annotate the genes, then cluster the proteins within its corresponding families, after that it is necessary to create the pangenome
> >graph, partition it and finally with the draw command create the U-shaped plot.
> {: .solution}
{: .challenge}

> ## Exercise 3: Exploring the pangenome graph.
> 1. In your terminal, execute the following command:
>
> ~~~
> $ ppanggolin write -p pangenome.h5 --gexf --output gexf
> ~~~
> {: .language-bash}
>  
> 2. With `scp` copy the produced file in your local computer.
> 3. Open the file in the Gephi program.
> 4. Go to the layout section and in the selection bar choose the ForceAtlas2.
> 5. In Tunning section mark the stronger gravity box and set the scale in 4000.
> <a href="../fig/01-06-05.png">
  <img src="../fig/01-06-05.png" width="512" height="512" alt="Gephi visualization" />
</a>
> 6. Finally color the nodes according to:
>
> a) Partition.
>  
> b) Number of organisms.
>
> c) Number of genes.
>  
> d) Proteins function (product).
>  
> e) Gene neighborhood (edges).
>  
> > ## Solution
> > a) <a href="../fig/01-04-06.png">
> > <img src="../fig/01-04-06.png" alt="Gephi visualization with the nodes colored according to PPanGGOLiN partitions, displaying persistent genome in pink and cloud genome in green." />
> > </a>
> >
> > b) <a href="../fig/01-04-07.png">
> > <img src="../fig/01-04-07.png" alt="Gephi visualization with the nodes colored according to the number of organisms in the analysis." />
> > </a>
> >
> > c) <a href="../fig/01-04-08.png">
> > <img src="../fig/01-04-08.png" alt="Gephi visualization with the nodes colored according to the number of genes in the analysis." />
> > </a>
> >
> > d) <a href="../fig/01-04-09.png">
> > <img src="../fig/01-04-09.png" alt="Gephi visualization with the nodes colored according to the proteins function." />
> > </a>
> >
> > e) <a href="../fig/01-04-10.png">
> > <img src="../fig/01-04-10.png" alt="Gephi visualization with the nodes colored according to the gene neighborhood." />
> > </a>
> {: .solution}
{: .challenge}

### Presence/absence files

You can also export some tables with the summary of the genes and genes families that are in each partition of the pangenome. 

With `ppanggolin write -p pangenome.h5 --stats` you can obtain all the statistics of the organisms. With the command `ppanggolin write -p pangenome.h5 --Rtab` you can obtain a presence/absence of the genes in each partition, if there is a 1 then the gene family is present in the genome and 0 otherwise. The command `ppanggolin write -p pangenome.h5 --csv` produce a `csv` file with the matrix associated to the presence/absence genes, it follows the same format that the roary `gene_presence_absence.csv` and it also works with scoary.

> ## References:
> For more details you can check this article:
>
> Gautreau G et al. (2020) PPanGGOLiN: Depicting microbial diversity via a partitioned pangenome graph. PLOS Computational Biology 16(3): e1007732. [https://doi.org/10.1371/journal.pcbi.1007732](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007732).
{: .callout}

{% include links.md %}



