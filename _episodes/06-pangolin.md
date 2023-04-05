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

PPanGGOLiN has a criterion to make the partitions that is more complex than the usual partitions. FIXME

## Building a Pangenome

PPanGGOLiN analysis can start from genomic DNA sequences ([.fasta](https://raw.githubusercontent.com/paumayell/pangenomics/gh-pages/files/ExampleFASTA.fasta)) or annotated genomes ([.gbk](https://raw.githubusercontent.com/paumayell/pangenomics/gh-pages/files/ExampleGBK.gbk)) of whole genomes, Metagenomic Assembled Genomes (MAG), and Single-cell Amplified Genomes (SAG), useful for large-scale environmental studies, including the non-cultivable species pangenome.  It is designed to scale up to tens of thousands of genomes.

In addition, PPanGGOLiN includes the panRGP method (Bazin et al. 2020) that predicts Regions of Genomic Plasticity (RGP) for each genome. RGPs are groups of genes made of shell and cloud genomes in the pangenome chart, most of which arise from horizontal gene transfer and correspond to genomic islands. RGPs from different genomes are then grouped into insertion sites based on their conserved persistent flanking genes.


Before starting using PPanGGOLiN, activate the Pangenomics environment.

~~~
$ conda activate Pangenomics_Global
~~~
{: .language-bash}

~~~
(Pangenomics_Global) ~$
~~~
{: .output}


### Genome annotation


~~~
$ mkdir -p  ~/pan_workshop/results/pangenome/ppanggolin
~~~
{: .language-bash}

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

Create a symbolic link with all the `.gbk` files in the directory of the PPanGGolin analysis to have easier access to them.
~~~
$ cd ~/pan_workshop/results/pangenome/ppanggolin
$ ln -s ~/pan_workshop/results/annotated/*.gbk .
~~~
{: .language-bash}

Each line of this file represents one organism, the first column contains a unique organism name and the second column contains the path to the associate `.gbk` file.

~~~
$ ls Streptococcus_agalactiae* | cut -d'.' -f1|while read line; do echo $line$'\t'$line.gbk >> organisms.gbk.list; done
~~~
{: .language-bash}

~~~
$ ls
$ cat organisms.gbk.list
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

Using the organisms list, the annotation of genomes is made with the `annotate` module of PPanGGOLiN.

~~~
$ ppanggolin annotate --anno organisms.gbk.list --output pangenome
~~~
{: .language-bash}

Now a new directory named `pangenome/`  was created, let's move into it and explore it. PPanGGolin created a special file `pangenome.h5` that will be 
used as input and output in all of the steps. Since it will be getting enriched we can monitor it's increase in size.
~~~
$ cd pangenome/
$ ls -lh
~~~
{: .language-bash}

~~~
total 8.9M
-rw-r--r-- 1 user user 8.9M mar 31 09:41 pangenome.h5
~~~
{: .output}

### Gene clustering

PPanGGolin uses by default MMseqs2  algoritm to run the clustering in all the proteins. You can provide your own gene families or cluster 
adding the flag `--clusters` and providing a `tsv` file with your families but previously in step 4 you need to provide the annotations. 

~~~
$ ppanggolin cluster --pangenome pangenome.h5 --cpu 8
~~~
{: .language-bash}

You can change the cluster parameters adding the flags `--coverage`(default 0.8) `--identity`(default 0.8) to the previous command. You can also provide your 

We can now notice that the size of out file has increased.

~~~
$ ls -lh
~~~
{: .language-bash}

~~~
total 9.6M
-rw-r--r-- 1 user user 9.6M mar 31 09:47 pangenome.h5
~~~
{: .output}

### Build a pangenome graph and partitions

In order to obtain the partitios of the pangenome, you need to construct the pangenome graph. You can specify if you want to eliminate genes families that are too duplicated with the option `remove-high-copy-number`.
~~~
$ ppanggolin graph --pangenome pangenome.h5 --cpu 8
~~~
{: .language-bash}

This is the step that will assign gene families to the 'persistent', 'shell', or 'cloud' partitions.
The one parameter that might be of importance is the `-K`, or `--nb_of_partitions` parameter. This will define the number of classes used to partition the pangenome. This may be of use if you expect to have well-defined subpopulations in your pangenome and you know exactly how many. If not, that number is detected automatically through an Integrated Completed Likelihood (ICL) criterion. The idea is that the most present partition will be 'persistent', the least present will be 'cloud', and all the others will be 'shell'. The number of partitions corresponding to the shell will be the number of expected subpopulations in your pangenome. (So if you expect 5 subpopulations, you could use `-K 7`).

In most cases, you should let the statistical criterion used by PPanGGOLiN find the optimal number of partitions for you.

~~~
$ ppanggolin partition --pangenome pangenome.h5 --cpu 8
~~~
{: .language-bash}

## Extracting results

### Print presence/absence files

You can also export some tables with the summary of the genes and genes families that are in each partition of the pangenome. 

With `ppanggolin write -p pangenome.h5 --stats` you can obtain all the statistics of the organisms. With the command `ppanggolin write -p pangenome.h5 --Rtab` you can obtain a presence/absence of the genes in each partition, if there is a 1 then the gene family is present in the genome and 0 otherwise. The command `ppanggolin write -p pangenome.h5 --csv` produce a `csv` file with the matrix associated to the presence/absence genes, it follows the same format that the roary `gene_presence_absence.csv` and it also works with scoary.
To make a file with the summary statistics: `ppanggolin info -p pangenome.h5 --content > stats/summary_statistics.txt`

### Draw plots and interactive graph

PPanGGOLiN provides multiple outputs to describe a pangenome. Let's create all of the plots and then we will move them to our 
local machines to view them.

* **U-shaped plot**

The U-shaped plot represents the number of gene families (y axis) per number of organisms (x axis). It is an `.html` file that can be opened 
with any browser and with which you can interact, zoom, move around, mouseover to see numbers in more detail, and you can save what you are 
seeing as a `.png` image file.

~~~
$ ppanggolin draw --pangenome pangenome.h5 --ucurve --output draw_ucurve
~~~
{: .language-bash}

* **Tile plot**

A tile plot is a heatmap representing the gene families (y axis) in the organisms (x axis) making up your pangenome. The tiles on the graph 
will be colored if the gene family is present in an organism and uncolored if absent. The gene families are ordered by partition, and the 
genomes are ordered by a hierarchical clustering based on their shared gene families (basically two genomes that are close together in 
terms of gene family composition will be close together in the figure).

This plot is quite helpful to observe potential structures in your pangenome, and can also help you to identify eventual outliers. You 
can interact with it, and mousing over a tile in the plot indicates which is the gene identifier(s), the gene family and the organism 
that corresponds to the tile.

If you build your pangenome using the 'workflow' subcommand and you have more than 500 organisms, only 
the 'shell' and the 'persistent' partitions will be drawn, leaving out the 'cloud' as the figure tends to be too heavy for a browser to open it otherwise.

~~~
$ ppanggolin draw --pangenome pangenome.h5 --tile_plot --output draw_tile
~~~
{: .language-bash}

If you do not want the 'cloud' gene families to be displayed, as it is a lot of data and can be hard to open with a browser, you can use the following option:

~~~
$ ppanggolin draw --pangenome pangenome.h5 --tile_plot --nocloud --output draw_tile_nocloud
~~~
{: .language-bash}

~~~
...
2023-03-27 20:06:23 tile_plot.py:l57 INFO       done with making the dendrogram to order the organisms on the plot
2023-03-27 20:06:23 tile_plot.py:l92 INFO       Getting the gene name(s) and the number for each tile of the plot ...
2023-03-27 20:06:23 tile_plot.py:l101 INFO      Done extracting names and numbers. Making the heatmap ...
2023-03-27 20:06:23 tile_plot.py:l157 INFO      Drawing the figure itself...
2023-03-27 20:06:24 tile_plot.py:l159 INFO      Done with the tile plot : 'draw_tile_nocloud/tile_plot.html'
~~~
{: .output}

Draw the interactive graph

~~~
$ ppanggolin write -p pangenome.h5 --gexf --output gexf
~~~
{: .language-bash}


### Visualize the results

Let's see what new files were created and move them to out local machine to view them.

~~~
$ tree
~~~
{: .language-bash}

~~~
.
├── draw_tile
│   └── tile_plot.html
├── draw_tile_nocloud
│   └── tile_plot.html
├── draw_ucurve
│   └── Ushaped_plot.html
├── gexf
│   └── pangenomeGraph.gexf
├── pangenome.h5
├── rgp
│   └── plastic_regions.tsv
└── spots
    ├── spots.tsv
    └── summarize_spots.tsv
~~~
{: .output}

Open a new terminal locally. Then move to the desired directories where the images and graph file will be downloaded.
~~~
$ cd /Desktop/Workshop/
~~~
{: .language-bash}

Copy the directories with the plots and graph with `scp`.
~~~
$ scp -r user@server-address:~/pan_workshop/results/pangenome/ppanggolin/pangenome/draw* .
$ scp -r user@server-address:~/pan_workshop/results/pangenome/ppanggolin/pangenome/gexf .

~~~
{: .language-bash}


To view the plots you can open the `html` files locally in the browser of your choice.

<a href="../fig/01-06-02.png">
  <img src="../fig/01-06-02.png" width="960" height="438" alt="Bar graph depicting the gene family frequency distribution, represented by a U-shaped plot.
                                                           	The number of organisms is plotted in the x axis and the number of gene families in the y axis." />
</a>

> ## Discussion: Partitions
>   Why do are two partitions in the same bar? FIXME
> > ## Solution
> > Because in PPanGGOLiN partitions not only depend on the number of genomes a gene ir present in, but also the conservation of the neighborhood of the gene.
> {: .solution}
{: .challenge}


<a href="../fig/01-06-03.png">
  <img src="../fig/01-06-03.png" width="956.5" height="453.5" alt="Tile plot displaying the gene families present within six strains of Streptococcus agalactiae, including the cloud gene families" />
</a>

<a href="../fig/01-06-04.png">
  <img src="../fig/01-06-04.png" width="956.5" height="434.5" alt="Tile plot displaying the gene families present within six strains of Streptococcus agalactiae, in the absence of cloud gene families" />
</a>


To view the interactive graph we will use the software **gephi**.

> ## Gephi setup
> Install gephi from it [web page](https://gephi.org/).   
> Open gephi:   
> > ## Linux
> > Go to the directory where you installed the program and type:
> > ~~~
> > ./gephi-0.10.1/bin/gephi
> > ~~~
> > {: .language-bash}  
> > If you do not see the graph properly you may have problems with the video driver. Open it this way instead:  
> > ~~~
> >  LIBGL_ALWAYS_SOFTWARE=1 ./gephi-0.10.1/bin/gephi
> >  ~~~
> >  {: .language-bash}  
> {: .solution}
> 
> > ## Windows
> > Windos way of oppening gephi
> {: .solution}
> >   
> Change the language to english to make it easier to find the options to choose.  
{: .prereq}
  
Open the file `pangenomeGraph.gexf`.
Go to the layout section and in the selection bar choose the ForceAtlas2.
In Tunning section mark the stronger gravity box and set the scale in 4000.
<a href="../fig/01-06-05.png">
  <img src="../fig/01-06-05.png" width="512" height="512" alt="Gephi visualization" />
</a>

> ## Exercise 3: Exploring the pangenome graph.
> Explore the options of visualization for the pangenome graph, try yo color the nodes according to:
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

> ## Exercise 1: PPanGGolin pipeline.
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


## Plasticity regions and spots of insertion

Region of Genome Plasticity (RGP) correspond to genomic islands, plasmid and regions that have been lost in multiples strains. You can do this analysis directly from your fasta files using the command `ppanggolin panrgp`. To predict the RGP after we perfom the partition we use the following command.

~~~
$ ppanggolin rgp --pangenome pangenome.h5 --cpu 8
~~~
{: .language-bash}

To obtain a file with the list of the plastic regions (RGPs) for each genome you can use the write module.

~~~
$ ppanggolin write -p pangenome.h5 --regions --output rgp
~~~
{: .language-bash}

Use the tree command to see everything that was created in our directory.
~~~
$ tree
~~~
{: .language-bash}

~~~
.
├── pangenome.h5
└── rgp
    └── plastic_regions.tsv
~~~
{: .output}
We now have a directory named `rgp` and the file inside, let's view its contents
~~~
$ head rgp/plastic_regions.tsv
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


The RGPs that was found in the same area in different genomes can be gather into spots of insertions. Those spots are groups of RGPs that 
have similar bordering persistent genes but not necessarialy the same gene content. This analysis allows us to study the dynamic of gene 
turnover of large regions in bacterial genomes. 

~~~
$ ppanggolin spot --pangenome pangenome.h5 --cpu 8
~~~
{: .language-bash}

You can obtain a file with the list of the spots for each genome by using the module write.

~~~
$ ppanggolin write -p pangenome.h5 --spots --output spots
~~~
{: .language-bash}

Explore the spots results.
~~~
$ tree
~~~
{: .language-bash}

~~~
.
├── pangenome.h5
├── rgp
│   └── plastic_regions.tsv
└── spots
    ├── spots.tsv
    └── summarize_spots.tsv
~~~
{: .output}

Let's explore the two files that were created in the `spots/` directory.
~~~
$ head spots/spots.tsv
~~~
{: .language-bash}

~~~
spot_id	rgp_id
spot_26	NZ_HG939456.1_RGP_2
spot_23	NZ_HG939456.1_RGP_4
spot_15	NC_007432.1_RGP_1
spot_31	NC_004368.1_RGP_1
spot_34	NC_004368.1_RGP_12
spot_25	NZ_HG939456.1_RGP_5
spot_16	NC_007432.1_RGP_5
spot_16	NZ_AAJQ01000021.1_RGP_0
spot_33	NC_004368.1_RGP_4
~~~
{: .output}

~~~
$ head spots/summarize_spots.tsv
~~~
{: .language-bash}

~~~
spot	nb_rgp	nb_families	nb_unique_family_sets	mean_nb_genes	stdev_nb_genes	max_nb_genes	min_nb_genes
spot_8	7	6	1	6	0.0	6	6
spot_2	6	15	4	4.5	0.837	6	4
spot_5	6	40	5	13.833	4.262	20	8
spot_3	5	93	5	36.6	17.315	54	16
spot_4	4	47	3	14.75	8.617	26	8
spot_19	3	40	3	22.667	2.082	25	21
spot_16	2	19	2	18	0.0	18	18
spot_0	2	107	2	66	1.414	67	65
spot_18	2	7	1	7.5	0.707	8	7
~~~
{: .output}

> ## Discussion
> What is the difference between RGP regions and spots of insertion?
>
> How can you use this information?
> > ## Solution
> > The RGPs are genomic islands, plasmid and regions that have been lost in multiple strains and the spots of insertions are groups of RGPs. 
> > Those analysis are useful to study the dynamic of gene turnover of large regions in bacterial genomes. Then, spots of the same pangenome can be compared and if we compare the different metrics together we can stablish the dynamic.
> {: .solution}
{: .discussion}

> ## References:
> For more details you can check this article:
>
> Gautreau G et al. (2020) PPanGGOLiN: Depicting microbial diversity via a partitioned pangenome graph. PLOS Computational Biology 16(3): e1007732. [https://doi.org/10.1371/journal.pcbi.1007732](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007732).
{: .callout}

{% include links.md %}



