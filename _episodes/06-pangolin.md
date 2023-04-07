---
title: "Neighboring Analysis of Gene Families in a Pangenome"
teaching: 20
exercises: 30
questions:
- "How can I build a pangenome of thousands of genomes?"
- "How can I visualize the spacial relationship between gene families?"

objectives:
- "Understand the fundamentals of the PPanGGOLiN tool."
- "Identify the principal differences between PPanGGOLiN and other pangenome tools."
- "Conduct a basic workflow with PPanGGOLiN."
- "Interpret the main results of PPanGGOLiN."
keypoints:
- "PPanGGOLiN is a sotfware to create and manipulate prokaryotic pangenomes."
- "PPanGGOLiN integrates gene families and their genomic neighborhood to build a graph and define the partitions."
- "PPanGGOLiN is designed to scale up to tens of thousands of genomes."
---

<a href="../fig/01-04-01.png">
  <img src="../fig/01-04-01.png" alt="PPanGGOLiN's logo." />
</a>

## PPanGGOLiN: Partitioned PanGenome Graph Of Linked Neighbors

[PPanGGOLiN](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007732) is a software designed to create and manipulate 
prokaryotic pangenomes in its own special way. It gives a lot more information than just 
the presence and absence of the gene families and it can be used with [thousands of genomes](https://github.com/labgem/PPanGGOLiN/wiki/Basic-usage-and-practical-information). This program **represents the pangenome in a graph** where each node is a gene family, and two gene 
families are connected by an edge if they are neighbors (their sequences are next to each other) in at least one genome. The width of this edge 
represents the number of genomes in which these two families are neighbors. And the size of the nodes represents the number of genes that are part of the
gene family.

To classify a gene family into a partition, PPanGGOLiN not only takes into account the percentage of genomes the family is present in, but it also 
considers its neighbors. If two gene families are consistently linked across the genomes, they will more likely belong to the same partition. The software uses statistical approach to decide how many partitions should be made and in which partition a gene family should be placed. 


|    	Classes   		 |                           	Definition                         		 |
|:---------------------:    |:---------------------------------------------------------------------:    |
| **Persistent genome**     |      	For gene families present in almost all genomes.    		 |
|	**Shell genome**  	 | For gene families present at intermediate frequencies in the genomes. There can be multiple shells.    |
|	**Cloud genome**  	 |   	For gene familes present at low frequency in the species.  		 |

The PPanGGOLiN pipeline can be divided into **building the pangenome**, and **extracting results**. You can also perform some **special analyses** and
extract the corresponding results.

## Building a Pangenome

The required steps to build a pangenome can be achieved with a single command `ppanggolin workflow`,
or can be run individually if you want to make adjustments to each of them. We will run them separately to learn what the workflow consists of.

The pangenome will be built in a single HDF-5 file that will be the input and output of all the commands and will get enriched with each of them.

### Genome annotation

Before starting using PPanGGOLiN, we need to activate the Pangenomics environment.
~~~
$ conda activate Pangenomics_Global
~~~
{: .language-bash}
~~~
(Pangenomics_Global) ~$
~~~
{: .output}

Now we should make a special directory for this analysis.
~~~
$ mkdir -p  ~/pan_workshop/results/pangenome/ppanggolin
~~~
{: .language-bash}

PPanGGOLiN analysis can start from genomic DNA sequences in FASTA format or annotated genomes in GBK or GFF formats. The first step is getting this 
genomic information into the HDF-5 file and annotate it if it is not already.  
We already have our Prokka annotations in `~/pan_workshop/results/annotated/`, so let's use them for our pangenome.
Create symbolic links of the `.gbk` files in the directory of the PPanGGolin analysis to have easier access to them.
~~~
$ cd ~/pan_workshop/results/pangenome/ppanggolin
$ ln -s ~/pan_workshop/results/annotated/*.gbk .
~~~
{: .language-bash}

To use the GBKs in the annotation step we need to create a text file with the unique name of each organism in one column and the path to the
corresponding `.gbk` in another one.

~~~
$ ls Streptococcus_agalactiae* | cut -d'.' -f1 | while read line # Obtain the name of the file without the file extension and open a loop
> do 
> echo $line$'\t'$line.gbk >> organisms.gbk.list # Print the name, a tab separator and the name with the extension and redirect it to a file
> done
~~~
{: .language-bash}

~~~
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

PPanGGolin uses by default [MMseqs2](https://github.com/soedinglab/MMseqs2) to run the clustering in all the proteins. You can adjust 
the clustering parameters adding the flags `--mode`, `--coverage`(default 0.8) `--identity`(default 0.8) to the command.

If you provided your annotations in the preovious step you can also provide clusters made previously with a different program 
with the `--clusters` option. 

We will use the default parameters.
~~~
$ ppanggolin cluster --pangenome pangenome.h5 --cpu 8
~~~
{: .language-bash}

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

### Build a pangenome graph

Now that the clustering step identified all the gene families we can build the graph with them. 

~~~
$ ppanggolin graph --pangenome pangenome.h5 --cpu 8
~~~
{: .language-bash}

### Partition the graph

Finally we can assign the gene families to the persistent, shell, or cloud partitions. PPanGGOLiN can find the optimal number of partitions, if it is
larger than three it will make more shell partitions. You can also specify how many partitions you want with the option `-K`.

Besides this partitions, PPanGGOLiN will also calculate the exact core (families in 100% of genomes) and exact accessory 
(families in less than 100% of genomes) and the soft core (families in more than 95% of genomes) and soft accessory (families in less than 95% of genomes).

~~~
$ ppanggolin partition --pangenome pangenome.h5 --cpu 8
~~~
{: .language-bash}

## Extracting results

### Print files with information of the pangenome

For a first glimpse of our pangenome we can obtain the summary statistics with the command `ppanggolin info`.
~~~
ppanggolin info -p pangenome.h5 --content
~~~
{: .language-bash}
~~~
Genes : 16439
Organisms : 8
Families : 2867
Edges : 3221
Persistent ( min:0.62, max:1.0, sd:0.08, mean:0.96 ): 1785
Shell ( min:0.38, max:0.62, sd:0.05, mean:0.51 ): 19
Cloud ( min:0.12, max:0.62, sd:0.13, mean:0.22 ): 1063
Number of partitions : 3
~~~
{: .output}

If we want to have this in a file we can redirect this output adding `> summary_statistics.txt` to the command.

With the `ppanggolin write` command you can extract many text files and tables with a lot of information. For this you need to provide the 
`pangenome.h5` file and the name of the directory to store the files. Each of the additional flags indicates wich file or files to write. Let's 
use all of the flags that will give us basic information of our analysis. And then see what was generated.

~~~
ppanggolin write -p pangenome.h5 --output files --stats --csv --Rtab --partitions --projection --families_tsv
~~~
{: .language-bash}
~~~
tree
~~~
{: .language-bash}
~~~
.
├── files
│   ├── gene_families.tsv
│   ├── gene_presence_absence.Rtab
│   ├── matrix.csv
│   ├── mean_persistent_duplication.tsv
│   ├── organisms_statistics.tsv
│   ├── partitions
│   │   ├── cloud.txt
│   │   ├── exact_accessory.txt
│   │   ├── exact_core.txt
│   │   ├── persistent.txt
│   │   ├── shell.txt
│   │   ├── soft_accessory.txt
│   │   ├── soft_core.txt
│   │   ├── S.txt
│   │   └── undefined.txt
│   └── projection
│       ├── Streptococcus_agalactiae_18RS21_prokka.tsv
│       ├── Streptococcus_agalactiae_2603V_prokka.tsv
│       ├── Streptococcus_agalactiae_515_prokka.tsv
│       ├── Streptococcus_agalactiae_A909_prokka.tsv
│       ├── Streptococcus_agalactiae_CJB111_prokka.tsv
│       ├── Streptococcus_agalactiae_COH1_prokka.tsv
│       ├── Streptococcus_agalactiae_H36B_prokka.tsv
│       └── Streptococcus_agalactiae_NEM316_prokka.tsv
└── pangenome.h5
~~~
{: .language-bash}


### Draw interactive plots

We can also extract two interactive plots with the command `ppanggolin draw`, follwing a similar syntax than with the command `ppanggolin write`.

~~~
$ ppanggolin draw --pangenome pangenome.h5 --output plots --ucurve --tile_plot
~~~
{: .language-bash}
~~~
tree plots
~~~
{: .language-bash}
~~~
plots/
├── tile_plot.html
└── Ushaped_plot.html
~~~
{: .output}

Let's download them to our local machines to explore them. Open a new local terminal, navigate to the directory where you want the files 
and use the `scp` command.
~~~
$ scp -r user@server-address:~/pan_workshop/results/pangenome/ppanggolin/pangenome/plots .
~~~
{: .language-bash}

Open both `html` files in a browser and explore them. They are interactive so play with the plots for a while to see what information can be 
obtained from them.

* **U-shaped plot**

The U-shaped plot is a bar chart where you can see how many gene families (y-axis) are in how many genomes (x-axis). They are colored according to the 
partition they belong to.
<a href="../fig/01-06-02.png">
  <img src="../fig/01-06-02.png" width="960" height="438" alt="Bar graph depicting the gene family frequency distribution, represented by a U-shaped plot. The number of organisms is plotted in the x axis and the number of gene families in the y axis." />

* **Tile plot**

 The tile plot is a presence/absence heatmap of the gene families (y-axis) in each genome (x-axis) ordered by a hierarchical clustering and showing the 
  multycopy families.
<a href="../fig/01-06-03.png">
  <img src="../fig/01-06-03.png" width="956.5" height="453.5" alt="Tile plot displaying the gene families present within six strains of Streptococcus agalactiae, including the cloud gene families" />

> ## Difficulty showing the plots?
> If your browser is taking too long to open your tile plot repeat the `pangolin draw` command with the option `--nocloud`, this will remove the 
> families in the cloud partition from the plot to make it easier to handle.  
{: .callout}
  
> ## Discussion: Partitions
>   Why are there bars with two partitions in the U plot?
> > ## Solution
> > Because in PPanGGOLiN the partitions not only depend on the number of genomes a gene family is present in, but also on the conservation of the 
> > neighborhood of the gene. So a gene that is in only a few genomes but has the same shell neighbors in all or most of them, will more likely be placed 
> > in the shell than in the cloud. 
> {: .solution}
{: .challenge}

### Draw the pangenome graph

Now its time to see the actual partitioned graph. For that we can use the `ppangolin write` command.
~~~
$ ppanggolin write -p pangenome.h5 --gexf --output gexf
~~~
{: .language-bash}

Let's download the `gexf` file to our local machines to explore it. Open a new local terminal, navigate to the directory where you want the file 
and use the `scp` command.
~~~
$ scp -r user@server-address:~/pan_workshop/results/pangenome/ppanggolin/pangenome/gexf .
~~~
{: .language-bash}
  
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
> 
> Change the language to English to make it easier to find the options to choose.  
{: .prereq}
  
Open the file `pangenomeGraph.gexf`.
Go to the layout section and in the selection bar choose the ForceAtlas2.
In Tunning section mark the stronger gravity box and set the scale in 4000.
<a href="../fig/01-06-05.png">
  <img src="../fig/01-06-05.png" width="512" height="512" alt="Gephi visualization" />
</a>

Each node is a gene family, which are labled with the name of a representative gene of the family.

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

## Special analyses
  
### Plasticity regions and spots of insertion

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
> Go to the PPanGGOLiN [GitHub Wiki](https://github.com/labgem/PPanGGOLiN/wiki) for the complete collection of instructions and posibilities.  
> And read the original PPanGGOLiN article to understand the details.
>
> Gautreau G et al. (2020) PPanGGOLiN: Depicting microbial diversity via a partitioned pangenome graph. PLOS Computational Biology 16(3): e1007732. [https://doi.org/10.1371/journal.pcbi.1007732](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007732).
{: .callout}

{% include links.md %}



