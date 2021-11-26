---
title: "Anvi'o"
teaching: 15 min
exercises: 40 min
questions:
- "What is Anvi'o?"
- "How can I obtain a pangenome analysis in Anvi'o?"
objectives:
- "Establish a dataset of genomes to obtain their pangenome"
- "Perform a basic workflow to obtain a pangenome in Anvi'o"
- "Understand and interpret the pangenome results"
keypoints:
- "Anvi’o is an open-source, community-driven analysis and visualization platform for microbial ‘omics. "
---
## Anvi'o

Anvi’o is an open-source, community-driven analysis and visualization platform for microbial ‘omics.
It brings together many aspects of today's cutting-edge strategies including **genomics, metagenomics, metatranscriptomics, phylogenomics, microbial population genetis, pangenomics and metapangenomis** in an *integrated* and *easy-to-use* fashion thorugh extensive interactive visualization capabilities. 


![Figure 1. Anvi'o network representation](../fig/anvio-network.png)



## Installation of Anvi'o

* Install a stable release of anvi’o on a Mac, Linux, or Windows running computer (best option for end users).
* Use anvi’o from the active development branch (best option for developers).
* Run anvi’o through its Docker containers without any installation (best option for the lazy).



## Selection of genomes 

Selection of proper data to construct a pangenome analysis is a crucial task of the process. By definition, a pangenome represents the entire set of genes from all strains in a clade. Thus, let's define and understand what a clade is. 

A **clade** (*kládos*, 'branch'), also known as *monophyletic group*, is a group of organisms that share the same common ancestor. In a phylogenetic tree, clades can be distinguished by identified the common ancestor of branch subgroups. For instance, in Fig.2, the blue and red subgroups are considered clades but green subgruop is not a clade, it is a *parahyletic group*. Why is that? Becuse it excludes the blue clade which has descended from the same common ancestor. Instead, the green and blue subgroups together form a clade.



![Figure 2. Cladogram representation](../fig/Cladogram.png)



For this lesson purposes, we selected 20 representative *Mycobacterium tuberculosis* genomes. 



## Brief description of *Mycobacterium tuberculosis* dataset

Tuberculosis (TB) is the most lethal infection disease worldwide with 1.4 million causalities annually. Members of the *Mycobacterium tuberculosis* complex (MTBC) are recognized as the causative agents of TB in humans and other mammals. The MTBC is composed by nine *Mycobacterium* species affecting wild and domesticated animals, and two species that infect principally humans *Mycobacterium tuberculosis sensu stricto* and *Mycobacterium africanum*. Interestingly, all MTBC members form a clade when compared with other *Mycobacterium* spp. They share 99.99% of their nucleotide sequences and it has been decreted a maximum distance of ~2,500 SNPs between any pair of genomes belonging to the MTBC. No horizontal gene transfer or genomic recombination have been reported in these complex of bacteria, suggesting that these organisms are 'clones'. A pangenome analysis of these bacteria can help us to determine if MTBC has an open or a closed pengenome. 


**Selected genomes**


| Organism                | Host    | Lineage    | Country     | SRA accesion number |
|-------------------------|---------|------------|-------------|---------------------|
|*M. tuberculosis*  N0004 | Human   | L3         |USA          | ERR2704697          |
|*M. tuberculosis*  N0031 | Human   | L2         |USA          | ERR2704676          |
|*M. tuberculosis*  N0052 | Human   | L2         |USA          | ERR2704699          |
|*M. tuberculosis*  N0054 | Human   | L3         |USA          | ERR2704678          |
|*M. tuberculosis*  N0069 | Human   | L1         |China        | ERR2704679          |
|*M. tuberculosis*  N0072 | Human   | L1         |USA          | ERR2704680          |
|*M. tuberculosis*  N0091 | Human   | L6         |Gambia       | ERR2704681          |
|*M. tuberculosis*  N0136 | Human   | L4         |USA          | ERR2704700          |
|*M. tuberculosis*  N0145 | Human   | L2         |USA          | ERR2704702          |
|*M. tuberculosis*  N0155 | Human   | L2         |USA          | ERR2704684          |
|*M. tuberculosis*  N0157 | Human   | L1         |USA          | ERR2704704          |
|*M. tuberculosis*  N1176 | Human   | L5         |Ghana        | ERR2704686          |
|*M. tuberculosis*  N1201 | Human   | L6         |Ghana        | ERR2704687          |
|*M. tuberculosis*  N1202 | Human   | L6         |Ghana        | ERR2704688          |
|*M. tuberculosis*  N1216 | Human   | L4         |Ghana        | ERR2704689          |
|*M. tuberculosis*  N1268 | Human   | L5         |Sierra Leone | ERR2704690          |
|*M. tuberculosis*  N1272 | Human   | L5         |Germany      | ERR2704692          |
|*M. tuberculosis*  N1274 | Human   | L3         |Germany      | ERR2704693          |
|*M. tuberculosis*  N1283 | Human   | L4         |Germany      | ERR2704709          |
|*M. tuberculosis*  N3913 | Human   | L7         |Australia    | ERR2704711          |


## The basic process to construct a pangenome starting with genbank files

Connect to server through ssh command and enter the password
~~~
ssh betterlab@132.248.196.38
~~~

**Note.** When typing the password, it will not show in the terminal for security reasons. Don't panic, is normal!


Create a new directory into the Pangenomics directory
~~~
cd Pangenomics/
~~~
~~~
mkdir GBK-MTBC

Move into the new directory GBK-MTBC and to obtain the complete path to this directory

~~~
cd GBK-MTBC/
~~~
~~~
pwd
~~~


Disconect to server
~~~
exit
~~~

In your computer, move to the direcotry which contains the genomes of interest and verify they are correct
~~~
cd /home/arya/Documents/Paulina/2Pangenome/Annotation_GBK/
~~~

~~~
ls
~~~

Copy your genomes to the new directory you created into the server 
~~~
scp *.gbk betterlab@132.248.196.38:/home/betterlab/Pangenomics/GBK-MTBC/.

Conect to server again
~~~
ssh betterlab@132.248.196.38
~~~

Verify you uploads
~~~
cd Pangenomics/GBK-MTBC/
~~~
~~~
ls
~~~

To start using ANVIO, activate the conda environment which contain it

~~~
conda activate Pangenomics
~~~~

**Important note:** Avoid including "-" symbol within the genbank file names


**STEP 1.** Process the GBK files with anvi-script-process-genbank script (in batch). This script takes a GenBank file, and outputs a FASTA file, as well as two additional TAB-delimited output files for external gene calls and gene functions that can be used with the programs anvi-gen-contigs-database and anvi-import-functions.

~~~
ls *gbk | cut -d'.' -f1 | while read line; do echo anvi-script-process-genbank -i GENBANK --input-genbank $line\.gbk -O $line; anvi-script-process-genbank -i GENBANK --input-genbank $line\.gbk -O $line; done
~~~

**STEP 2.** Reformat the fasta files 
~~~
ls *fa | while read line; do anvi-script-reformat-fasta $line -o $line\.fasta; done
~~~
**STEP 3.** Create the database 
~~~
ls *fasta | while read line; do anvi-gen-contigs-database -T 4 -f $line -o $line-contigs.db; done
~~~

**STEP 4.** Create a list of the genomes 
~~~
ls *fa |cut -d'-' -f1 |while read line; do echo $line$'\t'$line-contigs.db >>external-genomes.txt; done
~~~

**STEP 5.** Modify the headers of the list external-genomes.txt
~~~
nano external-genomes.txt
~~~
~~~
name	contigs_db_path
~~~

**STEP 6.** Rename the .db files
~~~
rename s'/.fa.fasta-contigs.db/.db/' *db
~~~

**STEP 7.** Create a database of the genomes of interest
~~~
anvi-gen-genomes-storage -e external-genomes.txt -o MTBC-rep_GENOMES.db
~~~

**STEP 8.** Create a pangenome with the list of genomes created before
~~~
anvi-pan-genome -g MTBC-rep_GENOMES.db \
                --project-name "Pan_MTBC-rep" \
                --output-dir Pangenome-MTBC \
                --num-threads 6 \
                --minbit 0.5 \
                --mcl-inflation 10 \
                --use-ncbi-blast
~~~

**Other example**
~~~
anvi-pan-genome -g MTB-133_GENOMES.db \
                --project-name "Pangenome_MTB133" \
                --output-dir Pangenome-MTB133 \
                --num-threads 6 \
                --minbit 0.5 \
                --mcl-inflation 10 \
                --use-ncbi-blast
~~~

**STEP 9.** Create the imagen of the results
~~~
anvi-display-pan -g MTBC-rep_GENOMES.db \
    -p Pangenome-MTBC/Pan_MTBC-rep-PAN.db
~~~

**STEP 10.**Visualize the results in a server page


{: .source}
{% include links.md %}
