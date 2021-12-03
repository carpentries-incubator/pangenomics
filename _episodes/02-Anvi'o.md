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



## The basic process to construct a pangenome starting with genbank files


Connect to the working server using the *ssh* command and enter the password
~~~
ssh betterlab@132.248.196.38
~~~
{: .source}

~~~
betterlab@132.248.196.38's password:
~~~
{: .output}

Type the password. **Note.** When typing the password, it will not show in the terminal for security reasons. Don't panic, is normal!
Once conexion has been established correctly, the header of your terminal will content the server information instead your personal computer information. 
~~~
(base) betterlab@betterlabub:~$
~~~

{: .output}

Create a new directory into the Pangenomics directory
~~~
cd Pangenomics/
~~~

~~~
mkdir GBK-MTBC
~~~


Move into the new directory GBK-MTBC and to obtain the complete path to this directory

~~~
cd GBK-MTBC/
pwd
~~~


Disconect to server
~~~
exit
~~~

In your computer, move to the direcotry which contains the genomes of interest and verify they are correct
~~~
cd /home/arya/Documents/Paulina/2Pangenome/Annotation_GBK/
ls
~~~

Copy your genomes to the new directory you created into the server 
~~~
scp *.gbk betterlab@132.248.196.38:/home/betterlab/Pangenomics/GBK-MTBC/.
~~~

Conect to server again
~~~
ssh betterlab@132.248.196.38
~~~

Verify you uploads
~~~
cd Pangenomics/GBK-MTBC/
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
{: .output}
{% include links.md %}
