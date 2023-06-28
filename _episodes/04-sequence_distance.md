---
title: "Distance between sequences"
teaching: 30
exercises: 15
questions:
- "How con we measure differences in gene sequences?"
objectives:
- "Calculate a score between two sequences using BLAST."
keypoints:
- "To build a pangenome you need to compare the genes and build gene families."
- "BLAST gives a score of similarity between two sequences."
---

## Finding gene families

In the [previous episode](https://paumayell.github.io/pangenomics/03-annotation-with-Prokka/index.html), we annotated all of our genomes, so now we know the genes that each individual genome has (and their protein sequences). 
To build a pangenome we 
need to figure out which genes to compare between genomes. For this, we need to build **gene families**, which are groups of homologous genes (i.e. 
genes with a common ancestor). Homology between genes is found through sequence similarity, and sequence similarity is measured by aligning the
sequences and measuring the percentage of identity. The process of building gene families is called clustering.

> ## Pizza pangenomics
> Do Roma Tomatoes and Cherry Tomatoes go in the same family? 
>  
> If two genes of different species come from a gene in an ancestral species, they are **orthologs**. And if a gene
> duplicates within a species, the two resulting genes are **paralogs**. Depending on your research questions you may want to have the
> paralogs separated into different families or in the same family with duplications.
> Paralogs tend to have a higher percentage of identity than orthologs, so if you want to separate the paralogs you can use an algorithm
> that makes the families using an identity threshold and set a high threshold.
>
> Since you want to offer the most variety of ingredients in your pizza restaurant it may be a better idea to use a higher identity threshold
> (however you are measuring the similarity between ingredients)
> that separates the Roma Tomatoes and the
> Cherry Tomatoes into different families, instead of having only one family of just _tomatoes_. 
{: .discussion}

In this episode, we will demonstrate how we measure the similarity of genes using [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
and use an algorithm to group the genes into 
gene families. This is usually done by software that automates these steps, but we will do it step by step with a reduced version of 
four of our genomes to understand how this process works. Later on in the lesson, 
we will repeat these steps but in an automated way with pangenomics software using the complete genomes.


## Aligning the protein sequences to each other with BLASTp

To do our small pangenome "by hand" we will use only some of the protein sequences for these four genomes A909, 2603V, NEM316, and 515. 
In the folder `data/annotated_mini` you have the 4 reduced genomes in amino acid FASTA format. 
~~~
$ cd ~/pan_workshop/data/annotated_mini/
$ ls
~~~
{: .language-bash}
~~~
Streptococcus_agalactiae_2603V_mini.faa  Streptococcus_agalactiae_515_mini.faa  Streptococcus_agalactiae_A909_mini.faa  Streptococcus_agalactiae_NEM316_mini.faa
~~~
{: .output}

First we need to put a label on each protein to know to which genome it belongs to, this will be important later. If we explore our annotated genomes, we have amino acid sequences with a header that has the sequence ID and the functional annotation.

~~~
$ head -n1 Streptococcus_agalactiae_A909_mini.faa
~~~
{: .langauge-bash}

~~~
>MGIDGNCP_01408 30S ribosomal protein S16
~~~
{: .output}

Let's run the following to put the name of the genome in the header of each sequence.

~~~
$ ls *.faa | while read line
do 
name=$(echo $line | cut -d'_' -f3) # Take the name of the genome from the file name and remove the file extension.
sed -i "s/\s*>/>${name}|/" $line # Substitute the symbol > for > followed by the name of the genome and a | symbol.
done

$ head -n1 Streptococcus_agalactiae_A909_mini.faa
~~~
{: .language-bash}

~~~
>A909|MGIDGNCP_01408 30S ribosomal protein S16
~~~
{: .output}

Now, we need to create one dataset with the sequences from all of our genomes. We will use it to generate a database, which is a set of files that have the information of our FASTA file but in a format that BLAST can use to align the query sequences to sequences in the database.

~~~
$ cat *.faa > mini-genomes.faa
~~~
{: .language-bash}

Now let's create the folders for the BLAST database and for the `blastp` run, and move the new file `mini-genomes.faa` to this new directory.

~~~
$ mkdir -p ~/pan_workshop/results/blast/mini/output_blast/
$ mkdir -p ~/pan_workshop/results/blast/mini/database/
$ mv mini-genomes.faa ~/pan_workshop/results/blast/mini/.
$ cd ~/pan_workshop/results/blast/mini/
~~~
{: .language-bash}

Now we will make the protein database from our FASTA.
~~~
$ makeblastdb -in mini-genomes.faa -dbtype prot -out database/mini-genomes 
~~~
{: .language-bash}

~~~
Building a new DB, current time: 05/23/2023 21:26:31
New DB name:   /home/dcuser/pan_workshop/results/blast/mini/database/mini-genomes
New DB title:  /home/dcuser/pan_workshop/results/blast/mini/mini-genomes.faa
Sequence type: Protein
Keep MBits: T
Maximum file size: 1000000000B
Adding sequences from FASTA; added 43 sequences in 0.00112104 seconds.
~~~
{: .output}

Now that we have all the sequences of all of our genomes in a BLAST database we can align each of the sequences (queries) to all of the other ones  (subjects) using `blastp`.

We will ask `blastp` to align the queries to the database and give the result in the format "6", which is a tab-separated file, with the fields Query 
Sequence-ID, Subject Sequence-ID, and E-value. 
BLAST aligns the query sequence to all of the sequences in the database, it measures the percentage of identity, the percentage of the query sequence that is covered by the subject sequence, and uses these measures to give a score of how good the match is between your query and each 
subject sequence. The [E-value](https://blast.ncbi.nlm.nih.gov/doc/blast-help/FAQ.html) represents the possibility of finding a match with a similar 
score in a database of a certain size by chance. So the lower the E-value, the more significant the match between our query and the subject sequences 
is.

~~~
$ blastp -query mini-genomes.faa -db database/mini-genomes -outfmt "6 qseqid sseqid evalue" > output_blast/mini-genomes.blast
~~~
{: .language-bash}

~~~
$ head -n4 output_blast/mini-genomes.blast
~~~
 {: .language.bash}
 
~~~
2603V|GBPINHCM_01420	NEM316|AOGPFIKH_01528	4.11e-67
2603V|GBPINHCM_01420	A909|MGIDGNCP_01408	4.11e-67
2603V|GBPINHCM_01420	515|LHMFJANI_01310	4.11e-67
2603V|GBPINHCM_01420	2603V|GBPINHCM_01420	4.11e-67
~~~
{: .output}

