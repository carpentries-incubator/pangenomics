---
title: "Measuring Sequence Similarity"
teaching: 30
exercises: 15
questions:
- "How can we measure differences in gene sequences?"
objectives:
- "Calculate a score between two sequences using BLAST."
keypoints:
- "To build a pangenome you need to compare the genes and build gene families."
- "BLAST gives a score of similarity between two sequences."
---

## Finding gene families

In the previous episode, we annotated all of our genomes, so now we know the genes that each individual genome has (and their protein sequences). 
To build a pangenome we 
need to figure out which genes to compare between genomes. For this, we need to build **gene families**, which are groups of homologous genes (i.e. 
genes with a common ancestor). Homology between genes is found through sequence similarity, and sequence similarity is measured by aligning the
sequences and measuring the percentage of identity. The process of building gene families is called clustering.

> ## Pizza pangenomics
> Do Roma Tomatoes and Cherry Tomatoes go in the same family? 
>  
>> ## Solution
>> 
>> If two genes of different species come from a gene in an ancestral species, they are **orthologs**. And if a gene
>> duplicates within a species, the two resulting genes are **paralogs**. Depending on your research questions you may want to have the
>> paralogs separated into different families or in the same family with duplications.
>> Paralogs tend to have a higher percentage of identity than orthologs, so if you want to separate the paralogs you can use an algorithm
>> that uses an identity threshold, and set a high threshold.
>>
>> Since you want to offer the most variety of ingredients in your pizza restaurant it may be a better idea to separate the ingredients into as many families as possible.
>> To achieve that you should use a **higher identity threshold** (whatever that means when you are comparing ingredients), this way you would separate the Roma Tomatoes and the
>> Cherry Tomatoes into two families, instead of having one family of just _tomatoes_. 
>{: .solution}
{: .discussion}

In this episode, we will demonstrate how we measure the similarity of genes using [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
and in the next one, we will use an algorithm to group the genes into 
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

First, we need to put a label on each protein to know to which genome it belongs to, this will be important later. If we explore our annotated
genomes, we have amino acid sequences with a header that has the sequence ID and the functional annotation.

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

Now, we need to create one dataset with the sequences from all of our genomes. We will use it to generate a database, 
which is a set of files that have the information of our FASTA file but in a format that BLAST can use to align the 
query sequences to sequences in the database.

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
score in a database of a certain size by chance. So the lower the E-value, the more significant the match between our query and the subject sequences is.

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

> ## Exercise 1: Remote blast search
> We already know how to perform a BLAST search of one FASTA file with many sequences to a custom database of the same sequences.  
> What if we want to make a search against the available NCBI databases?  
> 1) Search on the help page of `blastp` how you can do a remote search.  
> 2) Search on the help page of `blastp` which other fields can be part of your tabular output.  
> 3) Create a small FASTA file with only one sequence of one of our mini genomes.  
> 4) Run `blastp` remotely against the `refseq_protein` database for the created FASTA file and add more fields to the output.  
> (Note that adding the `qseqid` field will not be necessary because we are searching only one protein.)  
> 
> > ## Solution
> > Use the command-line manual of `blastp`
> > ~~~
> > $ blastp -help
> > ~~~
> >  {: .language.bash}
> > ~~~
> > -remote
> >   Execute search remotely?
> >    * Incompatible with:  gilist, seqidlist, taxids, taxidlist,
> >   negative_gilist, negative_seqidlist, negative_taxids, negative_taxidlist,
> >   subject_loc, num_threads
> > ~~~
> >  {: .output}
> > ~~~
> > Options 6, 7 and 10 can be additionally configured to produce
> >    a custom format specified by space delimited format specifiers,
> >    or by a token specified by the delim keyword.
> >     E.g.: "10 delim=@ qacc sacc score".
> >    The delim keyword must appear after the numeric output format
> >    specification.
> >    The supported format specifiers are:
> >    	    qseqid means Query Seq-id
> >    	       qgi means Query GI
> >    	      qacc means Query accesion
> >    	   qaccver means Query accesion.version
> >    	      qlen means Query sequence length
> > .
> > .
> > .
> > ~~~
> >  {: .output}
> > Print the sequence to know the identifier.
> > ~~~
> > $ head -n2 ~/pan_workshop/data/annotated_mini/Streptococcus_agalactiae_A909_mini.faa 
> > ~~~
> >  {: .language.bash}
> > ~~~
> > >A909|MGIDGNCP_01408 30S ribosomal protein S16
> > MAVKIRLTRMGSKKKPFYRINVADSRAPRDGRFIETVGTYNPLVAENQVTIKEERVLEWL
> > ~~~
> >  {: .output}
> > Create the new FASTA file with the sequence and put the identifier of the sequence in the name of the file.
> > ~~~
> > $ head -n2 ~/pan_workshop/data/annotated_mini/Streptococcus_agalactiae_A909_mini.faa > Streptococcus_agalactiae_A909_MGIDGNCP_01408.faa
> > ~~~
> >  {: .language.bash}
> > Run blast using the `-remote` flag against the `refseq_protein` database and and use different fields in the `-outfmt` option.
> > ~~~
> > $ blastp -query Streptococcus_agalactiae_A909_MGIDGNCP_01408.faa -db refseq_protein -remote -outfmt "6 sseqid evalue bitscore" > output_blast/Streptococcus_agalactiae_A909_MGIDGNCP_01408.blast
> > ~~~
> > ~~~
> > $ head output_blast/Streptococcus_agalactiae_A909_MGIDGNCP_01408.blast
> > ~~~
> >  {: .language.bash}
> > ~~~
> > ref|WP_109910314.1|	2.23e-36	126
> > ref|WP_278043300.1|	2.30e-36	126
> > ref|WP_000268757.1|	2.72e-36	126
> > ref|WP_017645295.1|	3.13e-36	126
> > ref|WP_120033169.1|	3.20e-36	126
> > ref|WP_136133384.1|	4.17e-36	125
> > ref|WP_020833411.1|	4.55e-36	125
> > ref|WP_195675206.1|	6.68e-36	125
> > ref|WP_004232185.1|	6.83e-36	125
> > ref|WP_016480974.1|	7.54e-36	125
> > ~~~
> > {: .output}
> > 
> {: .solution}
{: .challenge}
