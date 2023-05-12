---
title: "Understanding Pangenomes with BLAST"
teaching: 30
exercises: 15
questions:
- "What elements do I compare to calculate a Pangenome?"
- "How are genome families conformed?"
objectives:
- "Calculate an score between two sequences using BLAST"
- "Cluster gene families according to a treshold"
keypoints:
- "Genes are the elements compared in pangenomes"
- "BLAST gives a score between two sequences"
- "Gene families must be clustered according to a distance-treshold "
---

Curate the `.faa` files.

~~~
ls *.faa | while read line ; do name=$(echo $line | cut -d'_' -f3); sed -i "s/\s*>/>${name}|/" $line; done
~~~
{: .language-bash}

Create one data set with all files `.faa`.
~~~
cat ~/pan_workshop/results/annotated/*.faa > all-genomes.faa
~~~
{: .language-bash}

Create the folders for `blast`.

~~~
mkdir -p ~/pan_workshop/results/blast/output-blast/
mkdir -p ~/pan_workshop/results/blast/database/
mv ~/pan_workshop/results/annotated/all-genomes.faa ~/pan_workshop/results/blast/.
~~~
{: .language-bash}

Make the BLAST database.

~~~
makeblastdb -in ~/pan_workshop/results/blast/all-genomes.faa -dbtype prot -out ~/pan_workshop/results/blast/database/all-genomes 
~~~
{: .language-bash}

~~~
Building a new DB, current time: 05/05/2023 00:13:13
New DB name:   /home/haydee/pan_workshop/results/subset/blast/database/all-genomes
New DB title:  /home/haydee/pan_workshop/results/subset/blast/all-genomes.faa
Sequence type: Protein
Keep MBits: T
Maximum file size: 1000000000B
Adding sequences from FASTA; added 16 sequences in 0.000644922 seconds.
~~~
{: .output}


Run blastp. This step takes about 20 minutes.

~~~
blastp -query ~/pan_workshop/results/blast/all-genomes.faa -db ~/pan_workshop/results/blast/database/all-genomes -outfmt "6" > ~/pan_workshop/results/blast/output-blast/all-genomes.blast
~~~
{: .language-bash}

Now we will explore the blast matrix in R.

First we need to reed the database.

`data_blast <- read.delim("~/pan_workshop/results/blast/output-blast/all-genomes.blast", header=FALSE)`

The names of the columns.

`names_data=c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")`
`colnames(data_blast) <- names_data`
  
 We will fix the evalue that we will use to form the families.
 
`evalue <-1e-5`

The uniques genes are the following.

`gen_qseqid <- unique(data_blast$qseqid)`
`gen_sseqid <- unique(data_blast$sseqid)`


Count genes in genomes.

~~~
ls *.faa | while read line ; do name=$(echo $line | cut -d'_' -f3); count=$(grep -c $name ~/pan_workshop/results/subset/blast2/all-genomes.faa); echo $count $name; done |sort -nr
~~~
{: .language-bash}

