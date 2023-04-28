---
title: "Understanding Pangenomes with blast"
teaching: 30
exercises: 15
questions:
- "What elements do I compare to calculate a Pangenome?"
- "How are genome families conformed?"
objectives:
- "Calculate an score between two sequences using blast"
- "Cluster gene families according to a treshold"
keypoints:
- "Genes are the elements compared in pangenomes"
- "BLAST gives an score between two sequences"
- "Gene families must be clustered according to a distance-treshold "
---

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

Make the blast database.

~~~
makeblastdb -in ~/pan_workshop/results/blast/all-genomes.faa -dbtype prot -out ~/pan_workshop/results/blast/database/all-genomes 
~~~
{: .language-bash}

~~~
Building a new DB, current time: 04/27/2023 18:58:55
New DB name:   ~/pan_workshop/results/blast/database/all-genomes
New DB title: ~/pan_workshop/results/blast/all-genomes.faa
Sequence type: Protein
Keep MBits: T
Maximum file size: 1000000000B
Adding sequences from FASTA; added 16439 sequences in 0.3211 seconds.
~~~
{: .output}

Run blastp. This step takes about 20 minutes.

~~~
blastp -query ~/pan_workshop/results/blast/all-genomes.faa -db ~/pan_workshop/results/blast/database/all-genomes -outfmt "6" > ~/pan_workshop/results/blast/output-blast/all-genomes.blast
~~~
{: .language-bash}



