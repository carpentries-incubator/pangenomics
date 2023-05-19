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

# Make blastp

We are going to practice with four fasta files with reduced genomes from A909, 2603V, NEM316 and 515. 

In the folder `anottated/subset` you have the 4 reduced genomes. Fist we need to put a label on each gene that tells us which genome it is from, this will be important later.  For that you need to run the following.

~~~
cd ~/pan_workshop/results/annotated/subset/
ls *.faa | while read line ; do name=$(echo $line | cut -d'_' -f3); sed -i "s/\s*>/>${name}|/" $line; done
~~~
{: .language-bash}

Now, we need to create one data set with all files `.faa`.

~~~
cat ~/pan_workshop/results/annotated/*.faa > all-genomes.faa
~~~
{: .language-bash}

We will create the folders to make the blast dataset and to run blastp. Also we move the file `all-genomes.faa` to this new directory.

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

Run blastp.
~~~
blastp -query ~/pan_workshop/results/blast/all-genomes.faa -db ~/pan_workshop/results/blast/database/all-genomes -outfmt "6" > ~/pan_workshop/results/blast/output-blast/all-genomes.blast
~~~
{: .language-bash}

# Explore blast matrix

Now we will explore the blast matrix in Python. First we need to reed the database.

~~~
import os 
os.getcwd()
blast0 = pd.read_csv( '~/pan_workshop/results/blast/output_blast/all-genomes.blast', sep = '\t',names = ['qseqid','sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])  
~~~
{: .language-bash}

We want to work with the e-values, therefore we select that column.

~~~
blastE = pd.DataFrame(blast0,columns=['qseqid','sseqid','evalue'])
blastE
~~~
{: .language-python}


We want a list of the unique genes in our dataset.

~~~
qseqid_unique=pd.unique(blastE['qseqid'])
sseqid_unique=pd.unique(blastE['sseqid'])
genes = pd.unique(np.append(qseqid_unique, sseqid_unique))
~~~
{: .language-python}

We will fix the evalue that we will use to form the families.

~~~
evalue= 1e-5
~~~
{: .language-bash}
 

Now, we want to know what is the biggest genome to make the comparisions. In this case we wil chose the one with more genes, that happen to be the A909.  

~~~
ls *.faa | while read line ; do name=$(echo $line | cut -d'_' -f3); count=$(grep -c $name ~/pan_workshop/results/blast/all-genomes.faa); echo $count $name; done |sort -nr
~~~
{: .language-bash}

 

# Count genes in genomes.

Lets explore the small genomes content. 

Here we have the functional families provided by prokka for the A909 genome
~~~~
>A909|MGIDGNCP_01408 30S ribosomal protein S16
>A909|MGIDGNCP_00096 50S ribosomal protein L16
>A909|MGIDGNCP_01343 Replication protein RepB
>A909|MGIDGNCP_01221 Glycine betaine transporter OpuD
>A909|MGIDGNCP_01268 hypothetical protein
>A909|MGIDGNCP_00580 UDPNacetylglucosamineNacetylmuramyl(pentapeptide) pyrophosphorylundecaprenol Nacetylglucosamine transferase
>A909|MGIDGNCP_00352 Glutamate 5kinase 1
>A909|MGIDGNCP_00064 Putative Nacetylmannosamine6phosphate 2epimerase
>A909|MGIDGNCP_00627 hypothetical protein
>A909|MGIDGNCP_01082 Periplasmic murein peptidebinding protein
>A909|MGIDGNCP_00877 hypothetical protein
>A909|MGIDGNCP_00405 Ribosome hibernation promotion factor
~~~~
{: .output}

First we select the biggest genome, and in this case we wil chose the one with more genes, that happen to be the A909.  
~~~
ls *.faa | while read line ; do name=$(echo $line | cut -d'_' -f3); count=$(grep -c $name ~/pan_workshop/results/subset/blast2/all-genomes.faa); echo $count $name; done |sort -nr
~~~
{: .language-bash}


### Core  
Ribosomal proteins 30S (01408) and 50S (00096) stays in its own cluster (1 gene per genome)
01268 hypothetical protein in an extended family (2 copies per genome) 

### Shell
01343 Replication protein RepB is only shared by another genome
01221 Glycine betaine transporter OpuD is only shared by another genome

### cloud
00580 UDPNacetylglucosamineNacetylmuramyl(pentapeptide) pyrophosphorylundecaprenol Nacetylglucosamine transferase
00352 Glutamate 5kinase 1
00064 Putative Nacetylmannosamine6phosphate 2epimerase
00627 hypothetical protein
01082 Periplasmic murein peptidebinding protein
00877 hypothetical protein
00405 Ribosome hibernation promotion factor
