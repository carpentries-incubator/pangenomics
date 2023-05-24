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

In the folder `anottated/mini` you have the 4 reduced genomes. Fist we need to put a label on each gene that tells us which genome it is from, this will be important later.  For that you need to run the following.

~~~
cd ~/pan_workshop/results/annotated/mini/
ls *.faa | while read line ; do name=$(echo $line | cut -d'_' -f3); sed -i "s/\s*>/>${name}|/" $line; done
~~~
{: .language-bash}

Now, we need to create one data set with all files `.faa`.

~~~
cat ~/pan_workshop/results/annotated/mini/*.faa > mini-genomes.faa
~~~
{: .language-bash}

We will create the folders to make the blast dataset and to run blastp. Also we move the file `all-genomes.faa` to this new directory.

~~~
mkdir -p ~/pan_workshop/results/blast/mini/output-blast/
mkdir -p ~/pan_workshop/results/blast/mini/database/
mv ~/pan_workshop/results/annotated/mini/mini-genomes.faa ~/pan_workshop/results/blast/mini/.
cd ~/pan_workshop/results/blast/mini/
~~~
{: .language-bash}

Make the BLAST database.
~~~
makeblastdb -in ~/pan_workshop/results/blast/mini/mini-genomes.faa -dbtype prot -out ~/pan_workshop/results/blast/mini/database/mini-genomes 
~~~
{: .language-bash}

~~~
Building a new DB, current time: 05/23/2023 21:26:31
New DB name:   /home/haydee/pan_workshop/results/blast/mini/database/mini-genomes
New DB title:  /home/haydee/pan_workshop/results/blast/mini/mini-genomes.faa
Sequence type: Protein
Keep MBits: T
Maximum file size: 1000000000B
Adding sequences from FASTA; added 43 sequences in 0.00112104 seconds.
~~~
{: .output}

Run blastp.
~~~
blastp -query ~/pan_workshop/results/blast/mini/mini-genomes.faa -db ~/pan_workshop/results/blast/mini/database/mini-genomes -outfmt "6" > ~/pan_workshop/results/blast/mini/output-blast/mini-genomes.blast
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
>MGIDGNCP_01408 30S ribosomal protein S16
>MGIDGNCP_00096 50S ribosomal protein L16
>MGIDGNCP_01343 Replication protein RepB
>MGIDGNCP_01221 Glycine betaine transporter OpuD
>MGIDGNCP_01268 glycosyltransferase
>MGIDGNCP_00580 UDP-N-acetylglucosamine--N-acetylmuramyl-(pentapeptide) pyrophosphoryl-undecaprenol N-acetylglucosamine transferase
>MGIDGNCP_00352 Glutamate 5-kinase 1
>MGIDGNCP_00064 Putative N-acetylmannosamine-6-phosphate 2-epimerase
>MGIDGNCP_00627 bifunctional DNA primase/polymerase
>MGIDGNCP_01082 Periplasmic murein peptide-binding protein
>MGIDGNCP_00877 peptidase U32 family protein
>MGIDGNCP_00405 Ribosome hibernation promotion factor
~~~~
{: .output}

First we select the biggest genome, and in this case we wil chose the one with more genes, that happen to be the A909.  
~~~
ls *.faa | while read line ; do name=$(echo $line | cut -d'_' -f3); count=$(grep -c $name ~/pan_workshop/results/subset/blast2/all-genomes.faa); echo $count $name; done |sort -nr
~~~
{: .language-bash}


### Core  
Ribosomal proteins 30S (01408) and 50S (00096) stays in its own cluster (1 gene per genome)
01268 glycosyltransferase protein in an extended family (2 copies per genome) 

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
