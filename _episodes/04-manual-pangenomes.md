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

## Make blastp

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

Finally, we need to run the blastp.
~~~
blastp -query ~/pan_workshop/results/blast/mini/mini-genomes.faa -db ~/pan_workshop/results/blast/mini/database/mini-genomes -outfmt "6" > ~/pan_workshop/results/blast/mini/output-blast/mini-genomes.blast
~~~
{: .language-bash}

## Explore blast matrix

For this section, we will use Python. First we need to read the blast file that we produced.  The libraries that we will use are the following.

~~~
import os 
import pandas as pd
from matplotlib import cm
import numpy as np
~~~
{: .language-python}

The columns of the blast file are 'qseqid','sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'.

~~~
os.getcwd()
blast0 = pd.read_csv( '~/pan_workshop/results/blast/mini/output_blast/mini-genomes.blast', sep = '\t',names = ['qseqid','sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])  
~~~
{: .language-python}

We want to work with the e-values, therefore we select that column.

~~~
blastE = pd.DataFrame(blast0,columns=['qseqid','sseqid','evalue'])
blastE
~~~
{: .language-python}


|	 |qseqid	|sseqid	|evalue|
|-----|------|------|------|
|0	| `2603V|GBPINHCM_01420`	| `NEM316|AOGPFIKH_01528`	| 4.110000e-67 |
|1	| `2603V|GBPINHCM_01420`	| `A909|MGIDGNCP_01408`	| 4.110000e-67 |
|2 |	`2603V|GBPINHCM_01420` |	`515|LHMFJANI_01310` |	4.110000e-67 |
|3 |	`2603V|GBPINHCM_01420` |	`2603V|GBPINHCM_01420` |	4.110000e-67 |
|4 |	`2603V|GBPINHCM_01420` |	`A909|MGIDGNCP_01082`	| 1.600000e+00 |
|...|	...	|...	|...|
|389	| `NEM316|AOGPFIKH_00403`	| `A909|MGIDGNCP_01343`	 | 6.700000e+00 |


We will modify this data frame to obtain two new columns, one for the genomes of the `qseqid` gene and one for the `sseqid` gene. First, we obtain the genome of each gene in the `qseqid`.

~~~
qseqid = pd.DataFrame(blast0,columns=['qseqid'])
sseqid = pd.DataFrame(blast0,columns=['sseqid'])

newqseqid = qseqid["qseqid"].str.split("|", n = 1, expand = True)
newqseqid.columns= ["Genome1", "Gen"]
newqseqid["qseqid"]= qseqid
dfqseqid =newqseqid[['Genome1','qseqid']]

dfqseqid
~~~
{: .language-python}

The data frame `dfqseqid` have all the genomes associated to the genes in the column `qseqid`. Now we repeat the same for the `sseqid` column.

~~~
newsseqid = sseqid["sseqid"].str.split("|", n = 1, expand = True)
newsseqid.columns= ["Genome2", "Gen"]
newsseqid["sseqid"]= sseqid 
dfsseqid = newsseqid[['Genome2','sseqid']]

dfsseqid
~~~
{: .language-python}

We combine these two data frames and the `evalue` of the `blast0` data frame.

~~~
evalue = pd.DataFrame(blast0, columns=['evalue'])
df = dfqseqid
df['Genome2']=dfsseqid['Genome2']
df['sseqid']=sseqid
df['evalue']=evalue
df
~~~
{: .language-python}

| 	| Genome1	| qseqid	| Genome2	| sseqid	| evalue |
|-----|------|------|------|-------|-------|
|0	|2603V	| `2603V|GBPINHCM_01420`	| NEM316	| `NEM316|AOGPFIKH_01528`	| 4.110000e-67 |
|1	|2603V	| `2603V|GBPINHCM_01420`	| A909	| `A909|MGIDGNCP_01408`	| 4.110000e-67 |
|2	|2603V	| `2603V|GBPINHCM_01420`	| 515	| `515|LHMFJANI_01310`	| 4.110000e-67 |
|3	|2603V	| `2603V|GBPINHCM_01420`	| 2603V	| `2603V|GBPINHCM_01420`	| 4.110000e-67 |
|4	|2603V	| `2603V|GBPINHCM_01420`	| A909	| `A909|MGIDGNCP_01082`	| 1.600000e+00 |
|...	|...	| ...	|...	|...	|...|
|389	|NEM316	| `NEM316|AOGPFIKH_00403`	| A909	| `A909|MGIDGNCP_01343`	| 6.700000e+00 |

We want a list of the unique genes in our dataset.
~~~
qseqid_unique=pd.unique(df['qseqid'])
sseqid_unique=pd.unique(df['sseqid'])
genes = pd.unique(np.append(qseqid_unique, sseqid_unique))
~~~
{: .language-python}

We can check that we only have 43 genes with `len(genes)`.

We will fix the evalue that we will use to form the families.
~~~
evalue= 1e-5
~~~
{: .language-python}
 
Now, we want to know what is the biggest genome to make the comparisions. In this case we wil chose the one with more genes, that happen to be the A909.  
First, we compute the unique genomes.

~~~
genomes=pd.unique(df['Genome1'])
genomes=list(genomes)
genomes
~~~
{: .language-python}

Now, we will create a dictionary with the genes associated to each genome.

~~~
dic_gen_genomes={}
for a in genomes:
    temp=[]
    for i in range(len(genes)):
        if a in genes[i]:
            gen=genes[i]
            temp.append(gen)
    dic_gen_genomes[a]=temp
~~~
{: .language-python}

To compute the genome size, we do the following.
 
~~~
genome_temp=[]
size_genome=[]
for i in dic_gen_genomes.keys():
    size=len(dic_gen_genomes[i])
    genome_temp.append(i)
    size_genome.append(size)

genomes_sizes = pd.DataFrame(genome_temp, columns=['Genome'])
genomes_sizes['Size']=size_genome

genome_sizes_df = genomes_sizes.sort_values('Size', ascending=False)
genomes_sizes_df
~~~
{: .language-python}

| 	| Genome	| Size |
|------|--------|--------|
| 2	| A909 |	12 |
| 0	| 2603V |	11 |
| 1	| 515 |	10 |
| 3	| NEM316 |	10 |

So, the biggest genome es `A909` and we will start our algorithm with this genome.

The algorith that we will use is called bidirectional best-hit (BDBH).

<a href="{{ page.root }}/fig/bdbh.png">
   <img src="{{ page.root }}/fig/bdbh.png" alt=" Bidirectional best-hit algorithm" />
  </a>


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
cd ~/pan_workshop/results/annotated/mini/
ls *.faa | while read line ; do name=$(echo $line | cut -d'_' -f3); count=$(grep -c $name ~/pan_workshop/results/blast/mini/mini-genomes.faa); echo $count $name; done |sort -nr
~~~
{: .language-bash}

~~~
12 A909
11 2603V
10 NEM316
10 515
~~~
{: .output}


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
