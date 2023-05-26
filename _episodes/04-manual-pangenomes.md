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

In the folder `anottated/mini` you have the 4 reduced genomes. Fist we need to put a label on each gene that tells us which genome it is from, this will be important later.  For that you need to run the following. If we explore our genomes, each gene do not say from which genome it belongs.

~~~
$ cd ~/pan_workshop/data/annotated_mini/
$ head -n1 Streptococcus_agalactiae_A909_mini.faa
~~~
{: .langauge-bash}

~~~
>MGIDGNCP_01408 30S ribosomal protein S16
~~~
{: .output}

If we run the following, we will put a label in each gene that says from which genome it belongs.

~~~
$ ls *.faa | while read line
do 
name=$(echo $line | cut -d'_' -f3) 
sed -i "s/\s*>/>${name}|/" $line 
done

$ head -n1 Streptococcus_agalactiae_A909_mini.faa
~~~
{: .language-bash}

~~~
>A909|MGIDGNCP_01408 30S ribosomal protein S16
~~~
{: .output}

Now, we need to create one data set with all files `.faa`.

~~~
$ cat *.faa > mini-genomes.faa
~~~
{: .language-bash}

We will create the folders to make the blast dataset and to run blastp. Also we move the file `mini-genomes.faa` to this new directory.

~~~
$ mkdir -p ~/pan_workshop/results/blast/mini/output-blast/
$ mkdir -p ~/pan_workshop/results/blast/mini/database/
$ mv mini-genomes.faa ~/pan_workshop/results/blast/mini/.
$ cd ~/pan_workshop/results/blast/mini/
~~~
{: .language-bash}

Make the BLAST database.
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

Finally, we need to run the blastp.
~~~
$ blastp -query mini-genomes.faa -db database/mini-genomes -outfmt "6" > output-blast/mini-genomes.blast
$ head -n4 output-blast/mini-genomes.blast
~~~
{: .language-bash}

~~~
2603V|GBPINHCM_01420    NEM316|AOGPFIKH_01528   100.000 90      0       0       1       90      1       90      4.11e-67        187
2603V|GBPINHCM_01420    A909|MGIDGNCP_01408     100.000 90      0       0       1       90      1       90      4.11e-67        187
2603V|GBPINHCM_01420    515|LHMFJANI_01310      100.000 90      0       0       1       90      1       90      4.11e-67        187
2603V|GBPINHCM_01420    2603V|GBPINHCM_01420    100.000 90      0       0       1       90      1       90      4.11e-67        187
~~~
{: .output}

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
genome_sizes_df
~~~
{: .language-python}

| 	| Genome	| Size |
|------|--------|--------|
| 2	| A909 |	12 |
| 0	| 2603V |	11 |
| 1	| 515 |	10 |
| 3	| NEM316 |	10 |

~~~
genomes=genome_sizes_df['Genome'].tolist()
genomes
~~~
{: .language-python}

~~~
['A909', '2603V', '515', 'NEM316']
~~~
{: .output}

So, the biggest genome es `A909` and we will start our algorithm with this genome.

The algorith that we will use is called bidirectional best-hit (BDBH).

<a href="{{ page.root }}/fig/bdbh.png">
   <img src="{{ page.root }}/fig/bdbh.png" alt=" Bidirectional best-hit algorithm" />
  </a>


This algorithm identify homologous DNA sequences between pairs of genomes. We need to find for each gene in the biggest genome the bidirectional best hit in each genome, i.e., the gene that has the smallest evalue with respect to the fixed gene in the biggest genome. 

We will use two function, one to obtain the best hit for a fixed genome (the biggest) and other to otain the bidirectional best-hits.

~~~
def besthit(gen,genome,data):
    # gen: a fixed gen in the list of unique genes
    # genome: the genome in which we will look the best hit
    # df: the data frame with the evalues
    filtro=(data['qseqid']==gen) & (data['Genome2']==genome) & (data['Genome1']!=genome)
    if (len(data[filtro]) == 0 ):
        gen_besthit = "NA"
    else:
        gen_besthit = data.loc[filtro,'sseqid'].at[data.loc[filtro,'evalue'].idxmin()]
   
    return(gen_besthit)
~~~
{: .language-python}


~~~
def besthit_bdbh(gengenome,listgenomes,genome,data):
    # gengenome: a list with all the genes of the biggest genome.
    # listgenomes: the list with all the genomes in order.
    # genome: the genome to which the genes in `gengenome` belongs.
    # data: the data frame with the evalues.
    
    dic_besthits = {}
    for a in gengenome:
        temp=[]
        for b in listgenomes:
            temp2=besthit(a,b,data)
            temp3=besthit(temp2,genome,data)
            if temp3 == a:
                temp.append(temp2)
            else:
                temp.append('NA')
        dic_besthits[a]=temp
        
    return(dic_besthits)
~~~
{: .language-python}

In one of the previuos step we create a dicctionary with all the genes associated to each genome, as we know that the biggest genome is `A909` we will obtain the genes that belongs to `A909` and we will almacenated in a list.

~~~
genome_A909 = dic_gen_genomes['A909']
~~~
{: .language-python}

Now, we will apply the function `besthit_bdbh` to the previous list, `genomes` and  the genome `A909` that is `genomes[0]`.

~~~
g_A909_bdbh=besthit_bdbh(genome_A909,genomes,genomes[0],df)
~~~
{: .language-python}

We can convert this dicctionary to a data frame.

~~~
family_A909=pd.DataFrame(g_A909_bdbh).transpose()
family_A909.columns = ['g_A909','g_2603V','g_515','g_NEM316']
family_A909.g_A909 = family_A909.index
family_A909
~~~
{: .language-python}

In this step, we have all the families that contain one gene from the biggest genome. The following step is repeat this for the second biggest genome but before we need to remove from the `genes` list the genes that appears in the families that we obtained.


~~~
lista=[]
for elemt in g_A909_bdbh.keys():
    lista.append(elemt)
    for aaa in g_A909_bdbh[elemt]:
        lista.append(aaa)
~~~
{: .language-python}

~~~
genes2=genes
genes2=genes2.tolist()
genesremove=pd.unique(lista).tolist()
genesremove.remove('NA')
for a in genesremove:
    genes2.remove(a)

genes2
~~~
{: .language-python}

~~~
['2603V|GBPINHCM_00748', '2603V|GBPINHCM_01226', '515|LHMFJANI_01625', 'NEM316|AOGPFIKH_01842']
~~~
{: .output}

For this 4 genes we will repeat the algorithm. First, we create the list with the genes that belongs to the second biggest genome `2603V`.

~~~
genome_2603V=[]
for i in range(len(genes2)):
    if "2603V" in genes2[i]:
        gen = genes2[i]
        genome_2603V.append(gen)
        
genome_2603V
~~~
{: .language-python}

We apply the function `besthit_bdbh` to this list.

~~~
g_2603V_bdbh=besthit_bdbh(genome_2603V,genomes,genomes[1],df)
g_2603V_bdbh
~~~
{: .language-python}

We create the data frame.

~~~
family_2603V=pd.DataFrame(g_2603V_bdbh).transpose()
family_2603V.columns = ['g_A909','g_2603V','g_515','g_NEM316']
family_2603V.g_2603V = family_2603V.index
family_2603V
~~~
{: .language-python}

Again, we eliminated the genes from the last list and repeat the algorithm.

~~~
for a in genome_2603V:
    genes2.remove(a)

genes2
~~~
{: .language-python}

~~~
['515|LHMFJANI_01625', 'NEM316|AOGPFIKH_01842']
~~~
{: .output}

~~~
genome_515=[]
for i in range(len(genes2)):
    if "515" in genes2[i]:
        gen = genes2[i]
        genome_515.append(gen)
        
genome_515
~~~
{: .language-python}

~~~
g_515_bdbh=besthit_bdbh(genome_515,genomes,genomes[2],df)
g_515_bdbh
~~~
{: .language-python}

~~~
family_515=pd.DataFrame(g_515_bdbh).transpose()
family_515.columns = ['g_A909','g_2603V','g_515','g_NEM316']
family_515.g_515 = family_515.index
family_515
~~~
{: .language-python}


As in this last step we use all the genes, then we finish our algorithm. We will only create a final data frame.

~~~
families_bdbh=pd.concat([family_A909,family_2603V,family_515])
families_bdbh.to_csv('families_bdbh.csv')
families_bdbh
~~~
{: .language-python}

| Family_id | g_A909	| g_2603V |	g_515	| g_NEM316 |
|----|----|----|----|----|
| `A909|MGIDGNCP_01408` |	`A909|MGIDGNCP_01408`	| `2603V|GBPINHCM_01420` |	`515|LHMFJANI_01310`	| `NEM316|AOGPFIKH_01528` |
| `A909|MGIDGNCP_00096`	| `A909|MGIDGNCP_00096` |	`2603V|GBPINHCM_00097` | `515|LHMFJANI_00097` |	`NEM316|AOGPFIKH_00098` |
| `A909|MGIDGNCP_01343`	| `A909|MGIDGNCP_01343`	| NA |	NA | `NEM316|AOGPFIKH_01415` |
| `A909|MGIDGNCP_01221`	| `A909|MGIDGNCP_01221`	| NA	| `515|LHMFJANI_01130`	| NA |
| `A909|MGIDGNCP_01268` |	`A909|MGIDGNCP_01268` |	`2603V|GBPINHCM_01231` | `515|LHMFJANI_01178`	| `NEM316|AOGPFIKH_01341` |
| `A909|MGIDGNCP_00580` |	`A909|MGIDGNCP_00580` |	`2603V|GBPINHCM_00554`	| `515|LHMFJANI_00548` |	`NEM316|AOGPFIKH_00621` |
| `A909|MGIDGNCP_00352` |	`A909|MGIDGNCP_00352` |	`2603V|GBPINHCM_00348` |	`515|LHMFJANI_00342` |	`NEM316|AOGPFIKH_00350` |
| `A909|MGIDGNCP_00064` |	`A909|MGIDGNCP_00064` |	`2603V|GBPINHCM_00065` |	`515|LHMFJANI_00064` |	`NEM316|AOGPFIKH_00065` |
| `A909|MGIDGNCP_00627` |	`A909|MGIDGNCP_00627` |	NA	| NA	| NA |
| `A909|MGIDGNCP_01082` |	`A909|MGIDGNCP_01082` |	`2603V|GBPINHCM_01042` | NA | NA |
| `A909|MGIDGNCP_00877` |	`A909|MGIDGNCP_00877`	| `2603V|GBPINHCM_00815` | `515|LHMFJANI_00781`	| `NEM316|AOGPFIKH_00855` |
| `A909|MGIDGNCP_00405` |	`A909|MGIDGNCP_00405`	| `2603V|GBPINHCM_00401` | `515|LHMFJANI_00394` | `NEM316|AOGPFIKH_00403` |
| `2603V|GBPINHCM_00748` |	NA |	`2603V|GBPINHCM_00748` |	NA | NA |
| `2603V|GBPINHCM_01226` |	NA	| `2603V|GBPINHCM_01226` |	NA |NA |
| `515|LHMFJANI_01625` |	NA |	NA | `515|LHMFJANI_01625`	| `NEM316|AOGPFIKH_01842` |


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
