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

## Aligning the protein sequences to each other with BLASTp

In the previous episode we annotated all of our genomes and obtained the annotations in different formats, like GFF, GBK and FASTA. 
Now we want to understand how to go from annotations to a pangenome, so we will use a reduced version of 
the annotations to make a pangenome step by step.

To make a pangenome, first we need to know how similar are our sequences from each other, we will do this by aligning each of 
our sequences to all the rest with Protein [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi).
Then we need to use these results to cluster the sequences into gene families. 

Later on in the lesson 
we will repeat these steps but in an automated way with pangenomics software using the complete genomes.

To be able to do a pangenome "by hand" we will use only some of the protein sequences for these four genomes A909, 2603V, NEM316 and 515. 
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

Now that we have all the sequences of all of our genomes in a BLAST database we can align each of the sequences (queries) to all of other ones  (subjects) using `blastp`.
~~~
$ blastp -query mini-genomes.faa -db database/mini-genomes -outfmt "6 qseqid sseqid evalue" > output_blast/mini-genomes.blast
~~~
{: .language-bash}
Here we asked `blastp` to align the queries to the database and give the result in the format "6", which is a tab separated file, with the fields Query Sequence-ID, Subject Sequence-ID and E-value, which is the measure of similarity between sequences that we need. 

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

## Processing the BLAST results

For this section we will use Python. Let's open the notebook and start by importing the libraries that we will need.
~~~
import os 
import pandas as pd
from matplotlib import cm
import numpy as np
~~~
{: .language-python}
First we need to read the `mini-genomes.blast` file that we produced.
Let's import the BLAST results to Python using the column names: `qseqid`,`sseqid`, `evalue`.

~~~
os.getcwd()
blastE = pd.read_csv( '~/pan_workshop/results/blast/mini/output_blast/mini-genomes.blast', sep = '\t',names = ['qseqid','sseqid','evalue'])
blastE.head()
~~~
{: .language-python}

~~~
  qseqid	               sseqid	               evalue
0	2603V|GBPINHCM_01420	NEM316|AOGPFIKH_01528	4.110000e-67
1	2603V|GBPINHCM_01420	A909|MGIDGNCP_01408	4.110000e-67
2	2603V|GBPINHCM_01420	515|LHMFJANI_01310	4.110000e-67
3	2603V|GBPINHCM_01420	2603V|GBPINHCM_01420	4.110000e-67
4	2603V|GBPINHCM_01420	A909|MGIDGNCP_01082	1.600000e+00
~~~
{: .output}

Now we want to make two columns that have the name of the genomes of the queries, and the name of the genomes of the subjects. We will take this information from the query and subject IDs (the label that we added at the beggining of the episode).

First, let's obtain the genome of each query gene.
~~~
qseqid = pd.DataFrame(blastE,columns=['qseqid'])

newqseqid = qseqid["qseqid"].str.split("|", n = 1, expand = True)
newqseqid.columns= ["Genome1", "Gen"]
newqseqid["qseqid"]= qseqid
dfqseqid =newqseqid[['Genome1','qseqid']]

dfqseqid.head()
~~~
{: .language-python}

~~~
  Genome1	qseqid
0	2603V	2603V|GBPINHCM_01420
1	2603V	2603V|GBPINHCM_01420
2	2603V	2603V|GBPINHCM_01420
3	2603V	2603V|GBPINHCM_01420
4	2603V	2603V|GBPINHCM_01420
~~~
{: .output}

Now let's repeat the same for the `sseqid` column.
~~~
sseqid = pd.DataFrame(blastE,columns=['sseqid'])

newsseqid = sseqid["sseqid"].str.split("|", n = 1, expand = True)
newsseqid.columns= ["Genome2", "Gen"]
newsseqid["sseqid"]= sseqid 
dfsseqid = newsseqid[['Genome2','sseqid']]
~~~
{: .language-python}

Now that we have two dataframes with the new columns that we wanted, let's combine them with the `evalue` of the `blastE` dataframe into a new one called `df`.

~~~
evalue = pd.DataFrame(blastE, columns=['evalue'])
df = dfqseqid
df['Genome2']=dfsseqid['Genome2']
df['sseqid']=sseqid
df['evalue']=evalue
df.head()
~~~
{: .language-python}

~~~
  Genome1	qseqid	Genome2	sseqid	evalue
0	2603V	2603V|GBPINHCM_01420	NEM316	NEM316|AOGPFIKH_01528	4.110000e-67
1	2603V	2603V|GBPINHCM_01420	A909	A909|MGIDGNCP_01408	4.110000e-67
2	2603V	2603V|GBPINHCM_01420	515	515|LHMFJANI_01310	4.110000e-67
3	2603V	2603V|GBPINHCM_01420	2603V	2603V|GBPINHCM_01420	4.110000e-67
4	2603V	2603V|GBPINHCM_01420	A909	A909|MGIDGNCP_01082	1.600000e+00
~~~
{: .output}

Now we want a list of the unique genes in our dataset.
~~~
qseqid_unique=pd.unique(df['qseqid'])
sseqid_unique=pd.unique(df['sseqid'])
genes = pd.unique(np.append(qseqid_unique, sseqid_unique))
~~~
{: .language-python}

We can check that we have 43 genes in total with `len(genes)`.

Now, we want to know which one is the biggest genome (the one with more genes) to make the comparisions.  

First, we compute the unique genomes.

~~~
genomes=pd.unique(df['Genome1'])
genomes=list(genomes)
genomes
~~~
{: .language-python}
~~~
['2603V', '515', 'A909', 'NEM316']
~~~
{: .output}

Now, we will create a dictionary that shows which genes are in each genome.

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

We can now use this dictionary to know how many genes does each genome has, and therefore identify the biggest genome.
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

~~~
Genome	Size
2	A909	12
0	2603V	11
1	515	10
3	NEM316	10
~~~
{: .output}

Now we can sort our genomes by their size.
~~~
genomes=genome_sizes_df['Genome'].tolist()
genomes
~~~
{: .language-python}

~~~
['A909', '2603V', '515', 'NEM316']
~~~
{: .output}

So the biggest genome es `A909` and we will start our clustering algorithm with this genome.

## Finding gene families with BDBH algorithm

To make a gene family, we first need to identify the most similar genes between genomes. The 
Bidirectional Best Hit algorithm will allow us to find the pairs of genes that are the most similar 
(lowest e-value) to each other in each pair of genomes.  

<a href="{{ page.root }}/fig/bdbh.png">
   <img src="{{ page.root }}/fig/bdbh.png" alt=" Bidirectional best-hit algorithm" />
  </a>

For this we will define a function to find in each genome the gene that is most similar to each
gene in our biggest genome A909. 

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

Now we will define a second function, that uses the previous one, to obtain the bidirectional
best hits.  

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

In one of the previous steps, we created a dictionary with all the genes present in each genome.
Since we know that the biggest genome is `A909`, we will obtain the genes belonging to `A909` and 
almacenate them in a list.  

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
family_A909.head()
~~~
{: .language-python}

~~~
                      g_A909	g_2603V	g_515	g_NEM316
A909|MGIDGNCP_01408	A909|MGIDGNCP_01408	2603V|GBPINHCM_01420	515|LHMFJANI_01310	NEM316|AOGPFIKH_01528
A909|MGIDGNCP_00096	A909|MGIDGNCP_00096	2603V|GBPINHCM_00097	515|LHMFJANI_00097	NEM316|AOGPFIKH_00098
A909|MGIDGNCP_01343	A909|MGIDGNCP_01343	NA	NA	NEM316|AOGPFIKH_01415
A909|MGIDGNCP_01221	A909|MGIDGNCP_01221	NA	515|LHMFJANI_01130	NA
A909|MGIDGNCP_01268	A909|MGIDGNCP_01268	2603V|GBPINHCM_01231	515|LHMFJANI_01178	NEM316|AOGPFIKH_01341
~~~
{: .output}

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
family_2603V.head()
~~~
{: .language-python}

~~~
                      g_A909	g_2603V	g_515	g_NEM316
2603V|GBPINHCM_00748	NA	2603V|GBPINHCM_00748	NA	NA
2603V|GBPINHCM_01226	NA	2603V|GBPINHCM_01226	NA	NA
~~~
{: .output}

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
family_515.head()
~~~
{: .language-python}

~~~
                    g_A909	g_2603V	g_515	g_NEM316
515|LHMFJANI_01625	NA	NA	515|LHMFJANI_01625	NEM316|AOGPFIKH_01842
~~~
{: .output}


As in this last step we use all the genes, then we finish our algorithm. We will only create a final data frame.

~~~
families_bdbh=pd.concat([family_A909,family_2603V,family_515])
families_bdbh.to_csv('families_bdbh.csv')
families_bdbh.head()
~~~
{: .language-python}

~~~
	                  g_A909	              g_2603V	              g_515	              g_NEM316
A909|MGIDGNCP_01408	A909|MGIDGNCP_01408	2603V|GBPINHCM_01420	515|LHMFJANI_01310	NEM316|AOGPFIKH_01528
A909|MGIDGNCP_00096	A909|MGIDGNCP_00096	2603V|GBPINHCM_00097	515|LHMFJANI_00097	NEM316|AOGPFIKH_00098
A909|MGIDGNCP_01343	A909|MGIDGNCP_01343	NA	NA	NEM316|AOGPFIKH_01415
A909|MGIDGNCP_01221	A909|MGIDGNCP_01221	NA	515|LHMFJANI_01130	NA
A909|MGIDGNCP_01268	A909|MGIDGNCP_01268	2603V|GBPINHCM_01231	515|LHMFJANI_01178	NEM316|AOGPFIKH_01341
~~~
{: .output}


Finaly, we will export to a csv file.

~~~
families_bdbh.to_csv('~/pan_workshop/results/blast/mini/families_minis.csv')
~~~
{: .language-python}

## Explore functional annotation families

The unique functional families that our mini genomes have are the following.

~~~
$ cat mini-genomes.faa | grep '>' | cut -d' ' -f2,3 | sort | uniq
~~~
{: .language-bash}

~~~
30S ribosomal
50S ribosomal
bifunctional DNA
Glutamate 5-kinase
Glycine betaine
glycosyltransferase
Glycosyltransferase GlyG
peptidase U32
Periplasmic murein
PII-type proteinase
Putative N-acetylmannosamine-6-phosphate
Replication protein
Ribosome hibernation
UDP-N-acetylglucosamine--N-acetylmuramyl-(pentapeptide) pyrophosphoryl-undecaprenol
Vitamin B12
~~~
{: .output}

Now we will look for this functional families in the families that we obtain with the BDBH algorith to see if we recover this functional families.


~~~
$ cat mini-genomes.faa | grep '>' | cut -d' ' -f2 | sort | uniq | while read function
do 
grep $function mini-genomes.faa | cut -d' ' -f1 | cut -d'>' -f2 | while read line
do 
family=$(grep $line families_minis.csv| cut -d',' -f1)
echo $function$'\t'$line$'\t'$family
done
done
~~~
{: .language-bash}


~~~
30S     2603V|GBPINHCM_01420    A909|MGIDGNCP_01408
30S     515|LHMFJANI_01310      A909|MGIDGNCP_01408
30S     A909|MGIDGNCP_01408     A909|MGIDGNCP_01408
30S     NEM316|AOGPFIKH_01528   A909|MGIDGNCP_01408
...
glycosyltransferase     2603V|GBPINHCM_01231    A909|MGIDGNCP_01268
glycosyltransferase     515|LHMFJANI_01178      A909|MGIDGNCP_01268
glycosyltransferase     A909|MGIDGNCP_01268     A909|MGIDGNCP_01268
glycosyltransferase     NEM316|AOGPFIKH_01341   A909|MGIDGNCP_01268

Glycosyltransferase     2603V|GBPINHCM_01226    2603V|GBPINHCM_01226
...
Vitamin 515|LHMFJANI_01625      515|LHMFJANI_01625
Vitamin NEM316|AOGPFIKH_01842   515|LHMFJANI_01625
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
