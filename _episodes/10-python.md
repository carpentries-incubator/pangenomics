---
title: "Introduction to Python"
teaching: 30
exercises: 15
questions:
- ""
objectives:
- ""
keypoints:
- ""
---

## Libraries and builtin functions

~~~
import pandas as pd
from matplotlib import cm
import gudhi
import os  
~~~
{: .language-python}

~~~
len(genes)
~~~
{: .language-python}
~~~
43
~~~
{: .output}

## Assigning variables. Getting attributes. Joining strings.
~~~

url = "https://raw.githubusercontent.com/paumayell/pangenomics/gh-pages/files/mini-genomes.blast"
blastE = pd.read_csv(url, sep='\t', names=['qseqid', 'sseqid', 'evalue'])
 
~~~
{: .language-python}

~~~
qseqid_unique=pd.unique(blastE['qseqid'])
sseqid_unique=pd.unique(blastE['sseqid'])
genes = pd.unique(np.append(qseqid_unique, sseqid_unique))
~~~
{: .language-python}

## Vectors, dataframes, lists, dictionaries, arrays
~~~
df_genes=pd.DataFrame(genes, columns=['Genes'])
df_genome_genes=df_genes["Genes"].str.split("|", n = 1, expand = True)
df_genome_genes.columns= ["Genome", "Gen"]
genomes=pd.unique(df_genome_genes['Genome'])
genomes=list(genomes)
genomes
~~~
{: .language-python}

~~~
name_columns = matrixE2.columns
name_columns
~~~
{: .language-python}

 `numpy` array.
~~~
DistanceMatrix = matrixE2.to_numpy()
DistanceMatrix
~~~
{: .language-python}

~~~
array([[1.24e-174, 5.00e+000, 5.00e+000, ..., 5.00e+000, 5.00e+000,
        5.00e+000],
       [5.00e+000, 9.58e-100, 5.00e+000, ..., 5.00e+000, 5.00e+000,
        5.00e+000],
       [5.00e+000, 5.00e+000, 0.00e+000, ..., 5.00e+000, 5.00e+000,
        5.00e+000]
~~~
{: .output}

## Functions

~~~
def dimension(list):
    return (len(list[0])-1, list[1])
~~~
{: .language-python}

## For loops, if statements

~~~
d_simplex_time = dict()
d_simplex_const = dict()
names = []
for tuple_simple in all_simplex_sorted_dim_1:
    list_aux = []
    if len(tuple_simple[0])-1 == simplexTree.dimension(): 
        t_birth = tuple_simple[1]
        t_death = max_edge_length
        d_simplex_time[tuple(tuple_simple[0])] = (t_birth,t_death)
        list_aux = tuple([name_columns[tuple_simple[0][i]] for i in range(len(tuple_simple[0]))])
        d_simplex_const[list_aux] = (t_birth,t_death)
    else:
        t_birth = tuple_simple[1] 
        t_death = max_edge_length
        for simplex in d_simplex_time.keys():
            if set(tuple_simple[0]).issubset(set(simplex)):
                t_death = d_simplex_time[simplex][0] 
        d_simplex_time[tuple(tuple_simple[0])] = (t_birth,t_death)
        list_aux = tuple([name_columns[tuple_simple[0][i]] for i in range(len(tuple_simple[0]))])
        d_simplex_const[list_aux] = (t_birth,t_death) 
~~~
{: .language-python}

> ## Exercise 1: 
>
>  
> > ## Solution
> >
> > 
> {: .solution}
{: .challenge}

