---
title: "TDA in Pangenomes"
teaching: 30
exercises: 15
questions:
- "How can I apply TDA to describe Pangenomes"
objectives:
- "Describe Pangenomes using Gudhi"
keypoints:
- "Pangenomes can be described using TDA"
---

## Persistent approach to pangenomics

Import libraries
~~~
import pandas as pd
from matplotlib import cm
import numpy as np
import gudhi
import time 
~~~
{: .language-python}

Read blastp matrices from episode 4.
~~~
import os 
os.getcwd()
blast0 = pd.read_csv( '~/pan_workshop/results/blast/mini/output-blast/mini-genomes.blast', sep = '\t',names = ['qseqid','sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])  
~~~
{: .language-python}

Obtain evalues.

~~~
blastE = pd.DataFrame(blast0,columns=['qseqid','sseqid','evalue'])
blastE.head()
~~~
{: .language-python}

Obtain a list with genomes and genes unique.
~~~
qseqid_unique=pd.unique(blastE['qseqid'])
sseqid_unique=pd.unique(blastE['sseqid'])
genes = pd.unique(np.append(qseqid_unique, sseqid_unique))
~~~
{: .language-python}


Obtain distance matrix.

~~~
distance_list = blastE[ blastE['qseqid'].isin(genes) & blastE['sseqid'].isin(genes)]
distance_list.head()
~~~
{: .language-python}

Process the matrix.

~~~
# Define maximun biological distance between genes.
MaxDistance = 5.0000000

# reshape long to wide
#matrixE = pd.pivot_table(distance_list,index = "qseqid",values = "evalue",columns = 'sseqid', #fill_value= MaxDistance )

matrixE2 = pd.pivot_table(distance_list,index = "qseqid",values = "evalue",columns = 'sseqid')
matrixE=matrixE2.fillna(5)


# convert to a data frame
DistanceMatrix = matrixE.to_numpy()
DistanceMatrix.head()
~~~
{: .language-python}


Construct the Rips complex.

~~~
max_edge_length = 2
# Rips complex with the distance matrix
start_time = time.time()
ripsComplex = gudhi.RipsComplex(
    distance_matrix = DistanceMatrix, 
    max_edge_length = max_edge_length
)
print("The Rips complex was created in %s" % (time.time() - start_time) )
~~~
{: .language-python}

Create the filtration.

~~~
start_time = time.time()
simplexTree = ripsComplex.create_simplex_tree(
    max_dimension = 6)
print("The filtration of the Rips complex was created in %s" % (time.time() - start_time))
~~~
{: .language-python}

Create the persistence of simplices.
~~~
start_time = time.time()
persistence = simplexTree.persistence()         #Parametros : homology_coeff_field = 11 default, min_persistence , persistence_dim_max
print("The persistente diagram of the Rips complex was created in %s" % (time.time() - start_time))
~~~
{: .language-python}

Print the birth time of the simplices.

~~~
result_str = 'Rips complex of dimension ' + repr(simplexTree.dimension())
print(result_str)
fmt = '%s -> %.2f'
for filtered_value in simplexTree.get_filtration():
    print(tuple(filtered_value))
~~~
{: .language-python}

~~~
simplexTree.dimension(), simplexTree.num_vertices(), simplexTree.num_simplices()
~~~
{: .language-python}

~~~
start_time = time.time()
gudhi.plot_persistence_barcode(
    persistence = persistence, 
    alpha = 0.5,
    colormap = cm.Set2.colors
)
print("Bar code diagram was created in %s" % (time.time() - start_time))
~~~
{: .language-python}

Function for the dimension of the simplices.

~~~
def dimension(list):
    return (len(list[0])-1, list[1])
~~~
{: .language-python}


We filter according to the dimension function: it orders us from largest dimension to smallest and then from longest birth time to smallest.

~~~
all_simplex_sorted_dim_1 = sorted(simplexTree.get_filtration(), key = dimension, reverse = True)
all_simplex_sorted_dim_1 
~~~
{: .language-python}

Obtain the persistence of each simplex.

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


~~~
simplices = list()
simplices = list(d_simplex_const.keys())
~~~
{: .language-python}


~~~
pre_genes_numbers = []
        
for simplex in simplices:
    for i in range(len(simplex)):
        j = 0
        genString = ''
        while simplex[i][j] != '.':
            j = j+1
        if simplex[i][j] == '.':
            j = j+1
            while simplex[i][j] != '.':
                genString = genString + simplex[i][j]
                j = j+1
            pre_genes_numbers.append(genString)
pre_genes_numbers
genes_numbers = []
for gen in pre_genes_numbers:
    if gen not in genes_numbers:
        genes_numbers.append(gen)
genes_numbers
~~~
{: .language-python}

~~~
genes_numbers=genomas
~~~
{: .language-python}


~~~
bool_gen = dict()
genes_contains = dict()
num_new_columns = len(genes_numbers)
for simplex in simplices:
    genes_contains = dict()
    for i in range(len(simplex)):
        for gen in genes_numbers:
            if gen in simplex[i]:
                genes_contains[gen] = 1
    for gen in genes_numbers:
        if gen not in genes_contains.keys():
            genes_contains[gen] = 0
    bool_gen[simplex] = genes_contains
~~~
{: .language-python}


~~~
births = []
deaths = []
persistent_times = []
for values in d_simplex_time.values():
    births.append(values[0])
    deaths.append(values[1])
    persistent_times.append(values[1]-values[0])
~~~
{: .language-python}


~~~
data = {
    't_birth': births,
    't_death': deaths,
    'persistence': persistent_times
}
simplex_list = pd.DataFrame(index = simplices, data = data)
#simplex_list.head(10)
simplex_list.to_csv('~/pan_workshop/results/blast/mini/persistent_simplices.csv')
~~~
{: .language-python}


~~~
for gen in genes_numbers:
    data = dict()
    dataFrame_aux = []
    for simplex in simplices:
        data[simplex] = bool_gen[simplex][gen]
    dataFrame_aux = pd.DataFrame.from_dict(data, orient='index', columns = [str(gen)])
    pd.concat([simplex_list, dataFrame_aux], axis = 1)
    dataFrame_aux.to_csv('~/pan_workshop/results/blast/mini/persistent_simplices_'+str(gen)+'.csv')
~~~
{: .language-python}

Joint the 5 data frames.


This code attempts to generate a pangenome with an specified percentage of core/shell/cloud gene families. Neverthelles caiution is recommended because in a real pangenome, data should have more hierarchical structure. This is temporary data to exemplify TDA concepts.
~~~
import pandas as pd
import numpy as np
import random

import warnings

warnings.filterwarnings("ignore")  # Deactivate Python warnings

m =5 # The number of genomes (Cols)
n = 10  # Specify the number of genes (rows)
shellsize=int(m/2)+1

core = 0.3  # Probability of being core filling a row with ones
shell = 0.6  # Probability of being shell filling the remaining zeros with ones
cloud =1-core-shell

# Create an empty DataFrame
# Create an empty DataFrame
df = pd.DataFrame(columns=[f'Column{i+1}' for i in range(m)])


# Fill each row with ones based on probabilities P1 and P2
for _ in range(n):
    P1=random.random()  # Generate a random value between 0 and 1

    #print("P1",P1)
    if (core>1-P1): 
      row = [1] * m
      #print("core",row)
  
    elif(core+shell>1-P1): 
        #If the row is not already filled with ones, randomly distribute m/2 ones based on P2
        genome_row = np.random.choice([shellsize, m-1])  # Fill the row with ones based on P1
        #print("number of genomes in which this family is present",genome_row)
        # Create an empty list representing the row
        row = [0] * m
        # Generate k unique random indices
        indices = random.sample(range(m), genome_row)
        # Set the values at the random indices to 1
        for index in indices:
          row[index] = 1
        #print("shell",row)
    else:
        row = [0] * m
        num_ones = 1  # Calculate the number of ones to distribute
        # Choose a random index within the row
        random_index = random.randint(0, m - 1)
        #print("rNDOM INDEX", random_index)  
        # Set the value at the random index to 1
        row[random_index] = 1
        #print("cloud",row)

    # row van a ser las columnas del dataframe
    df = df.append(pd.Series(row, index=df.columns), ignore_index=True)  # Append the row to the DataFrame

print(df)
~~~
{: .python}
