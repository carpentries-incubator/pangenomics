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

We will work with the four mini genomes of episode 4. First we need to import all the libraries that we will use.

~~~
import pandas as pd
from matplotlib import cm
import numpy as np
import gudhi
import time
import os  
~~~
{: .language-python}

Now, we need to read the `mini-genomes.blast` file that we produce in episode 4. 
Read blastp matrices from episode 4.

~~~
os.getcwd()
blastE = pd.read_csv( '~/pan_workshop/results/blast/mini/output-blast/mini-genomes.blast', sep = '\t',names = ['qseqid','sseqid', 'evalue'])  
~~~
{: .language-python}

Obtain a list with the unique genes.
~~~
qseqid_unique=pd.unique(blastE['qseqid'])
sseqid_unique=pd.unique(blastE['sseqid'])
genes = pd.unique(np.append(qseqid_unique, sseqid_unique))
~~~
{: .language-python}

We have 43 unique genes, we can check it as follows.

~~~
len(genes)
~~~
{: .language-python}


To use the `gudhi` packages, we need a distance matrix. In this case we will use the `evalue` as the mesure of how similar the genes are. First, we will process the `blastE` data frame to a list and then we will convert in a matrix object.

~~~
distance_list = blastE[ blastE['qseqid'].isin(genes) & blastE['sseqid'].isin(genes)]
distance_list.head()
~~~
{: .language-python}


~~~
  qseqid	              sseqid	              evalue
0	2603V|GBPINHCM_01420	NEM316|AOGPFIKH_01528	4.110000e-67
1	2603V|GBPINHCM_01420	A909|MGIDGNCP_01408	4.110000e-67
2	2603V|GBPINHCM_01420	515|LHMFJANI_01310	4.110000e-67
3	2603V|GBPINHCM_01420	2603V|GBPINHCM_01420	4.110000e-67
4	2603V|GBPINHCM_01420	A909|MGIDGNCP_01082	1.600000e+00
~~~
{: .output}


To convert the `distance_list` to a matrix object we will use the convention that the maximum biological distance between genes are `5`, so if we do not have the `evalue` between two genes, that imply that the evalua was too big, so we will fill that spots with the maximum biological distance.

~~~
MaxDistance = 5.0000000

# reshape long to wide
matrixE = pd.pivot_table(distance_list,index = "qseqid",values = "evalue",columns = 'sseqid')
matrixE.iloc[1:5,1:5]
~~~
{: .language-python}

~~~
sseqid	2603V|GBPINHCM_00065	2603V|GBPINHCM_00097	2603V|GBPINHCM_00348	2603V|GBPINHCM_00401
qseqid				
2603V|GBPINHCM_00065	1.240000e-174	NaN	NaN	NaN
2603V|GBPINHCM_00097	NaN	9.580000e-100	NaN	NaN
2603V|GBPINHCM_00348	NaN	NaN	0.0	NaN
2603V|GBPINHCM_00401	NaN	NaN	NaN	2.560000e-135
~~~
{: .output}


~~~~
matrixE2=matrixE.fillna(MaxDistance)
matrixE2.iloc[0:4,0:4]
~~~~
{: .language-python}

~~~
sseqid	2603V|GBPINHCM_00065	2603V|GBPINHCM_00097	2603V|GBPINHCM_00348	2603V|GBPINHCM_00401
qseqid				
2603V|GBPINHCM_00065	1.240000e-174	5.000000e+00	5.0	5.000000e+00
2603V|GBPINHCM_00097	5.000000e+00	9.580000e-100	5.0	5.000000e+00
2603V|GBPINHCM_00348	5.000000e+00	5.000000e+00	0.0	5.000000e+00
2603V|GBPINHCM_00401	5.000000e+00	5.000000e+00	5.0	2.560000e-135
~~~
{: .output}


Finaly, we need the distance matrix as a `numpy` array.

~~~~
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
        5.00e+000],
       ...,
       [5.00e+000, 5.00e+000, 5.00e+000, ..., 1.64e-143, 5.00e+000,
        5.00e+000],
       [5.00e+000, 5.00e+000, 5.00e+000, ..., 5.00e+000, 4.11e-067,
        5.00e+000],
       [5.00e+000, 5.00e+000, 5.00e+000, ..., 5.00e+000, 5.00e+000,
        0.00e+000]])
~~~
{: .output}


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
