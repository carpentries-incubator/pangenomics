---
title: "Examples TDA in genomics "
teaching: 30
exercises: 15
questions:
- "How can I apply TDA to describe Pangenomes"
objectives:
- "Describe Pangenomes using Gudhi"
keypoints:
- "Pangenomes can be described using TDA"
---

### **1. Library**
To begin, we will import the necessary packages.
~~~
import pandas as pd
import numpy as np
import umap
import gudhi as gd
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from itertools import combinations
import scipy.cluster.hierarchy as sch

from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo import draw
~~~
{: .language-python}

### Funtions
We  define funtions


~~~
##Calcula la secuencia media a partiar de secuencias
def compute_median_sequence(a, b, c):
    median = ""
    for i in range(len(a)):
        counts = {'0': 0, '1': 0}
        counts[a[i]] += 1
        counts[b[i]] += 1
        counts[c[i]] += 1
        majority = max(counts, key=counts.get)
        median += majority
    return median
##Calcula la secuencia media a partiar de vectores
def compute_median_vector(a, b, c):
    a = np.array(a)
    b = np.array(b)
    c = np.array(c)
    median = np.zeros_like(a)
    for i in range(len(a)):
        counts = np.bincount([a[i], b[i], c[i]])
        majority = np.argmax(counts)
        median[i] = majority
    return median.tolist()
##Crea un nuevo diccionario, a partir del anterior, y le agrega los puntos medios
def process_dict_elements(dictionary):
    keys = list(dictionary.keys())
    result = {}
    for combination in combinations(keys, 3):
        a, b, c = combination
        median = compute_median_vector(dictionary[a], dictionary[b], dictionary[c])
        new_key = f"{a}_{b}_{c}"  # Clave nueva basada en la combinación de claves originales
        result[new_key] = median
    dictionary.update(result)  # Agregar las medianas calculadas al diccionario original
    return dictionary
##calcula la matriz de distancia a partir de un dataframe
def distancia(df,metrica='hamming'):
    distances = pdist(df.values.T, metric=metrica)
    distance_matrix = squareform(distances)
    return(distance_matrix)
##Calcula el comprejo de Rips y regresa la persistencia.
def complejo(distance_matrix):
# Create the simplicial complex from the distance matrix
    rips_complex = gd.RipsComplex(distance_matrix)
    simplex_tree = rips_complex.create_simplex_tree(max_dimension=3)
# Compute persistence
    persistence = simplex_tree.persistence()
    return(persistence)


def hamming_distance(string1: str, string2: str) -> int:
    """Return the Hamming distance between two strings."""
    if len(string1) != len(string2):
        raise ValueException("Strings must be of equal length.")
    dist_counter = 0
    for n in range(len(string1)):
        if string1[n] != string2[n]:
            dist_counter += 1
    return dist_counter

def plot_cladogram(data):
    # Convert the dictionary into a matrix
    matrix = np.array([list(value) for value in data.values()])

    # Calculate the distance matrix
    dist_matrix = sch.distance.pdist(matrix)

    # Perform hierarchical clustering
    linkage_matrix = sch.linkage(dist_matrix)

    # Plot the dendrogram
    plt.figure(figsize=(8, 6))
    dendrogram = sch.dendrogram(linkage_matrix, labels=list(data.keys()), orientation='right')

    # Adjust the margins and labels of the x-axis
    plt.subplots_adjust(bottom=0.1)
    plt.xticks(rotation='vertical')

    # Show the cladogram
    plt.show()

def crear_archivo_txt(diccionario):
    with open('genomas.fa', 'w') as archivo:
        for clave, valores in diccionario.items():
            archivo.write(">" + clave + "\n")  # Escribe la clave con ">" al inicio
            cadena_valores = ''.join(str(valor) for valor in valores)  # Concatena los elementos de la lista
            archivo.write(cadena_valores + "\n")  # Escribe la cadena de valores    
~~~
{: .language-python}

### Example 1: From book
~~~
data = {'Genoma1': [0, 0],
        'Genoma2': [1, 0],
        'Genoma3': [0, 1],
        'Genoma4': [1, 1]}
df_libro = pd.DataFrame(data, index=['Gen1', 'Gen2'])
df_libro
~~~
{: .language-python}

~~~
Genoma1	Genoma2	Genoma3	Genoma4
Gen1	0	1	0	1
Gen2	0	0	1	1
~~~
{: .output}

The distancia function takes a DataFrame df and an optional parameter metrica (defaulting to 'hamming').

The resulting distance_matrix can be used for further analysis, such as clustering or dimensionality reduction, to explore the relationships and similarities between the variables/columns of the DataFrame. 
~~~
def distancia(df, metrica='hamming'):
    # Compute pairwise distances between columns of the DataFrame
    distances = pdist(df.values.T, metric=metrica)
    
    # Convert the condensed distance matrix to a squareform distance matrix
    distance_matrix = squareform(distances)

~~~
{: .language-python}
Calculate 
~~~
matrix_dintancia_libro=distancia(df_libro)
~~~
{: .language-python}

~~~
array([[0. , 0.5, 0.5, 1. ],
       [0.5, 0. , 1. , 0.5],
       [0.5, 1. , 0. , 0.5],
       [1. , 0.5, 0.5, 0. ]])
~~~
{: .output}

Define the function complejo to compute the persistence of a Rips simplicial complex from a distance matrix.
~~~
def complejo(distance_matrix):
    # Create the Rips simplicial complex from the distance matrix
    rips_complex = gd.RipsComplex(distance_matrix)
    # Create the simplex tree from the Rips complex with a maximum dimension of 3
    simplex_tree = rips_complex.create_simplex_tree(max_dimension=3)
    # Compute the persistence of the simplicial complex
    persistence = simplex_tree.persistence()
    # Return the persistence diagram or barcode
    return persistence
~~~
{: .language-python}

we used the previosuly function and calcultate de persitence and plot
~~~
persistence_libro=complejo(matrix_dintancia_libro)
gd.plot_persistence_barcode(persistence_libro)
~~~
{: .language-python}
<a href="../fig/tda_11_barcode_1.png">
  <img src="../fig/tda_11_barcode_1.png" alt="Persistence Barcode" width="50%" height="auto" />
  </a>
~~~
gd.plot_persistence_diagram(persistence_libro,legend=True)
~~~
{: .language-python}
<a href="../fig/tda_11_diagram_1.png">
  <img src="../fig/tda_11_diagram_1.png" alt="Persistence Diagram" width="50%" height="auto" />
  </a>




### Example 2: Genes
In this example, we want to use Topological Data Analysis to detect if there is Horizontal Gene Transfer within this group of genomes.

First, let's import the file `familias_mini.csv` which contains a table of gene presence and absence in 4 _Streptococcus_ genomes.
~~~
df = pd.read_csv("/home/shaday/GIT/pangenomics/files/familias_minis.csv", index_col=0)
df_filled = df.fillna(0)
df=df_filled.replace(to_replace=r'.+', value=1, regex=True)
df
~~~
{: .language-python}
~~~
	g_A909	g_2603V	g_515	g_NEM316
A909|MGIDGNCP_01408	1	1	1	1
A909|MGIDGNCP_00096	1	1	1	1
A909|MGIDGNCP_01343	1	0	0	1
A909|MGIDGNCP_01221	1	0	1	0
A909|MGIDGNCP_01268	1	1	1	1
A909|MGIDGNCP_00580	1	1	1	1
A909|MGIDGNCP_00352	1	1	1	1
A909|MGIDGNCP_00064	1	1	1	1
A909|MGIDGNCP_00627	1	0	0	0
A909|MGIDGNCP_01082	1	1	0	0
A909|MGIDGNCP_00877	1	1	1	1
A909|MGIDGNCP_00405	1	1	1	1
2603V|GBPINHCM_00748	0	1	0	0
2603V|GBPINHCM_01226	0	1	0	0
515|LHMFJANI_01625	0	0	1	1
~~~
{: .output}

Now we will use the file 'minigenomes_allig.fasta,' which contains the sequence of the previously aligned genomes, to build a phylogenetic tree among them.
~~~
url = "https://raw.githubusercontent.com/paumayell/pangenomics/gh-pages/files/minigenomes_allig.fasta"
response = requests.get(url)
response.raise_for_status()  # Check if any errors occurred during the download
# Save the downloaded content to a local file
with open("minigenomes_allig.fasta", "wb") as file:
    file.write(response.content)
sequences = list(SeqIO.parse("minigenomes_allig.fasta", "fasta"))
# Rest of your code that uses the sequences
alignment = MultipleSeqAlignment(sequences)
# Calculate the distance matrix
calculator = DistanceCalculator('identity')
distance_matrix = calculator.get_distance(alignment)

# Build the UPGMA tree
constructor = DistanceTreeConstructor()
upgma_tree = constructor.upgma(distance_matrix)

# Draw the UPGMA tree
draw(upgma_tree)

~~~
{: .language-python}
 <a href="../fig/tda_11_philo_tree.png">
  <img src="../fig/tda_11_philo_tree.png" alt="Phylogenetic tree" width="50%" height="auto" />
</a>

The phylogenetic tree groups the genomes into pairs, which does not help infer whether horizontal gene transfer occurred at some point during evolution among these species. Next, we will use persistent homology to try to detect this by identifying 1-hole structures.


~~~
matrix_dintancia_genes=distancia(df)
persistence_genes=complejo(matrix_dintancia_genes)
gd.plot_persistence_barcode(persistence_genes)
~~~
{: .language-python}
 <a href="../fig/tda_11_barcode_5.png">
  <img src="../fig/tda_11_barcode_5.png" alt="Persistence Barcode" width="50%" height="auto" />
</a>

In the persistence barcode code, we did not detect any 1-hole structures. We can explore various strategies to try to detect this.

### Select by triplets.
We start select the first four genes and repeat the previous calculations.
~~~
df_primera=df.iloc[:4,:]
matrix_dintancia_genes_primera=distancia(df_primera)
persistence_genes_primera=complejo(matrix_dintancia_genes_primera)
gd.plot_persistence_barcode(persistence_genes_primera)
~~~
{: .language-python}
 <a href="../fig/tda_11_barcode_4.png">
  <img src="../fig/tda_11_barcode_4.png" alt="Persistence Barcode" width="50%" height="auto" />
</a>
We can observe that we have a 1-hole, indicating the presence of horizontal gene transfer among our four genomes.

Now, if we make another selection by taking the last four genes, we can observe the following:

~~~
df_segunda=df.iloc[-4:,:]
df_segunda
~~~
{: .language-python}
~~~
	g_A909	g_2603V	g_515	g_NEM316
A909|MGIDGNCP_00405	1	1	1	1
2603V|GBPINHCM_00748	0	1	0	0
2603V|GBPINHCM_01226	0	1	0	0
515|LHMFJANI_01625	0	0	1	1
~~~
{: .output}
~~~
matrix_dintancia_genes_segunda=distancia(df_segunda)
persistence_genes_segunda=complejo(matrix_dintancia_genes_segunda)
gd.plot_persistence_barcode(persistence_genes_segunda)
~~~
{: .language-python}
 <a href="../fig/tda_11_barcode_5.png">
  <img src="../fig/tda_11_barcode_5.png" alt="Persistence Barcode" width="50%" height="auto" />
</a>

In this selection, we do not detect this 1-hole.

### The mediam complex
~~~
#crear un diccionario de cada genoma convertivo a "0" y "1" de presencia y auscnia de genes
genomas = {}
for columna in df.columns:
    genomas[columna] = list(np.array(df[columna]))
genomas
~~~
{: .language-python}
~~~
{'g_A909': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0],
 'g_2603V': [1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0],
 'g_515': [1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1],
 'g_NEM316': [1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1]}
~~~
{: .output}

~~~
genomas_mediam=process_dict_elements(genomas)
genomas_mediam
~~~
{: .language-python}
~~~
{'g_A909': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0],
 'g_2603V': [1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0],
 'g_515': [1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1],
 'g_NEM316': [1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1],
 'g_A909_g_2603V_g_515': [1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0],
 'g_A909_g_2603V_g_NEM316': [1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0],
 'g_A909_g_515_g_NEM316': [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1],
 'g_2603V_g_515_g_NEM316': [1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1]}
~~~
{: .output}


~~~
df_mediam = pd.DataFrame.from_dict(genomas_mediam)
df_mediam
~~~
{: .language-python}
~~~
	g_A909_g_2603V_g_NEM316	g_A909_g_515_g_NEM316	g_2603V_g_515_g_NEM316
0	1	1	1
1	1	1	1
2	1	1	0
3	0	1	0
~~~
{: .output}

~~~
matrix_dintancia_extendida=distancia(df_mediam)
persistence_extendida=complejo(matrix_dintancia_extendida)
matrix_dintancia_extendida
~~~
{: .language-python}
~~~
array([[0.        , 0.33333333, 0.26666667, 0.26666667, 0.13333333,
        0.13333333, 0.2       , 0.33333333],
       [0.33333333, 0.        , 0.33333333, 0.33333333, 0.2       ,
        0.2       , 0.4       , 0.26666667],
       [0.26666667, 0.33333333, 0.        , 0.13333333, 0.13333333,
        0.26666667, 0.06666667, 0.06666667],
       [0.26666667, 0.33333333, 0.13333333, 0.        , 0.26666667,
        0.13333333, 0.06666667, 0.06666667],
       [0.13333333, 0.2       , 0.13333333, 0.26666667, 0.        ,
        0.13333333, 0.2       , 0.2       ],
       [0.13333333, 0.2       , 0.26666667, 0.13333333, 0.13333333,
        0.        , 0.2       , 0.2       ],
       [0.2       , 0.4       , 0.06666667, 0.06666667, 0.2       ,
        0.2       , 0.        , 0.13333333],
       [0.33333333, 0.26666667, 0.06666667, 0.06666667, 0.2       ,
        0.2       , 0.13333333, 0.        ]])
~~~
{: .output}

~~~
gd.plot_persistence_barcode(persistence_extendida)
~~~
{: .language-python}
 <a href="../fig/tda_11_barcode_3.png">
  <img src="../fig/tda_11_barcode_3.png" alt="Persistence Barcode" width="50%" height="auto" />
</a>

~~~
gd.plot_persistence_diagram(persistence_extendida,legend=True)
~~~
{: .language-python}
 <a href="../fig/tda_11_diagram_3.png">
  <img src="../fig/tda_11_diagram_3.png" alt="Persistence Diagram" width="50%" height="auto" />
</a>

## Method 2

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
