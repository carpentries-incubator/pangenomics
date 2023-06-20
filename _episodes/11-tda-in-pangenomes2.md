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

### **1. Library**
To begin, we will import the necessary packages.
~~~
import pandas as pd
import numpy as np
import umap
import umap.plot
import gudhi as gd
import matplotlib.pyplot as plt
#from umap import UMAPTransformer
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
Calculate 
~~~
distancia(df_libro)
~~~
{: .language-python}

~~~
array([[0. , 0.5, 0.5, 1. ],
       [0.5, 0. , 1. , 0.5],
       [0.5, 1. , 0. , 0.5],
       [1. , 0.5, 0.5, 0. ]])
~~~
{: .output}

~~~
matrix_dintancia_libro=distancia(df_libro)
persistence_libro=complejo(matrix_dintancia_libro)
~~~
{: .language-python}


~~~
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

