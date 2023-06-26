---
title: "Computational Tools for TDA"
teaching: 30
exercises: 15
questions:
- "How can I computationally manipulate simplex"
objectives:
- "Operate simplex in a computational environment"
keypoints:
- "Ghudi is a TDA tool"
---
## Introduction to GUDHI and Simplicial Homology

Welcome to this lesson on using GUDHI and exploring simplicial homology. GUDHI (Geometry Understanding in Higher Dimensions) is an open-source that provides algorithms and data structures for the analysis of geometric data. It offers a wide range of tools for topological data analysis, including simplicial complexes and computations of their homology.  

To begin, we will import the necessary packages.
~~~
from IPython.display import Image
from os import chdir
import numpy as np
import gudhi as gd
import matplotlib.pyplot as plt
import argparse
import seaborn as sns
~~~
{: .language-python}


### **Example 1:** SimplexTree and Manual Filtration
The SimplexTree data structure in GUDHI allows efficient manipulation of simplicial complexes. You can demonstrate its usage by creating a SimplexTree object, adding simplices manually, and then filtering the complex based on a filtration value.   

#### **Create SimplexTree**
With the following command, you can create a SimplexTree object named `st`, which we will use to add the information of your filtered simplicial complex:

~~~
st = gd.SimplexTree()  ## 
~~~
{: .language-python}
#### **Insert simplex**
In GUDHI, you can use the `st.insert()` function to add simplices to a SimplexTree data structure. Additionally, you have the flexibility to specify the filtration level of each simplex. If no filtration level is specified, it is assumed to be added at filtration time 0.

~~~
#insert 0-simplex (the vertex), 
st.insert([0])
st.insert([1])

~~~
{: .language-python}
~~~
True
~~~
{: .output}

Now let's insert 1-simplices at different filtration levels. If adding a simplex requires a lower-dimensional simplex to be present, the missing simplices will be automatically completed.

Here's an example of inserting 1-simplices into the SimplexTree at various filtration levels:

~~~
# Insert 1-simplices at different filtration levels
st.insert([0, 1], filtration=0.5)
st.insert([1, 2], filtration=0.8)
st.insert([0, 2], filtration=1.2)

~~~
{: .language-python}

~~~
True
~~~
{: .output}

In the code snippet above, we created a SimplexTree object named `st`. We then inserted three 1-simplices into the SimplexTree at filtration levels of 0.5, 0.8, and 1.2, respectively. The 1-simplices are defined by specifying their vertices as lists [v1, v2], where v1 and v2 are the indices of the vertices.

> ## Note
> If any lower-dimensional simplices are missing, GUDHI's SimplexTree will automatically complete them. For example, when inserting
> the 1-simplex [1, 2] at filtration level 0.8, if the 0-simplex [1] or [2] was not already present, GUDHI will add it to the SimplexTree.
{: .callout}

This approach allows you to gradually build the simplicial complex by inserting simplices at different filtration levels, and GUDHI takes care of maintaining the necessary lower-dimensional simplices.
Remember to provide clear instructions and explanations for each step and encourage learners to experiment with different filtration levels and simplex insertions to gain a better understanding of how GUDHI's SimplexTree handles complex construction and completion.
> ## FIXME
> Esta última indicación creo que más bien iría en Instructor Notes, no en el episodio.
{: .caution}

Now, you can use the `st.num_vertices()` and `st.num_simplices()` commands to see the number of vertices and simplices, respectively, in your simplicial complex stored in the st SimplexTree object.

~~~
num_vertices = st.num_vertices()
num_simplices = st.num_simplices()

print("Number of vertices:", num_vertices)
print("Number of simplices:", num_simplices)
~~~
{: .language-python}

~~~
Number of vertices: 3
Number of simplices: 6
~~~
{: .output}

The `st.persistence()` function in GUDHI's SimplexTree is used to compute the persistence diagram of the simplicial complex. The persistence diagram provides a compact representation of the birth and death of topological features as the filtration parameter varies.

Here's an example of how to use `st.persistence()`:

~~~
# Compute the persistence diagram
persistence_diagram = st.persistence()

# Print the persistence diagram
for point in persistence_diagram:
    birth = point[0]
    death = point[1]
    print("Birth:", birth)
    print("Death:", death)
    print()
~~~
{: .language-python}

~~~
Birth: 0
Death: (0.0, inf)

Birth: 0
Death: (0.0, 0.5)
~~~
{: .output}

Plot the persistence diagram

~~~
gd.plot_persistence_diagram(persistence_diagram,legend=True)
~~~
{: .language-python}

~~~
<AxesSubplot:title={'center':'Persistence diagram'}, xlabel='Birth', ylabel='Death'>
~~~
{: .output}

> ## FIXME
> Aquí me sale un plot diferente con menos rayas horizontelaes y sin las unidades en el eje Y:
{: .caution}

 <a href="../fig/tda_09_diagram_1.png">
  <img src="../fig/tda_09_diagram_1.png" alt="Persistence Diagram" width="50%" height="auto" />
</a>


Plot the barcode

~~~
gd.plot_persistence_barcode(persistence_diagram,legend=True)
~~~
{: .language-python}

~~~
<AxesSubplot:title={'center':'Persistence barcode'}>
~~~
{: .output}
 <a href="../fig/tda_09_barcode_1.png">
  <img src="../fig/tda_09_barcode_1.png" alt="Persistence Diagram" width="50%" height="auto" />
</a>

> ## FIXME
> En el ejercicio 1 el la K se ve con los signos y no como Latex. Hay que poner el código completo en la solución del ejercicio, en lugar del texto
> que dice qué funciones usar. Con el código que yo pude hacer no me salió la misma gráfica.
{: .caution}

> ## Exercise 1: Creating a Manually Filtered Simplicial Complex.
>  In the following graph, we have $K$ a simplicial complex filtered representation of simplicial complexes.
>  <a href="../fig/tda_09_filtracion_ex.png">
  <img src="../fig/tda_09_filtracion_ex.png" alt="Exercise 1 Filtration" width="100%" height="auto"/>
</a>
> 
> Perform persistent homology and plot the persistence diagram and barcode.
> > ## Solution  
> > 1. Create a SimplexTree with `gd.SimplexTree()`.
> >  ~~~
>> st = gd.SimplexTree()  
>> ~~~
>> {: .language-python}  
>> 2. Insert vertices at time 0 using `st.insert()`
> >  ~~~
>> #insert 0-simplex (the vertex), 
>> st.insert([0])
>>  st.insert([1])
>> st.insert([2])
>>  st.insert([3])
>>  st.insert([4])
>> ~~~
>> {: .language-python}   
>> 3. Insert the remaining simplices by setting the filtration time using `st.insert([0, 1], filtration=)`.
> >  ~~~
>> #insert 1-simplex level filtration 1 
>> st.insert([0, 2], filtration=1)
>> st.insert([3, 4], filtration=1)
>> #insert 1-simplex level filtration 2 
>> st.insert([0, 1], filtration=2)
>> #insert 1-simplex level filtration 3 
>>  st.insert([2, 1], filtration=3)
>> #insert 1-simplex level filtration 4 
>> st.insert([2, 1,0], filtration=4)
>> ~~~
>> {: .language-python}   
>> 4. Perform persistent homology using `st.persistence()`.
> >  ~~~
>># Compute the persistence diagram
>> persistence_diagram = st.persistence() 
>> ~~~
>> {: .language-python}  
>> 5. Plot the persistence diagram.
> >  ~~~
>># plot the persistence diagram
>> gd.plot_persistence_diagram(persistence_diagram,legend=True)
>> ~~~
>> {: .language-python}  
>> 6. Plot the barcode.
> >  ~~~
>> gd.plot_persistence_barcode(persistence_diagram,legend=True)
>> ~~~
>> {: .language-python}  
> >6. Get this output  
>> <a href="../fig/tda_09_diagram_F.png">
>>   <img src="../fig/tda_09_diagram_F.png" alt="Persistence Diagram" width="50%" height="auto" />
>> </a>  
> {: .solution}

{: .challenge}



### **Example 2:** Rips complex from datasets 
In this example, we will demonstrate an application of persistent homology using a dataset generated by us. Persistent homology is a technique used in topological data analysis to study the shape and structure of datasets.

The `make_circles()` function from scikit-learn's datasets module is used to generate synthetic circular data. We specify the number of points to generate (n_samples), the amount of noise to add to the data points (noise), and the scale factor between the inner and outer circle (factor).

The generated dataset consists of two arrays: circles and labels. The circles array contains the coordinates of the generated data points, while the labels array assigns a label to each data point (in this case, it will be 0 or 1 representing the two circles).

> ## FIXME
> Aquí falta más texto para decir que vas a generar una nube de puntos y qué es make_circles n_samples, noise y factor. Y en general en los siguientes pasos falta un poco de texto y decir qué hacen las funciones nuevas.
{: .caution}

~~~
from sklearn import datasets  # Import the datasets module from scikit-learn

# Generate synthetic data using the make_circles function
# n_samples: Number of points to generate
# noise: Standard deviation of Gaussian noise added to the data
# factor: Scale factor between inner and outer circle
circles, labels = datasets.make_circles(n_samples=100, noise=0.06, factor=0.5)

# Print the dimensions of the generated data
print('Data dimension: {}'.format(circles.shape))
~~~
{: .language-python}

~~~
Data dimension:(100, 2)
~~~
{: .output}

Plot dataset
> ## FIXME
> Aquí falta saber qué hace sns.set() y decir que vas a hacer un scatterplot que represente la nube de puntos
{: .caution}

~~~
fig = plt.figure()  # Create a new figure
ax = fig.add_subplot(111)  # Add a subplot to the figure
ax = sns.scatterplot(x=circles[:,0], y=circles[:,1], s=15)  # Create a scatter plot using seaborn
plt.title('Dataset with N=%s points'%(circles.shape[0]))  # Set the title of the plot
plt.show()  # Display the plot
~~~
{: .language-python}
 <a href="../fig/tda_circles.png">
  <img src="../fig/tda_circles.png" alt="Plot Circles" width="50%" height="auto" />
</a>

 The `RipsComplex()` function creates a one-skeleton graph from the point cloud.

> ## FIXME
> A partir de aquí el código no me funciona
{: .caution}

~~~
%%time
Rips_complex = gd.RipsComplex(circles, max_edge_length=0.6) 
~~~
{: .language-python}
~~~
CPU times: user 461 µs, sys: 88 µs, total: 549 µs
Wall time: 557 µs
~~~
{: .output}

The `create_simplex_tree()` method creates the filtered complex.
~~~
%%time

Rips_simplex_tree = Rips_complex.create_simplex_tree(max_dimension=3) 
~~~
{: .language-python}

~~~
CPU times: user 2.25 ms, sys: 0 ns, total: 2.25 ms
Wall time: 1.95 ms
~~~
{: .output}

The `get_filtration()` method computes the simplices of the filtration
~~~
%%time

filt_Rips = list(Rips_simplex_tree.get_filtration())
~~~
{: .language-python}
~~~
CPU times: user 108 ms, sys: 0 ns, total: 108 ms
Wall time: 121 ms
~~~
{: .output}

We can compute persistence on the simplex tree structure using the `persistence()` method
~~~
%%time

diag_Rips = Rips_simplex_tree.persistence()
~~~
{: .language-python}
~~~
CPU times: user 9.13 ms, sys: 41 µs, total: 9.17 ms
Wall time: 8.82 ms
~~~
{: .output}


~~~
%%time
gd.plot_persistence_diagram(diag_Rips,legend=True)
~~~
{: .language-python}
~~~
(array([-0.1       ,  0.        ,  0.1       ,  0.2       ,  0.3       ,
         0.4       ,  0.5       ,  0.62569893]),
 [Text(0, -0.1, '-0.100'),
  Text(0, 0.0, '0.000'),
  Text(0, 0.1, '0.100'),
  Text(0, 0.20000000000000004, '0.200'),
  Text(0, 0.30000000000000004, '0.300'),
  Text(0, 0.4, '0.400'),
  Text(0, 0.5000000000000001, '0.500'),
  Text(0, 0.6256989291775961, '$+\\infty$')])
~~~
{: .output}
 <a href="../fig/tda_09_persistence_example2.png">
  <img src="../fig/tda_09_persistence_example2.png" alt="Persistence diagram" width="50%" height="auto" />
</a>



~~~
%%time
gd.plot_persistence_barcode(diag_Rips,legend=True)
plt.grid(color = 'black', linestyle = '-', linewidth = 1)
plt.savefig('persistencebarcodeCircles' , dpi=600, transparent=True)
plt.xticks(size=15)
plt.yticks(size=15)
~~~
{: .language-python}
~~~
(array([  0.,  20.,  40.,  60.,  80., 100., 120.]),
 [Text(0, 0.0, '0'),
  Text(0, 20.0, '20'),
  Text(0, 40.0, '40'),
  Text(0, 60.0, '60'),
  Text(0, 80.0, '80'),
  Text(0, 100.0, '100'),
  Text(0, 120.0, '120')])
~~~
{: .output}
 <a href="../fig/tda_09_bardcode_example2.png">
  <img src="../fig/tda_09_bardcode_example2.png" alt="Bard Code" width="50%" height="auto" />
</a>



### **Example 3:** Rips complex from datasets 

> ## FIXME
> Aquí falta texto para describir la próxima actividad y luego texto en los pasos.
> Quitar o traducir los comentarios que están en el código en español
{: .caution}

~~~
from gudhi.datasets.generators import _points
from gudhi import AlphaComplex
~~~
{: .language-python}


~~~
import requests
#load the file spiral_2d.csv
url = 'https://raw.githubusercontent.com/paumayell/pangenomics/gh-pages/files/spiral_2d.csv'
# Obtener el contenido del archivo
response = requests.get(url)
content = response.text
# Cargar los datos en un arreglo de NumPy
data = np.loadtxt(content.splitlines(), delimiter=' ')
# Graficar los puntos
plt.scatter(data[:, 0], data[:, 1], marker='.', s=1)
plt.show()
~~~
{: .language-python}
<a href="../fig/tda_09_sperial.png">
  <img src="../fig/tda_09_sperial.png" alt="Plot Spiral" width="50%" height="auto" />
</a>

Define simplicial complex
~~~
alpha_complex = AlphaComplex(points = data)
simplex_tree = alpha_complex.create_simplex_tree()
~~~
{: .language-python}

> ## FIXME
> Aquí me sale un mensaje de warning grande, pero sí me sale la misma gráfica
{: .caution}

~~~
diag = simplex_tree.persistence()
diag = simplex_tree.persistence(homology_coeff_field=2, min_persistence=0)
print("diag=", diag)

gd.plot_persistence_diagram(diag)
~~~
{: .language-python}
<a href="../fig/tda_09_persistence_example3.png">
  <img src="../fig/tda_09_persistence_example3.png" alt="Persistence diagram" width="50%" height="auto" />
</a>

> ## FIXME
> Aquí hay que decir algo del código que está comentado o quitarlo
{: .caution}
~~~
gd.plot_persistence_barcode(diag)
#plt.savefig('persistence_barcodeSpiral.svg' , dpi=1200)
plt.show()
~~~
{: .language-python}
 <a href="../fig/tda_09_bardcode_example3.png">
  <img src="../fig/tda_09_bardcode_example3.png" alt="Bard Code" width="50%" height="auto" />
</a>

> ## FIXME
> Aquí también decir algo o quitar el código en el que se guarda la imagen. Aquí no está la imagen que se genera, que es igual a la anterior pero con > el formato un poquito diferente. Decir por qué se está haciendo.
{: .caution}

~~~
%%time
gd.plot_persistence_barcode(diag,legend=True)
plt.grid(color = 'black', linestyle = '-', linewidth = 1)
#plt.savefig('persistencebarcodeCircles' , dpi=600, transparent=True)
plt.xticks(size=15)
plt.yticks(size=15)
~~~
{: .language-python}

> ## FIXME
> En el ejercicio poner el código completo para llegar a la solución
{: .caution}

> ## Exercise 2: Torus.
>  To build a torus using the tadasets function and apply persistent homology.
> <a href="../fig/tda_09_torus.png">
  <img src="../fig/tda_09_torus.png" alt="Exercise 2 Torus" width="50%" height="auto"/>
</a>
> > ## Solution  
> > 1. `import tadasets`.
>> 2. `torus = tadasets.torus(n=100)`
>> 3. Create a Rips complex from the torus points `gd.RipsComplex(points=torus)`
>> 4. Obtain the simplicial complex `rips_complex.create_simplex_tree(max_dimension=2)`
>> 5. Compute the persistent homology of the simplicial complex `simplicial_complex.persistence()`
>> 6. Plots diagrams
> {: .solution}
{: .challenge}

> ## FIXME
> Poner algo más en los keypoints
{: .caution}

