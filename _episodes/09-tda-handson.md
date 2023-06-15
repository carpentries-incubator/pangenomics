---
title: "Computational Tools for TDA"
teaching: 30
exercises: 15
questions:
- "How can I computationaly manipulate simplex"
objectives:
- "Operate simplex in a computational environment"
keypoints:
- "Ghudi is a TDA tool"
---
# **Introduction to GUDHI and Simplicial Homology**

Welcome to this lesson on using GUDHI and exploring simplicial homology. GUDHI (Geometry Understanding in Higher Dimensions) is an open-source C++ library that provides algorithms and data structures for the analysis of geometric data. It offers a wide range of tools for topological data analysis, including simplicial complexes and computations of their homology.
## **1. Library**
~~~
from IPython.display import Image
from os import chdir
import numpy as np
import gudhi as gd
import matplotlib.pyplot as plt
import argparse
~~~
{: .language-python}


## **Example 1:** SimplexTree and Manual Filtration
The SimplexTree data structure in GUDHI allows efficient manipulation of simplicial complexes. You can demonstrate its usage by creating a SimplexTree object, adding simplices manually, and then filtering the complex based on a filtration value.
### **Create SimplexTree**
With the following command, you can create a SimplexTree object named ´st´, which we will use to add the information of your filtered simplicial complex:

~~~
st = gd.SimplexTree()  ## 
~~~
{: .language-python}
### **Insert simplex**
In GUDHI, you can use the 'st.insert()' function to add simplices to a SimplexTree data structure. Additionally, you have the flexibility to specify the filtration level of each simplex. If no filtration level is specified, it is assumed to be added at filtration time 0.

~~~
#insert 0-simplex (the vertex), 
st.insert([0])
st.insert([1])

~~~
{: .language-python}

Now let's insert 1-simplices at different filtration levels. If adding a simplex requires a lower-dimensional simplex to be present, the missing simplices will be automatically completed.

Here's an example of inserting 1-simplices into the SimplexTree at various filtration levels:

~~~
# Insert 1-simplices at different filtration levels
st.insert([0, 1], filtration=0.5)
st.insert([1, 2], filtration=0.8)
st.insert([0, 2], filtration=1.2)

~~~
{: .language-python}

In the code snippet above, we create a SimplexTree object named st. We then insert three 1-simplices into the SimplexTree at filtration levels of 0.5, 0.8, and 1.2, respectively. The 1-simplices are defined by specifying their vertices as lists [v1, v2], where v1 and v2 are the indices of the vertices.

> Note: If any lower-dimensional simplices are missing, GUDHI's SimplexTree will automatically complete them. For example, when inserting the 1-simplex [1, 2] at filtration level 0.8, if the 0-simplex [1] or [2] was not already present, GUDHI will add it to the SimplexTree.

This approach allows you to gradually build the simplicial complex by inserting simplices at different filtration levels, and GUDHI takes care of maintaining the necessary lower-dimensional simplices.
Remember to provide clear instructions and explanations for each step and encourage learners to experiment with different filtration levels and simplex insertions to gain a better understanding of how GUDHI's SimplexTree handles complex construction and completion.


Now, you can use the ´st.num_vertices()´ and ´st.num_simplices()´ commands to see the number of vertices and simplices, respectively, in your simplicial complex stored in the st SimplexTree object.

~~~
num_vertices = st.num_vertices()
num_simplices = st.num_simplices()

print("Number of vertices:", num_vertices)
print("Number of simplices:", num_simplices)
~~~
{: .language-python}


The st.persistence() function in GUDHI's SimplexTree is used to compute the persistence diagram of the simplicial complex. The persistence diagram provides a compact representation of the birth and death of topological features as the filtration parameter varies.

Here's an example of how to use st.persistence():

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


