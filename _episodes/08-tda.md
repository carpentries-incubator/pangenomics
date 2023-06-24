---
title: "Topological Data Analysis"
teaching: 30
exercises: 15
questions:
- "What is topological data analysis?"
objectives:
- "Understand that simplices comprise vertices, edges, triangles, and higher dimensions structures."
keypoints:
- "TDA describes data forms."
math: true
---

## Topological data analysis

Topological data analysis is a technique that uses concepts from topology to analyze complex data and find patterns and structures that are not apparent at first glance. This technique is based on constructing a simplicial complex composed of a collection of simple geometric objects called simplices. The topology of this complex is used to analyze and visualize the relationships between the data.

### Simplex
A simplex (plural simplices) is a simple geometric object of any dimension (point, line segment, triangle, tetrahedron, etc.). Simplices are used to construct simplicial complexes. In the following figure, some examples of simplices are shown.  

> ## Simplex
>
> > ## Mathematical definition
>> Given a set $ P=\\{p_0,...,p_k\\}\subset \mathbb{R}^d $ of $ k+1 $ affinely independent points, the **k-dimensional simplex** $\sigma$ (or **k-simplex** for short) spanned by $P$ is the set of convex combinatios 
>>
>> $ \sum_{i=0}^k\lambda_ip_i, \quad with \quad  \sum_{i=0}^k\lambda_i = 1 \quad \lambda_i \geq 0. $ 
>>
>> 
>> The points $p_0, ..., p_k$ are called the vertices of $\sigma$.
>> 
> {: .solution}
{: .discussion}

  <a href="../fig/tda_Vertices.png">
  <img src="../fig/tda_Vertices.png" alt="Example" width="70%" height="auto"/>
</a>

### Simplicial complex
A simplicial complex is a mathematical structure comprising a collection of simplices constructed from a data set. In other words, a simplicial complex is a collection of vertices, edges, triangles, tetrahedra, and other elements that follows specific rules. As a starting point, we can think of a simplicial complex as extending the notion of a graph only formed by vertices and edges. The following figure shows an example of a simplicial complex.  

  <a href="../fig/Tda_Simplicial_complex_example.png">
  <img src="../fig/Tda_Simplicial_complex_example.png" alt="Example Simplicial Complex" width="70%" height="auto"/>
</a>

> ## Simplicial complex
>
> > ## Mathematical definition
>> A simplicial complex $K$ in $ \mathbb{R}^d $ is a collection of simplices s.t:
>> 1. Any face of a simplex from $K$ is a simplex of $K$,
>> 2. The intersection of any two simplices of $K$ is either empty or a common face of both.
> > 
> {: .solution}
{: .discussion}

  <a href="../fig/tda_paste.png">
  <img src="../fig/tda_paste.png" alt="Example" width="70%" height="auto"/>
</a>

"In the figure, in a) we have an example of a simplicial complex, while b) does not satisfy the condition that the intersection of any pair of simplices is either empty or a face, in other words, it must fit together properly."


### Abstract simplex and simplicial complex

> ## Abstract Simplicial Complex
>
> > ## Mathematical definition  
>>Let $P= \\{p_1,...,p_n\\}$ be a (finite) set. An **abstract simplicial complex** $K$ with vertex set $P$ is a set of subsets of $P$ satisfying the two conditions:   
>> 1. the elements of $P$ belong to $K$,   
>> 2. if $\tau \in K$ and $\sigma \subset \tau$, then $\sigma \in K$.   
>>  
> > The elements of $K$ are the simplices.  
> {: .solution}
{: .discussion}
 
 Simplicial complexes can be seen at the same time as geometric/topological spaces (good for topological/geometrical inference) and as combinatorial objects (abstract simplicial complexes, good for computations)

> ## Exercise 1: Identify the simplices
>  In the following graph, we have 2 representations of simplicial complexes.
>  <a href="../fig/tda_08_exercise_1.png">
  <img src="../fig/tda_08_exercise_1.png" alt="Exercise 1" width="70%" height="auto"/>
</a>
> 
> How many simplices (0-simplices, 1-simplices, 2-simplices) do the simplicial complexes in the figure have?
> > ## Solution
> >   
> > |             | **Figure A** | **Figure B** |    
> > |:-----------:|:------------:|:------------:|   
> > | $0-simplex$ |       9      |      11      |    
> > | $1-simplex$ |      11      |      12      |    
> > | $2-simplex$ |       2      |       2      |    
> > 
> {: .solution}
{: .challenge}


### Čech and (Vietoris)-Rips complexes

The Vietoris-Rips complex and the Čech complex are two types of simplicial complexes used to construct discrete structures from sets of points in space.


The **Vietoris-Rips complex** is constructed from a set of points in a metric space. Given a set of points and a distance parameter called the "threshold," points within a distance less than or equal to the threshold are connected, forming the 1-simplices of the complex. Higher-dimensional simplices are then constructed by closing under combinations of 1-simplices that form a complete simplex, i.e., all fully connected subsets. The Vietoris-Rips complex captures the connectivity information between points and their topological structure at different scales.

Where the set of points represents your dataset, such as a group of genomes to study, and the distance could be a Hamming distance. 'Threshold' represents the selected distance for two genomes to be considered related."

> ## Vietoris complex
>
>> ## Mathematical definition
>> Given a point cloud $P=\\{p_1,...,p_n\\}\subset \mathbb{R}^d $, its **Rips complex** of radius $r>0$ is the simplicial complex $R(P,r)$ s.t. $vert(R(P,r))=P$  and
>>  $$ \sigma = [p_{i_0},p_{i_1},...,p_{i_k}] \in R(P,r) \quad iff \quad  \lVert p_{i_j} -p_{i_l}  \rVert  \leq 2r, \forall \leq j,l\leq k $$
>> 
> {: .solution}
{: .discussion}

On the other hand, the **Čech complex** is based on constructing simplicial cells rather than simply connecting points at specific distances. Given a set of points and a distance parameter, all sets of points whose balls of radius equal to the distance parameter have a non-empty intersection are considered. These sets of points become the simplices of the Čech complex. Similar to the Vietoris-Rips complex, higher-dimensional simplices can be constructed by closing under combinations of lower-dimensional simplices that form a complete simplex.

> ## Čech complex
>
>> ## Mathematical definition
>> Given a point cloud $P=\\{p_1,...,p_n\\}\subset \mathbb{R}^d$, its **Cech complex** of radius $r>0$ is the simplicial complex $C(P,r)$ s.t. $vert(C(P,r))=P$  and
>>  $$ \sigma = [p_{i_0},p_{i_1},...,p_{i_k}] \in C(P,r) \quad iff \quad \cap_{j=0}^k Bp_{i_j} \neq \emptyset $$
>> 
> {: .solution}
{: .discussion}
In the following figure, we have the representation of the Rips simplicial complex (left) and the Cech simplicial complex (right), for the same set of points $P={p_1,...,p_9}$.

 <a href="../fig/Tda_rips_cech.png">
  <img src="../fig/Tda_rips_cech.png" alt="Example Filtration" width="70%" height="auto" />
</a>

The following table shows the number of simplices for each simplicial complex.

 |             | $R(P,r)$ | $C(P,r)$ |    
 |:-----------:|:------------:|:------------:|   
 | $0-simplex$ |       9      |      9      |    
 | $1-simplex$ |      9      |      9      |    
 | $2-simplex$ |       2      |       1      |    


The difference lies in the $2-simplices$. Given 3 points, to construct a 2-simplex (triangle) in the Cech complex, the 3 balls of radius $r$ must intersect each other, while in the Rips complex, the balls must intersect pairwise."

Both the Vietoris-Rips complex and the Čech complex are tools used in topological analysis and computational geometry to study the structure and properties of sets of points in space. These complexes provide a discrete representation of the proximity and connectivity information of the points, enabling the analysis of their topology and geometric characteristics.

> ## Note
> Čech complexes can be quite hard to compute.
{: .callout}

## Simplicial homology

Simplicial homology is a technique used to quantify the topological structure of a simplicial complex. This technique is based on the identification of cycles and voids in the complex, which can be quantified by assigning integer values called "homology degrees". Simplicial homology is often used in topological data analysis to find patterns and structures in the data. Some definitions that we must keep in mind are the following.

**Holes:** Holes are empty regions or connected spaces in a simplicial complex. Simplicial homology allows for the detection and quantification of the presence of holes in the complex.

**Connected Components:** Connected components are sets of simplices in a simplicial complex that are connected to each other through shared simplices. Simplicial homology can identify and count the connected components in the complex.

**Betti Numbers:** Betti numbers are numerical invariants that measure the number of connected components and holes in a simplicial complex. Betti-0 ($\beta_0$) counts the number of connected components, while Betti-1 ($\beta_1$) counts the number of one-dimensional holes.


> ## Exercise 2:  Identify Betti number
>  In the following graph, we have 2 representations of simplicial complexes.
>  <a href="../fig/tda_08_exercise_1.png">
  <img src="../fig/tda_08_exercise_1.png" alt="Exercise 1" width="70%" height="auto" />
</a>
> 
> What are the $\beta_0$ and $\beta_1$ of these simplicial complexes in the image?.
>
> Hint: Remember that a painted triangle (color rose in this case) represents a 2-simplex, while a missing triangle forms a 1-hole.
> > ## Solution  
> >   
> >  |           | **Figure A** | **Figure B** |   
> >  |:---------:|:------------:|:------------:|   
> >  | $\beta_0$ |       1      |       3      |   
> >  | $\beta_1$ |       1      |       2      |    
> >   
> {: .solution}
{: .challenge}



### Filtration
  
A filtration of a simplicial complex is an ordered sequence of subcomplexes of the original complex, where each subcomplex contains its predecessor in the sequence. In other words, it is a way to decompose the complex into successive stages, where each stage adds or removes simplices compared to the previous stage.

> ## Filtration
>
> > ## Mathematical definition
>>A filtration of a simplicial complex $K$ is a collection $K_0 satisfying the two conditions: 
>>  1. $K_N=K$.
>> 2.  $K_i$ is a subcomplex of $K_{i+1}$, for $i=0,1,...,N-1$.
> {: .solution} 
{: .discussion}

An example of a filtration for $K=\langle a, b , c \rangle $ is the following:

  <a href="../fig/Tda-Filtacion1.png">
  <img src="../fig/Tda-Filtacion1.png" alt="Example Filtration" width="70%" height="auto"/>
</a>

In this case, we have 5 levels of filtration. The red color represents the simplices that are being added at each level of filtration, starting from $K_0$ with the vertices ${a,b,c}$ and ending at $K_4=\langle a, b, c \rangle $, the triangle (including the face) with vertices $a,b,c$.


Now, to apply persistent homology, we need to vary the parameters associated with the filtration. In the filtration graph, we have 5 distinct steps for the filtered simplicial complex. In each of these steps, we can have a different number of connected components and 1-holes. 


> ## Exercise 3:  Calculate the Betti numbers.
>  For the filtration shown in Figure 1, calculate the Betti numbers ($\beta_0 $  and $\beta_1$)for each level of the filtration
> 
> > ## Solution  
> >   
> >  |           | $\beta_0$ | $\beta_1$ |   
> >  |:---------:|:------------:|:------------:|   
> >  | $K_0$ |       3      |       0      |   
> >  | $K_1$ |       2      |       0      |    
> >  | $K_2$ |       1      |       0      |   
> >  | $K_3$ |       1      |       1      |    
> >  | $K_4$ |       1      |       0      |    
> >   
> {: .solution}
{: .challenge}

We can observe that at the beginning (filtration level 0), we have 3 connected components. At filtration level 1, we have 2 connected components, meaning that 2 connected components merged to form a single one, or we can also say that one connected component died. The same thing happens at filtration level 2, leaving us with only one connected component. At filtration level 3, we still have one connected component, but a 1-hole appears. Finally, at filtration level 4, we still have one connected component, and the 1-hole disappears.

To visualize these changes, we will use a representation with the following persistence diagrams and barcode.

### Persistence Diagram
The persistence diagram is a visual representation of the evolution of cycles and cavities in different dimensions as the simplicial complex is modified. It helps understand the persistence and relevance of topological structures in the complex.

Continuing with the example of the filtration of $K$, to construct the persistence diagram, we need to empty the information we obtained in the previous exercise.
  <a href="../fig/tda_08_diagrama.png">
  <img src="../fig/tda_08_diagrama.png" alt="Example Barcode Diagram" width="50%" height="auto"/>
</a>

In the persistence diagram, connected components are represented in red, and 1-holes are shown in blue. The $x$ and $y$ axes represent the filtration levels and are labeled as "birth" and "death" respectively. The coordinates of the points indicate when each connected component or 1-hole was born and died. Additionally, there is a connected component that has infinite death time, meaning it never dies during the filtration. This is referred to as "infinite persistence".


### Barcode Diagram
A barcode diagram is a graphical tool used to visualize the persistence diagram. It consists of bars that represent the persistence intervals of cycles and cavities, indicating their duration and relevance in the simplicial complex.
  <a href="../fig/tda_08_barcode.png">
  <img src="../fig/tda_08_barcode.png" alt="Example Persistence Diagram" width="50%" height="auto"/>
</a>

In the barcode diagram, the $X$ axis represents the filtration levels, and the beginning of a red bar indicates the birth of a connected component, while the end of that bar represents its death (i.e., when it merged with another component). For the blue bars, the idea is analogous, but they represent 1-holes instead.

## App to play

This app, created in GeoGebra by José María Ibarra Rodríguez (jose.ibarra@c3consensus.com), allows us to manipulate the parameter $r$  to construct a Rips simplicial complex for a set of 16 points in $\mathbb{R}^2$. Additionally, it enables us to view the barcode of the filtration and visualize the simplices as they are formed.

Steps to use:

1. Click on the square button in the bottom left corner to maximize the app.
2. In the top right corner, there is a slider for the values of $r$ that you can adjust according to your preferences.
3. At the bottom, you can select what you want to display, such as simplices, balls, and the barcode.

<iframe src="https://www.geogebra.org/classic/s7W7zbG4?embed" width="600" height="400" allowfullscreen style="border: 1px solid #e4e4e4;border-radius: 4px;" frameborder="0"></iframe>



<script>
  MathJax.Hub.Config({
    tex2jax: {
      inlineMath: [['$', '$'], ['\\(', '\\)']],
      displayMath: [['$$', '$$'], ['\\[', '\\]']],
      processEscapes: true,
      processEnvironments: true,
    }
  });
  MathJax.Hub.Queue(["Typeset", MathJax.Hub]);
</script>



> ## Exercise 4: Using the app example
>
> Using the application, could you please answer the following questions:
>  
> a. For what value of $r$ does the 1-simplex (edge) appear?  
> b. For what value of $r$ do we have only one connected component ($\beta_0 = 1$)?  
> c. How many 1-holes appear, and for what values of $r$ do they appear?  
> d. What is the persistence of the 1-holes?  
> e. What can that persistence tell us about the shape of the data?  
>  
> > ## Solution
> > a.  $r=0.5$  
> > b.  $r=0.83$  
> > c. Two. The first one $r=0.76$ and the second one $r=0.86$.  
> > d. The persistences are $0.12$ and $1.65$.  
> > e.  The second hole has a significant persistence of $1.6$, which suggests that the data is indeed arranged in a circular manner.  
> {: .solution}
{: .challenge}


