---
title: "Topological data analysis"
teaching: 30
exercises: 15
questions:
- "What is topological data analysis?"
objectives:
- "Understand a simplex"
keypoints:
- "TDA describes data forms"
---

## **Topological data analysis**

Topological data analysis is a technique that uses concepts from topology to analyze complex data and find patterns and structures that are not apparent at first glance. This technique is based on the construction of a structure called a simplicial complex, which is composed of a collection of simple geometric objects called simplices. The topology of this complex is used to analyze and visualize the relationships between the data.

### **Simplicial complex**


Given a set $P=\{p_0,...,p_k\}\subset \mathbb{R}^d$ of $k+1$ affinely independent points, the **k-dimensional simplex** $\sigma$ (or **k-simplex** for short) spanned by $P$ is the set of convez combinatios
$$ \sum_{i=0}^k\lambda_ip_i, \quad with \quad  \sum_{i=0}^k\lambda_i = 1 \quad \lambda_i \geq 0. $$ 
The points $p_0, ..., p_k$ are called the vertices of $\sigma$.
>Simplex: Simple geometric object of any dimension (point, line segment, triangle, tetrahedron, etc.) that is used to construct simplicial complexes.

  <a href="../fig/tda_Vertices.png">
  <img src="../fig/tda_Vertices.png" alt="Example" />
</a>

A **simplicial complex** $K$ in $\mathbb{R}^d$ is a collection of simplices s.t:
1. any face a simplex of $K$ is a simplex of K,
2. the intesection of any twho simplices of $K$ is ether empty or a common face of both.

  <a href="../fig/tda_paste.png">
  <img src="../fig/tda_paste.png" alt="Example" />
</a>

>Simplicial complex: Mathematical structure composed of a collection of simple geometric objects called simplices, which is constructed from a set of data.

### **Abstract simplex and simplical complex**
Let $P= \{p_1,...,p_n\}$ be a (finite) set. An **abstract simlicial complex** $K$ with vertex set $P$ is a set of subsets of $P$ satisfying the two conditions:
1. the elements of $P$ belong to $K$,
2. if $\tau \in K$ and $\sigma \subset \tau$, then $\sigma \in K$.

The elements of $K$ are the simplices.


> Note: Simplicial complexes can be seen at the same time as geometric/topological
spaces (good for topological/geometrical inference) and as combinatorial objects (abstract simplicial complexes, good for computations)

  ### **Filtration**
  A filtration of a simplicial complex $K$ is a collection $K_0 \subset K_1 \subset ... \subset K_N$ of complexes such that:
  1. $K_N=K$.
  2. $K_i$ is a subcomplex of $K_{i+1}$, for $i=0,1,...,N-1$.

  <a href="../fig/Tda-Filtacion1.png">
  <img src="../fig/Tda-Filtacion1.png" alt="Example Filtration" />
</a>

### **Cech and (Vietoris)-Rips complexes**

Given a point cloud $P=\{p_1,...,p_n\}\subset \mathbb{R}^d$, its **Cech complex** of radius $r>0$ is the simplicial complex $C(P,r)$ s.t. $vert(C(P,r))=P$ and
$$ \sigma = [p_{i_0},p_{i_1},...,p_{i_k}] \in C(P,r) \quad iff \quad \cap_{j=0}^k B(p_{i_j} \neq \emptyset $$

> Cech complexes can be quite hard to compute.

Given a point cloud $P=\{p_1,...,p_n\}\subset \mathbb{R}^d$, its **Rips complex** of radius $r>0$ is the simplicial complex $R(P,r)$ s.t. $vert(R(P,r))=P$ and
$$ \sigma = [p_{i_0},p_{i_1},...,p_{i_k}] \in R(P,r) \quad iff \quad  \lVert p_{i_j} -p_{i_l}  \rVert  \leq 2r, \forall \leq j,l\leq k $$


 <a href="../fig/Tda_rips_cech.png">
  <img src="../fig/Tda_rips_cech.png" alt="Example Filtration" />
</a>

## **Simplicial homology**

Simplicial homology is a technique used to quantify the topological structure of a simplicial complex. This technique is based on the identification of cycles and voids in the complex, which can be quantified by assigning integer values called "homology degrees". Simplicial homology is often used in topological data analysis to find patterns and structures in the data.

>Homology degree: Integer assigned to a cycle or cavity in a simplicial complex using the technique of simplicial homology. The homology degree is used to quantify the "amount" of topological structure present in the complex.

>Cycle: Topological object that closes in on itself and has no holes. A cycle in a simplicial complex can be represented as a linear combination of simplices.

>Cavity: Topological object that represents a hole in the simplicial complex. A cavity in a simplicial complex is represented as a linear combination of simplices.


<iframe src="https://www.geogebra.org/classic/s7W7zbG4?embed" width="800" height="600" allowfullscreen style="border: 1px solid #e4e4e4;border-radius: 4px;" frameborder="0"></iframe>
