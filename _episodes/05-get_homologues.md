---
title: "Get_Homologues"
teaching: 0
exercises: 0
questions:
- "Key question (FIXME)"
- "What is Get_Homologues?"
objectives:
- "First learning objective. (FIXME)"
keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---
# Get_Homologues
## Step 1
We need to create a new folder
~~~
mkdir get_homologues_clus
cd get_homologues_clus
~~~
{: .source}

~~~
betterlab@132.248.196.38's password:
~~~
{: .output}

## Step 2
Use clustering COGtriangle algorithm (COGS, PubMed=20439257)
~~~
get_homologues.pl -d /home/betterlab/GenomeMining/datos/gbk -e -G 
~~~
{: .source}

~~~
betterlab@132.248.196.38's password:
~~~
{: .output}
## Step 3
Use orthoMCL algorithm (OMCL, PubMed=12952885)
~~~
get_homologues.pl -d /home/betterlab/GenomeMining/datos/gbk -e -M 
~~~
{: .source}

~~~
betterlab@132.248.196.38's password:
~~~
{: .output}

{% include links.md %}
