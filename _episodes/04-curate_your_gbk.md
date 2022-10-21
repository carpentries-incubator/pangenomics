---
title: "Annotating Genomic Data"
teaching: 30
exercises: 15
questions:
- "How to download NCBI genomic data from the command line?"``
- "How to annotate genome FASTA files?"
objectives:
- "Learn how to use the Prokka genome annotation utility."
keypoints:
- "Prokka is a command line utility that provides rapid prokaryotic genome annotation."
---

You need to make the file `correctgbk.sh`, for this you use `nano correctgbk.sh` and paste the next script. 
~~~
file=$1
locus=$(grep -m 1 "LOCUS" $file |cut -d\  -f 8 |cut -b1-9)  #selecionas de el primer Locus, los primeros caracteres donde empiezan con N
perl -p -i -e 's/\n// if /ORGANISM/' $file  #cambiar 
perl -p -i -e 's/\s*Unclassified/ '"${locus}"'/' $file
~~~
{: .language-bash}
~~~
ls *.gbk | while read file
do 
bash correctgbk.sh $file
done
~~~
{: .language-bash}
