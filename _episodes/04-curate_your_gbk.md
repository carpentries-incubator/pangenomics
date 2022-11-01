---
title: "Curating output files"
teaching: 30
exercises: 15
questions:
- "How can I modify all my files?"
objectives:
- "Write a script to add the NCBI Id to prokka genbank files."
keypoints:
- "There are always some details that requires manual curation."
- "Scripts can help you to automatize your work."
---

You need to make the file `correctgbk.sh`, for this you use `nano correctgbk.sh` and paste the next script. 
~~~
file=$1
locus=$(grep -m 1 "LOCUS" $file |cut -d\  -f 8 |cut -b1-11)  #selecionas de el primer Locus, los primeros 11 caracteres donde empiezan con N
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
