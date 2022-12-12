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
## Curating Prokka output files

Prokka's output files need to be corrected before moving forward with additional analyses. 
Create the file `correctgbk.sh`. We suggest the use of nano text editor to create your file, `nano correctgbk.sh` and paste the next script. 
~~~
file=$1
locus=$(grep -m 1 "DEFINITION" $file |cut -d " " -f6,7) #if you have details the strain in your gbk files, use this line. Else use the next line.
#locus=$(grep -m 1 "LOCUS" $file |cut -d\  -f 8 |cut -b1-11)  #select the first 11 characters from the first "LOCUS"
perl -p -i -e 's/\n// if /ORGANISM/' $file  #cambiar 
perl -p -i -e 's/\s*Unclassified/ '"${locus}"'/' $file
~~~
{: .language-bash}

This script change "Unclassified" in the lines "Organisms" from the gbk files by the "strain" or "locus" depending on what you selected.
~~~
ls *.gbk | while read file
do 
bash correctgbk.sh $file
done
~~~
{: .language-bash}
