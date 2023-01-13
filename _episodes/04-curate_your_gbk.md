---
title: "Curating output files"
teaching: 30
exercises: 15
questions:
- "How can I modify all my files?"
objectives:
- "Write a script to add the NCBI Id to prokka genbank files."
keypoints:
- "There are always some details that require manual curation."
- "Scripts can help you to automatize your work."
---
## Curating Prokka output files

Prokka's output files need to be corrected before moving forward with additional analyses. 
Create the file `correctgbk.sh`. We suggest the use of nano text editor to create your file, `nano correctgbk.sh` and paste the following script. 
~~~
file=$1
locus=$(grep -m 1 "DEFINITION" $file |cut -d " " -f6,7) # separate the DENFITION row by spaces and save columns 6 and 7 in the locus variable.

#NOTE If you don't have the strain details in your gbk files, comment the previous line and uncomment the next one.
#locus=$(grep -m 1 "LOCUS" $file |cut -d\  -f 8 |cut -b1-11)  #select the first 11 characters from the first "LOCUS"
perl -p -i -e 's/\n// if /ORGANISM/' $file  #Put the Organism row on a single line 
perl -p -i -e 's/\s*Unclassified/ '"${locus}"'/' $file    #Change Unclassfied to the value of the locus variable
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
