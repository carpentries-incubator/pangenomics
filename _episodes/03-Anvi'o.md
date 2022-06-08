---
title: "Interactive pangenome graphics"
teaching: 15 min
exercises: 40 min
questions:
- "What is Anvi'o?"
- "How can I obtain a pangenome analysis in Anvi'o?"
objectives:
- "Establish a dataset of genomes to obtain their pangenome"
- "Perform a basic workflow to obtain a pangenome in Anvi'o"
- "Understand and interpret the pangenome results"
keypoints:
- "Anvi’o is an open-source, community-driven analysis and visualization platform for microbial ‘omics. "
---


![Figure 1. Anvi'o network representation](../fig/anvio-network.png)

## Anvi'o

Anvi’o is an open-source, community-driven analysis and visualization platform for microbial ‘omics.
It brings together many aspects of today's cutting-edge strategies including **genomics, metagenomics, metatranscriptomics, phylogenomics, microbial population genetis, pangenomics and metapangenomis** in an *integrated* and *easy-to-use* fashion thorugh extensive interactive visualization capabilities. 



## The basic process to construct a pangenome starting with genbank files

To start using Anvi'o, activate the conda environment used to installation.
Instead of (base), the beginning of the line code will be (Pangenomics) indicated that the environment is active. 
~~~
conda activate Pangenomics
~~~
{: .source}

~~~
(Pangenomics) betterlab@betterlabub:~/Pangenomics/Anvio/MTBC$
~~~
{: .output}


Ten steps workflow to construct a Pangenome in Anvi'o
===============================================
Move into the directory results and create a directory for the Anvi'o analysis
~~~
cd results
mkdir anvi-o
cd anvi-o
~~~
{: .source}

Create a directory that will be used to storage all the files necessary to construct the genomes database
~~~
mkdir genome-db
~~~
{: .source}


1. Process the genome files (.gbk) with the `anvi-script-process-genbank` function. Remember that the .gbk files are not in this working directory
~~~
ls ~/Pangenomics/Shaday/gbk_ncbi/*.gbk | cut -d'/' -f7 | cut -d '.' -f1 | while read line; do anvi-script-process-genbank -i GENBANK --input-genbank ~/dc_workshop/results/annotated/$line.gbk -O genome-db/$line; done
~~~
{: .source}

~~~
cd genome-db
ls
~~~
{: .source}

~~~
agalactiae_18RS21-contigs.fa               agalactiae_A909-contigs.fa                 agalactiae_COH1-contigs.fa
agalactiae_18RS21-external-functions.txt   agalactiae_A909-external-functions.txt     agalactiae_COH1-external-functions.txt
agalactiae_18RS21-external-gene-calls.txt  agalactiae_A909-external-gene-calls.txt    agalactiae_COH1-external-gene-calls.txt
agalactiae_515-contigs.fa                  agalactiae_CJB111-contigs.fa               agalactiae_H36B-contigs.fa
agalactiae_515-external-functions.txt      agalactiae_CJB111-external-functions.txt   agalactiae_H36B-external-functions.txt
agalactiae_515-external-gene-calls.txt     agalactiae_CJB111-external-gene-calls.txt  agalactiae_H36B-external-gene-calls.txt
~~~
{: .output}


2. Reformat the fasta files

~~~
ls *fa |while read line; do anvi-script-reformat-fasta --seq-type NT $line -o $line\.fasta; done
~~~
{: .source}

~~~
agalactiae_18RS21-contigs.fa               agalactiae_A909-contigs.fa                 agalactiae_COH1-contigs.fa
agalactiae_18RS21-contigs.fa.fasta         agalactiae_A909-contigs.fa.fasta           agalactiae_COH1-contigs.fa.fasta
agalactiae_18RS21-external-functions.txt   agalactiae_A909-external-functions.txt     agalactiae_COH1-external-functions.txt
agalactiae_18RS21-external-gene-calls.txt  agalactiae_A909-external-gene-calls.txt    agalactiae_COH1-external-gene-calls.txt
agalactiae_515-contigs.fa                  agalactiae_CJB111-contigs.fa               agalactiae_H36B-contigs.fa
agalactiae_515-contigs.fa.fasta            agalactiae_CJB111-contigs.fa.fasta         agalactiae_H36B-contigs.fa.fasta
agalactiae_515-external-functions.txt      agalactiae_CJB111-external-functions.txt   agalactiae_H36B-external-functions.txt
agalactiae_515-external-gene-calls.txt     agalactiae_CJB111-external-gene-calls.txt  agalactiae_H36B-external-gene-calls.txt
~~~
{: .output}

3. Create a database per genome with the `anvi-gen-contigs-database` script

~~~
ls *fasta | while read line; do anvi-gen-contigs-database -T 4 -f $line -o $line-contigs.db; done
ls
~~~
{: .source}

~~~
agalactiae_18RS21-contigs.fa                   agalactiae_A909-contigs.fa                     agalactiae_COH1-contigs.fa
agalactiae_18RS21-contigs.fa.fasta             agalactiae_A909-contigs.fa.fasta               agalactiae_COH1-contigs.fa.fasta
agalactiae_18RS21-contigs.fa.fasta-contigs.db  agalactiae_A909-contigs.fa.fasta-contigs.db    agalactiae_COH1-contigs.fa.fasta-contigs.db
agalactiae_18RS21-external-functions.txt       agalactiae_A909-external-functions.txt         agalactiae_COH1-external-functions.txt
agalactiae_18RS21-external-gene-calls.txt      agalactiae_A909-external-gene-calls.txt        agalactiae_COH1-external-gene-calls.txt
agalactiae_515-contigs.fa                      agalactiae_CJB111-contigs.fa                   agalactiae_H36B-contigs.fa
agalactiae_515-contigs.fa.fasta                agalactiae_CJB111-contigs.fa.fasta             agalactiae_H36B-contigs.fa.fasta
agalactiae_515-contigs.fa.fasta-contigs.db     agalactiae_CJB111-contigs.fa.fasta-contigs.db  agalactiae_H36B-contigs.fa.fasta-contigs.db
agalactiae_515-external-functions.txt          agalactiae_CJB111-external-functions.txt       agalactiae_H36B-external-functions.txt
agalactiae_515-external-gene-calls.txt         agalactiae_CJB111-external-gene-calls.txt      agalactiae_H36B-external-gene-calls.txt
~~~
{: .output}


4. When using external genomes in anvi'o, a list of the genome ids and their corresponding genome database is required. This list will 
~~~
ls *.fa | cut -d '-' -f1 | while read line; do echo $line$'\t'$line-contigs.db >>external-genomes.txt; done
head external-genomes.txt
~~~
{: .source}

~~~
agalactiae_18RS21	agalactiae_18RS21-contigs.db
agalactiae_515	agalactiae_515-contigs.db
agalactiae_A909	agalactiae_A909-contigs.db
agalactiae_CJB111	agalactiae_CJB111-contigs.db
agalactiae_COH1	agalactiae_COH1-contigs.db
agalactiae_H36B	agalactiae_H36B-contigs.db
~~~
{: .output}

5. Modify the headers of the list external-genomes.txt
$ nano external-genomes.txt
  GNU nano 4.8                                                             external-genomes.txt                                                                       
agalactiae_18RS21       agalactiae_18RS21-contigs.db
agalactiae_515  agalactiae_515-contigs.db
agalactiae_A909 agalactiae_A909-contigs.db
agalactiae_CJB111       agalactiae_CJB111-contigs.db
agalactiae_COH1 agalactiae_COH1-contigs.db
agalactiae_H36B agalactiae_H36B-contigs.db





^G Get Help     ^O Write Out    ^W Where Is     ^K Cut Text     ^J Justify      ^C Cur Pos      M-U Undo        M-A Mark Text   M-] To Bracket  M-Q Previous
^X Exit         ^R Read File    ^\ Replace      ^U Paste Text   ^T To Spell     ^_ Go To Line   M-E Redo        M-6 Copy Text   ^Q Where Was    M-W Next

head external-genomes.txt
name	contigs_db_path
agalactiae_18RS21	agalactiae_18RS21-contigs.db
agalactiae_515	agalactiae_515-contigs.db
agalactiae_A909	agalactiae_A909-contigs.db
agalactiae_CJB111	agalactiae_CJB111-contigs.db
agalactiae_COH1	agalactiae_COH1-contigs.db
agalactiae_H36B	agalactiae_H36B-contigs.db

6. Rename the .db files
$ rename s'/.fa.fasta-contigs.db/.db/' *db
$ ls *.db
agalactiae_18RS21-contigs.db  agalactiae_A909-contigs.db    agalactiae_COH1-contigs.db
agalactiae_515-contigs.db     agalactiae_CJB111-contigs.db  agalactiae_H36B-contigs.db

7. Create the Pangenome database
$ anvi-gen-genomes-storage -e external-genomes.txt -o AGALACTIAE_GENOMES.db
$ ls
AGALACTIAE_GENOMES.db                      agalactiae_A909-contigs.db                 agalactiae_COH1-contigs.fa
agalactiae_18RS21-contigs.db               agalactiae_A909-contigs.fa                 agalactiae_COH1-contigs.fa.fasta
agalactiae_18RS21-contigs.fa               agalactiae_A909-contigs.fa.fasta           agalactiae_COH1-external-functions.txt
agalactiae_18RS21-contigs.fa.fasta         agalactiae_A909-external-functions.txt     agalactiae_COH1-external-gene-calls.txt
agalactiae_18RS21-external-functions.txt   agalactiae_A909-external-gene-calls.txt    agalactiae_H36B-contigs.db
agalactiae_18RS21-external-gene-calls.txt  agalactiae_CJB111-contigs.db               agalactiae_H36B-contigs.fa
agalactiae_515-contigs.db                  agalactiae_CJB111-contigs.fa               agalactiae_H36B-contigs.fa.fasta
agalactiae_515-contigs.fa                  agalactiae_CJB111-contigs.fa.fasta         agalactiae_H36B-external-functions.txt
agalactiae_515-contigs.fa.fasta            agalactiae_CJB111-external-functions.txt   agalactiae_H36B-external-gene-calls.txt
agalactiae_515-external-functions.txt      agalactiae_CJB111-external-gene-calls.txt  external-genomes.txt
agalactiae_515-external-gene-calls.txt     agalactiae_COH1-contigs.db






Step 1
===============================================
**Process the genome files (.gbk)**

This script takes a GenBank file, and outputs a FASTA file, as well as two additional TAB-delimited output files for external gene calls and gene functions that can be used with the programs anvi-gen-contigs-database and anvi-import-functions.

~~~
ls *gbk | cut -d'.' -f1 | while read line; do echo anvi-script-process-genbank -i GENBANK --input-genbank $line\.gbk -O $line; anvi-script-process-genbank -i GENBANK --input-genbank $line\.gbk -O $line; done
~~~
{: .source}

~~~
anvi-script-process-genbank -i GENBANK --input-genbank Mtb_N0004_L3.gbk -O Mtb_N0004_L3
Num GenBank entries processed ................: 94
Num gene records found .......................: 4,202
Num genes reported ...........................: 4,202
Num genes with AA sequences ..................: 4,202
Num genes with functions .....................: 3,347
Num partial genes ............................: 0
Num genes excluded ...........................: 0

FASTA file path ..............................: Mtb_N0004_L3-contigs.fa
External gene calls file .....................: Mtb_N0004_L3-external-gene-calls.txt
TAB-delimited functions ......................: Mtb_N0004_L3-external-functions.txt

* Mmmmm ☘

.
.
.

anvi-script-process-genbank -i GENBANK --input-genbank Mtb_N3913_L7.gbk -O Mtb_N3913_L7
Num GenBank entries processed ................: 115
Num gene records found .......................: 4,312
Num genes reported ...........................: 4,312
Num genes with AA sequences ..................: 4,312
Num genes with functions .....................: 3,429
Num partial genes ............................: 0
Num genes excluded ...........................: 0

FASTA file path ..............................: Mtb_N3913_L7-contigs.fa
External gene calls file .....................: Mtb_N3913_L7-external-gene-calls.txt
TAB-delimited functions ......................: Mtb_N3913_L7-external-functions.txt

* Mmmmm ☘
~~~
{: .output}

Then, you can explore the files created 

~~~
ls
~~~
{: .source}

~~~
Mtb_N0004_L3-contigs.fa               Mtb_N0091_L6.gbk                      Mtb_N1202_L6-external-gene-calls.txt
Mtb_N0004_L3-external-functions.txt   Mtb_N0136_L4-contigs.fa               Mtb_N1202_L6.gbk
Mtb_N0004_L3-external-gene-calls.txt  Mtb_N0136_L4-external-functions.txt   Mtb_N1216_L4-contigs.fa
Mtb_N0004_L3.gbk                      Mtb_N0136_L4-external-gene-calls.txt  Mtb_N1216_L4-external-functions.txt
Mtb_N0031_L2-contigs.fa               Mtb_N0136_L4.gbk                      Mtb_N1216_L4-external-gene-calls.txt
Mtb_N0031_L2-external-functions.txt   Mtb_N0145_L2-contigs.fa               Mtb_N1216_L4.gbk
Mtb_N0031_L2-external-gene-calls.txt  Mtb_N0145_L2-external-functions.txt   Mtb_N1268_L5-contigs.fa
Mtb_N0031_L2.gbk                      Mtb_N0145_L2-external-gene-calls.txt  Mtb_N1268_L5-external-functions.txt
Mtb_N0052_L2-contigs.fa               Mtb_N0145_L2.gbk                      Mtb_N1268_L5-external-gene-calls.txt
Mtb_N0052_L2-external-functions.txt   Mtb_N0155_L2-contigs.fa               Mtb_N1268_L5.gbk
Mtb_N0052_L2-external-gene-calls.txt  Mtb_N0155_L2-external-functions.txt   Mtb_N1272_L5-contigs.fa
Mtb_N0052_L2.gbk                      Mtb_N0155_L2-external-gene-calls.txt  Mtb_N1272_L5-external-functions.txt
Mtb_N0054_L3-contigs.fa               Mtb_N0155_L2.gbk                      Mtb_N1272_L5-external-gene-calls.txt
Mtb_N0054_L3-external-functions.txt   Mtb_N0157_L1-contigs.fa               Mtb_N1272_L5.gbk
Mtb_N0054_L3-external-gene-calls.txt  Mtb_N0157_L1-external-functions.txt   Mtb_N1274_L3-contigs.fa
Mtb_N0054_L3.gbk                      Mtb_N0157_L1-external-gene-calls.txt  Mtb_N1274_L3-external-functions.txt
Mtb_N0069_L1-contigs.fa               Mtb_N0157_L1.gbk                      Mtb_N1274_L3-external-gene-calls.txt
Mtb_N0069_L1-external-functions.txt   Mtb_N1176_L5-contigs.fa               Mtb_N1274_L3.gbk
Mtb_N0069_L1-external-gene-calls.txt  Mtb_N1176_L5-external-functions.txt   Mtb_N1283_L4-contigs.fa
Mtb_N0069_L1.gbk                      Mtb_N1176_L5-external-gene-calls.txt  Mtb_N1283_L4-external-functions.txt
Mtb_N0072_L1-contigs.fa               Mtb_N1176_L5.gbk                      Mtb_N1283_L4-external-gene-calls.txt
Mtb_N0072_L1-external-functions.txt   Mtb_N1201_L6-contigs.fa               Mtb_N1283_L4.gbk
Mtb_N0072_L1-external-gene-calls.txt  Mtb_N1201_L6-external-functions.txt   Mtb_N3913_L7-contigs.fa
Mtb_N0072_L1.gbk                      Mtb_N1201_L6-external-gene-calls.txt  Mtb_N3913_L7-external-functions.txt
Mtb_N0091_L6-contigs.fa               Mtb_N1201_L6.gbk                      Mtb_N3913_L7-external-gene-calls.txt
Mtb_N0091_L6-external-functions.txt   Mtb_N1202_L6-contigs.fa               Mtb_N3913_L7.gbk
Mtb_N0091_L6-external-gene-calls.txt  Mtb_N1202_L6-external-functions.txt
~~~
{: .output}


Step 2
===============================================
**Reformat the fasta files** 


~~~
ls *fa | while read line; do anvi-script-reformat-fasta $line -o $line\.fasta; done
~~~
{: .source}

~~~
Input ........................................: Mtb_N0004_L3-contigs.fa
Output .......................................: Mtb_N0004_L3-contigs.fa.fasta
Minimum length ...............................: 0
Max % gaps allowed ...........................: 100.00%
Total num contigs ............................: 94
Total num nucleotides ........................: 4,249,707
Contigs removed ..............................: 0 (0.00% of all)
Nucleotides removed ..........................: 0 (0.00% of all)
Nucleotides modified .........................: 0 (0.00% of all)
Deflines simplified ..........................: False

.
.
.

Input ........................................: Mtb_N3913_L7-contigs.fa
Output .......................................: Mtb_N3913_L7-contigs.fa.fasta
Minimum length ...............................: 0
Max % gaps allowed ...........................: 100.00%
Total num contigs ............................: 115
Total num nucleotides ........................: 4,326,679
Contigs removed ..............................: 0 (0.00% of all)
Nucleotides removed ..........................: 0 (0.00% of all)
Nucleotides modified .........................: 0 (0.00% of all)
Deflines simplified ..........................: False
~~~
{: .output}

Now, you can see the formated files (*.fasta) created  

~~~
ls
~~~
{: .source}

~~~
Mtb_N0004_L3-contigs.fa               Mtb_N0091_L6.gbk                      Mtb_N1202_L6-external-gene-calls.txt
Mtb_N0004_L3-contigs.fa.fasta         Mtb_N0136_L4-contigs.fa               Mtb_N1202_L6.gbk
Mtb_N0004_L3-external-functions.txt   Mtb_N0136_L4-contigs.fa.fasta         Mtb_N1216_L4-contigs.fa
Mtb_N0004_L3-external-gene-calls.txt  Mtb_N0136_L4-external-functions.txt   Mtb_N1216_L4-contigs.fa.fasta
Mtb_N0004_L3.gbk                      Mtb_N0136_L4-external-gene-calls.txt  Mtb_N1216_L4-external-functions.txt
Mtb_N0031_L2-contigs.fa               Mtb_N0136_L4.gbk                      Mtb_N1216_L4-external-gene-calls.txt
Mtb_N0031_L2-contigs.fa.fasta         Mtb_N0145_L2-contigs.fa               Mtb_N1216_L4.gbk
Mtb_N0031_L2-external-functions.txt   Mtb_N0145_L2-contigs.fa.fasta         Mtb_N1268_L5-contigs.fa
Mtb_N0031_L2-external-gene-calls.txt  Mtb_N0145_L2-external-functions.txt   Mtb_N1268_L5-contigs.fa.fasta
Mtb_N0031_L2.gbk                      Mtb_N0145_L2-external-gene-calls.txt  Mtb_N1268_L5-external-functions.txt
Mtb_N0052_L2-contigs.fa               Mtb_N0145_L2.gbk                      Mtb_N1268_L5-external-gene-calls.txt
Mtb_N0052_L2-contigs.fa.fasta         Mtb_N0155_L2-contigs.fa               Mtb_N1268_L5.gbk
Mtb_N0052_L2-external-functions.txt   Mtb_N0155_L2-contigs.fa.fasta         Mtb_N1272_L5-contigs.fa
Mtb_N0052_L2-external-gene-calls.txt  Mtb_N0155_L2-external-functions.txt   Mtb_N1272_L5-contigs.fa.fasta
Mtb_N0052_L2.gbk                      Mtb_N0155_L2-external-gene-calls.txt  Mtb_N1272_L5-external-functions.txt
Mtb_N0054_L3-contigs.fa               Mtb_N0155_L2.gbk                      Mtb_N1272_L5-external-gene-calls.txt
Mtb_N0054_L3-contigs.fa.fasta         Mtb_N0157_L1-contigs.fa               Mtb_N1272_L5.gbk
Mtb_N0054_L3-external-functions.txt   Mtb_N0157_L1-contigs.fa.fasta         Mtb_N1274_L3-contigs.fa
Mtb_N0054_L3-external-gene-calls.txt  Mtb_N0157_L1-external-functions.txt   Mtb_N1274_L3-contigs.fa.fasta
Mtb_N0054_L3.gbk                      Mtb_N0157_L1-external-gene-calls.txt  Mtb_N1274_L3-external-functions.txt
Mtb_N0069_L1-contigs.fa               Mtb_N0157_L1.gbk                      Mtb_N1274_L3-external-gene-calls.txt
Mtb_N0069_L1-contigs.fa.fasta         Mtb_N1176_L5-contigs.fa               Mtb_N1274_L3.gbk
Mtb_N0069_L1-external-functions.txt   Mtb_N1176_L5-contigs.fa.fasta         Mtb_N1283_L4-contigs.fa
Mtb_N0069_L1-external-gene-calls.txt  Mtb_N1176_L5-external-functions.txt   Mtb_N1283_L4-contigs.fa.fasta
Mtb_N0069_L1.gbk                      Mtb_N1176_L5-external-gene-calls.txt  Mtb_N1283_L4-external-functions.txt
Mtb_N0072_L1-contigs.fa               Mtb_N1176_L5.gbk                      Mtb_N1283_L4-external-gene-calls.txt
Mtb_N0072_L1-contigs.fa.fasta         Mtb_N1201_L6-contigs.fa               Mtb_N1283_L4.gbk
Mtb_N0072_L1-external-functions.txt   Mtb_N1201_L6-contigs.fa.fasta         Mtb_N3913_L7-contigs.fa
Mtb_N0072_L1-external-gene-calls.txt  Mtb_N1201_L6-external-functions.txt   Mtb_N3913_L7-contigs.fa.fasta
Mtb_N0072_L1.gbk                      Mtb_N1201_L6-external-gene-calls.txt  Mtb_N3913_L7-external-functions.txt
Mtb_N0091_L6-contigs.fa               Mtb_N1201_L6.gbk                      Mtb_N3913_L7-external-gene-calls.txt
Mtb_N0091_L6-contigs.fa.fasta         Mtb_N1202_L6-contigs.fa               Mtb_N3913_L7.gbk
Mtb_N0091_L6-external-functions.txt   Mtb_N1202_L6-contigs.fa.fasta
Mtb_N0091_L6-external-gene-calls.txt  Mtb_N1202_L6-external-functions.txt
~~~
{: .output}

~~~
ls *.fasta
~~~
{: .source}

~~~
Mtb_N0004_L3-contigs.fa.fasta  Mtb_N0072_L1-contigs.fa.fasta  Mtb_N0157_L1-contigs.fa.fasta  Mtb_N1268_L5-contigs.fa.fasta
Mtb_N0031_L2-contigs.fa.fasta  Mtb_N0091_L6-contigs.fa.fasta  Mtb_N1176_L5-contigs.fa.fasta  Mtb_N1272_L5-contigs.fa.fasta
Mtb_N0052_L2-contigs.fa.fasta  Mtb_N0136_L4-contigs.fa.fasta  Mtb_N1201_L6-contigs.fa.fasta  Mtb_N1274_L3-contigs.fa.fasta
Mtb_N0054_L3-contigs.fa.fasta  Mtb_N0145_L2-contigs.fa.fasta  Mtb_N1202_L6-contigs.fa.fasta  Mtb_N1283_L4-contigs.fa.fasta
Mtb_N0069_L1-contigs.fa.fasta  Mtb_N0155_L2-contigs.fa.fasta  Mtb_N1216_L4-contigs.fa.fasta  Mtb_N3913_L7-contigs.fa.fasta
~~~
{: .output}

Step 3
===============================================
**Create a database per genome**

~~~
ls *fasta | while read line; do anvi-gen-contigs-database -T 4 -f $line -o $line-contigs.db; done
~~~
{: .source}

~~~
Input FASTA file .............................: /home/betterlab/Pangenomics/MTBC/Mtb_N0004_L3-contigs.fa.fasta

Anvi'o made things up for you
===============================================
You are generating a new anvi'o contigs database, but you are not specifying a
project name for it. FINE. Anvi'o, in desperation, will use the input file name
to set the project name for this contigs database (i.e.,
'Mtb_N0004_L3_contigs_fa'). If you are not happy with that, feel free to kill
and restart this process. If you are not happy with this name, but you don't
like killing things either, maybe next time you should either name your FASTA
files better, or use the `--project-name` parameter to set your desired name.

Name .........................................: Mtb_N0004_L3_contigs_fa
Description ..................................: No description is given
Num threads for gene calling .................: 4                                                                                                 

Finding ORFs in contigs
===============================================
Genes ........................................: /tmp/tmp80we2s58/contigs.genes
Amino acid sequences .........................: /tmp/tmp80we2s58/contigs.amino_acid_sequences
Log file .....................................: /tmp/tmp80we2s58/00_log.txt

CITATION
===============================================
Anvi'o will use 'prodigal' by Hyatt et al (doi:10.1186/1471-2105-11-119) to
identify open reading frames in your data. When you publish your findings,
please do not forget to properly credit their work.

Result .......................................: Prodigal (v2.6.3) has identified 3890 genes.                                                      

                                                                                                                                                  
CONTIGS DB CREATE REPORT
===============================================
Split Length .................................: 20,000
K-mer size ...................................: 4
Skip gene calling? ...........................: False
External gene calls provided? ................: False
Ignoring internal stop codons? ...............: False
Splitting pays attention to gene calls? ......: True
Contigs with at least one gene call ..........: 94 of 94 (100.0%)                                                                                 
Contigs database .............................: A new database, Mtb_N0004_L3-contigs.fa.fasta-contigs.db, has been created.
Number of contigs ............................: 94
Number of splits .............................: 209
Total number of nucleotides ..................: 4,249,707
Gene calling step skipped ....................: False
Splits broke genes (non-mindful mode) ........: False
Desired split length (what the user wanted) ..: 20,000
Average split length (what anvi'o gave back) .: 22,168

✓ anvi-gen-contigs-database took 0:00:14.507798

.
.
.

Input FASTA file .............................: /home/betterlab/Pangenomics/MTBC/Mtb_N3913_L7-contigs.fa.fasta

Anvi'o made things up for you
===============================================
You are generating a new anvi'o contigs database, but you are not specifying a
project name for it. FINE. Anvi'o, in desperation, will use the input file name
to set the project name for this contigs database (i.e.,
'Mtb_N3913_L7_contigs_fa'). If you are not happy with that, feel free to kill
and restart this process. If you are not happy with this name, but you don't
like killing things either, maybe next time you should either name your FASTA
files better, or use the `--project-name` parameter to set your desired name.

Name .........................................: Mtb_N3913_L7_contigs_fa
Description ..................................: No description is given
Num threads for gene calling .................: 4                                                                                                 

Finding ORFs in contigs
===============================================
Genes ........................................: /tmp/tmpxgytnn6r/contigs.genes
Amino acid sequences .........................: /tmp/tmpxgytnn6r/contigs.amino_acid_sequences
Log file .....................................: /tmp/tmpxgytnn6r/00_log.txt

CITATION
===============================================
Anvi'o will use 'prodigal' by Hyatt et al (doi:10.1186/1471-2105-11-119) to
identify open reading frames in your data. When you publish your findings,
please do not forget to properly credit their work.

Result .......................................: Prodigal (v2.6.3) has identified 3988 genes.                                                      

                                                                                                                                                  
CONTIGS DB CREATE REPORT
===============================================
Split Length .................................: 20,000
K-mer size ...................................: 4
Skip gene calling? ...........................: False
External gene calls provided? ................: False
Ignoring internal stop codons? ...............: False
Splitting pays attention to gene calls? ......: True
Contigs with at least one gene call ..........: 115 of 115 (100.0%)                                                                               
Contigs database .............................: A new database, Mtb_N3913_L7-contigs.fa.fasta-contigs.db, has been created.
Number of contigs ............................: 115
Number of splits .............................: 237
Total number of nucleotides ..................: 4,326,679
Gene calling step skipped ....................: False
Splits broke genes (non-mindful mode) ........: False
Desired split length (what the user wanted) ..: 20,000
Average split length (what anvi'o gave back) .: 21,928

✓ anvi-gen-contigs-database took 0:00:14.542257
~~~
{: .output}

Now, you must obtained the database files created (*.db)

~~~
ls
~~~
{: .source}

~~~
Mtb_N0004_L3-contigs.fa                   Mtb_N0091_L6.gbk                          Mtb_N1202_L6-external-gene-calls.txt
Mtb_N0004_L3-contigs.fa.fasta             Mtb_N0136_L4-contigs.fa                   Mtb_N1202_L6.gbk
Mtb_N0004_L3-contigs.fa.fasta-contigs.db  Mtb_N0136_L4-contigs.fa.fasta             Mtb_N1216_L4-contigs.fa
Mtb_N0004_L3-external-functions.txt       Mtb_N0136_L4-contigs.fa.fasta-contigs.db  Mtb_N1216_L4-contigs.fa.fasta
Mtb_N0004_L3-external-gene-calls.txt      Mtb_N0136_L4-external-functions.txt       Mtb_N1216_L4-contigs.fa.fasta-contigs.db
Mtb_N0004_L3.gbk                          Mtb_N0136_L4-external-gene-calls.txt      Mtb_N1216_L4-external-functions.txt
Mtb_N0031_L2-contigs.fa                   Mtb_N0136_L4.gbk                          Mtb_N1216_L4-external-gene-calls.txt
Mtb_N0031_L2-contigs.fa.fasta             Mtb_N0145_L2-contigs.fa                   Mtb_N1216_L4.gbk
Mtb_N0031_L2-contigs.fa.fasta-contigs.db  Mtb_N0145_L2-contigs.fa.fasta             Mtb_N1268_L5-contigs.fa
Mtb_N0031_L2-external-functions.txt       Mtb_N0145_L2-contigs.fa.fasta-contigs.db  Mtb_N1268_L5-contigs.fa.fasta
Mtb_N0031_L2-external-gene-calls.txt      Mtb_N0145_L2-external-functions.txt       Mtb_N1268_L5-contigs.fa.fasta-contigs.db
Mtb_N0031_L2.gbk                          Mtb_N0145_L2-external-gene-calls.txt      Mtb_N1268_L5-external-functions.txt
Mtb_N0052_L2-contigs.fa                   Mtb_N0145_L2.gbk                          Mtb_N1268_L5-external-gene-calls.txt
Mtb_N0052_L2-contigs.fa.fasta             Mtb_N0155_L2-contigs.fa                   Mtb_N1268_L5.gbk
Mtb_N0052_L2-contigs.fa.fasta-contigs.db  Mtb_N0155_L2-contigs.fa.fasta             Mtb_N1272_L5-contigs.fa
Mtb_N0052_L2-external-functions.txt       Mtb_N0155_L2-contigs.fa.fasta-contigs.db  Mtb_N1272_L5-contigs.fa.fasta
Mtb_N0052_L2-external-gene-calls.txt      Mtb_N0155_L2-external-functions.txt       Mtb_N1272_L5-contigs.fa.fasta-contigs.db
Mtb_N0052_L2.gbk                          Mtb_N0155_L2-external-gene-calls.txt      Mtb_N1272_L5-external-functions.txt
Mtb_N0054_L3-contigs.fa                   Mtb_N0155_L2.gbk                          Mtb_N1272_L5-external-gene-calls.txt
Mtb_N0054_L3-contigs.fa.fasta             Mtb_N0157_L1-contigs.fa                   Mtb_N1272_L5.gbk
Mtb_N0054_L3-contigs.fa.fasta-contigs.db  Mtb_N0157_L1-contigs.fa.fasta             Mtb_N1274_L3-contigs.fa
Mtb_N0054_L3-external-functions.txt       Mtb_N0157_L1-contigs.fa.fasta-contigs.db  Mtb_N1274_L3-contigs.fa.fasta
Mtb_N0054_L3-external-gene-calls.txt      Mtb_N0157_L1-external-functions.txt       Mtb_N1274_L3-contigs.fa.fasta-contigs.db
Mtb_N0054_L3.gbk                          Mtb_N0157_L1-external-gene-calls.txt      Mtb_N1274_L3-external-functions.txt
Mtb_N0069_L1-contigs.fa                   Mtb_N0157_L1.gbk                          Mtb_N1274_L3-external-gene-calls.txt
Mtb_N0069_L1-contigs.fa.fasta             Mtb_N1176_L5-contigs.fa                   Mtb_N1274_L3.gbk
Mtb_N0069_L1-contigs.fa.fasta-contigs.db  Mtb_N1176_L5-contigs.fa.fasta             Mtb_N1283_L4-contigs.fa
Mtb_N0069_L1-external-functions.txt       Mtb_N1176_L5-contigs.fa.fasta-contigs.db  Mtb_N1283_L4-contigs.fa.fasta
Mtb_N0069_L1-external-gene-calls.txt      Mtb_N1176_L5-external-functions.txt       Mtb_N1283_L4-contigs.fa.fasta-contigs.db
Mtb_N0069_L1.gbk                          Mtb_N1176_L5-external-gene-calls.txt      Mtb_N1283_L4-external-functions.txt
Mtb_N0072_L1-contigs.fa                   Mtb_N1176_L5.gbk                          Mtb_N1283_L4-external-gene-calls.txt
Mtb_N0072_L1-contigs.fa.fasta             Mtb_N1201_L6-contigs.fa                   Mtb_N1283_L4.gbk
Mtb_N0072_L1-contigs.fa.fasta-contigs.db  Mtb_N1201_L6-contigs.fa.fasta             Mtb_N3913_L7-contigs.fa
Mtb_N0072_L1-external-functions.txt       Mtb_N1201_L6-contigs.fa.fasta-contigs.db  Mtb_N3913_L7-contigs.fa.fasta
Mtb_N0072_L1-external-gene-calls.txt      Mtb_N1201_L6-external-functions.txt       Mtb_N3913_L7-contigs.fa.fasta-contigs.db
Mtb_N0072_L1.gbk                          Mtb_N1201_L6-external-gene-calls.txt      Mtb_N3913_L7-external-functions.txt
Mtb_N0091_L6-contigs.fa                   Mtb_N1201_L6.gbk                          Mtb_N3913_L7-external-gene-calls.txt
Mtb_N0091_L6-contigs.fa.fasta             Mtb_N1202_L6-contigs.fa                   Mtb_N3913_L7.gbk
Mtb_N0091_L6-contigs.fa.fasta-contigs.db  Mtb_N1202_L6-contigs.fa.fasta             
Mtb_N0091_L6-external-functions.txt       Mtb_N1202_L6-contigs.fa.fasta-contigs.db
Mtb_N0091_L6-external-gene-calls.txt      Mtb_N1202_L6-external-functions.txt
~~~
{: .output}

~~~
ls *.db
~~~
{: .source}

~~~
Mtb_N0004_L3-contigs.fa.fasta-contigs.db  Mtb_N0136_L4-contigs.fa.fasta-contigs.db  Mtb_N1216_L4-contigs.fa.fasta-contigs.db
Mtb_N0031_L2-contigs.fa.fasta-contigs.db  Mtb_N0145_L2-contigs.fa.fasta-contigs.db  Mtb_N1268_L5-contigs.fa.fasta-contigs.db
Mtb_N0052_L2-contigs.fa.fasta-contigs.db  Mtb_N0155_L2-contigs.fa.fasta-contigs.db  Mtb_N1272_L5-contigs.fa.fasta-contigs.db
Mtb_N0054_L3-contigs.fa.fasta-contigs.db  Mtb_N0157_L1-contigs.fa.fasta-contigs.db  Mtb_N1274_L3-contigs.fa.fasta-contigs.db
Mtb_N0069_L1-contigs.fa.fasta-contigs.db  Mtb_N1176_L5-contigs.fa.fasta-contigs.db  Mtb_N1283_L4-contigs.fa.fasta-contigs.db
Mtb_N0072_L1-contigs.fa.fasta-contigs.db  Mtb_N1201_L6-contigs.fa.fasta-contigs.db  Mtb_N3913_L7-contigs.fa.fasta-contigs.db
Mtb_N0091_L6-contigs.fa.fasta-contigs.db  Mtb_N1202_L6-contigs.fa.fasta-contigs.db
~~~
{: .output}


Step 4
===============================================
**Create a list of ids and their corresponding genome database**

~~~
ls *fa |cut -d'-' -f1 |while read line; do echo $line$'\t'$line-contigs.db >>external-genomes.txt; done
~~~
{: .source}

You can explore this created list

~~~
head external-genomes.txt
~~~
{: .source}

~~~
Mtb_N0004_L3	Mtb_N0004_L3-contigs.db
Mtb_N0031_L2	Mtb_N0031_L2-contigs.db
Mtb_N0052_L2	Mtb_N0052_L2-contigs.db
Mtb_N0054_L3	Mtb_N0054_L3-contigs.db
Mtb_N0069_L1	Mtb_N0069_L1-contigs.db
Mtb_N0072_L1	Mtb_N0072_L1-contigs.db
Mtb_N0091_L6	Mtb_N0091_L6-contigs.db
Mtb_N0136_L4	Mtb_N0136_L4-contigs.db
Mtb_N0145_L2	Mtb_N0145_L2-contigs.db
Mtb_N0155_L2	Mtb_N0155_L2-contigs.db
~~~
{: .output}

~~~
tail external-genomes.txt
~~~
{: .source}

~~~
Mtb_N0157_L1	Mtb_N0157_L1-contigs.db
Mtb_N1176_L5	Mtb_N1176_L5-contigs.db
Mtb_N1201_L6	Mtb_N1201_L6-contigs.db
Mtb_N1202_L6	Mtb_N1202_L6-contigs.db
Mtb_N1216_L4	Mtb_N1216_L4-contigs.db
Mtb_N1268_L5	Mtb_N1268_L5-contigs.db
Mtb_N1272_L5	Mtb_N1272_L5-contigs.db
Mtb_N1274_L3	Mtb_N1274_L3-contigs.db
Mtb_N1283_L4	Mtb_N1283_L4-contigs.db
Mtb_N3913_L7	Mtb_N3913_L7-contigs.db
~~~
{: .output}

Step 5
===============================================
**Modify the headers of the list external-genomes.txt**
~~~
nano external-genomes.txt
~~~
{: .source}

~~~
  GNU nano 4.8                                                   external-genomes.txt                                                             
Mtb_N0004_L3    Mtb_N0004_L3-contigs.db
Mtb_N0031_L2    Mtb_N0031_L2-contigs.db
Mtb_N0052_L2    Mtb_N0052_L2-contigs.db
Mtb_N0054_L3    Mtb_N0054_L3-contigs.db
Mtb_N0069_L1    Mtb_N0069_L1-contigs.db
Mtb_N0072_L1    Mtb_N0072_L1-contigs.db
Mtb_N0091_L6    Mtb_N0091_L6-contigs.db
Mtb_N0136_L4    Mtb_N0136_L4-contigs.db
Mtb_N0145_L2    Mtb_N0145_L2-contigs.db
Mtb_N0155_L2    Mtb_N0155_L2-contigs.db
Mtb_N0157_L1    Mtb_N0157_L1-contigs.db
Mtb_N1176_L5    Mtb_N1176_L5-contigs.db
Mtb_N1201_L6    Mtb_N1201_L6-contigs.db
Mtb_N1202_L6    Mtb_N1202_L6-contigs.db
Mtb_N1216_L4    Mtb_N1216_L4-contigs.db
Mtb_N1268_L5    Mtb_N1268_L5-contigs.db
Mtb_N1272_L5    Mtb_N1272_L5-contigs.db
Mtb_N1274_L3    Mtb_N1274_L3-contigs.db
Mtb_N1283_L4    Mtb_N1283_L4-contigs.db
Mtb_N3913_L7    Mtb_N3913_L7-contigs.db


^G Get Help     ^O Write Out    ^W Where Is     ^K Cut Text     ^J Justify      ^C Cur Pos      M-U Undo        M-A Mark Text   M-] To Bracket
^X Exit         ^R Read File    ^\ Replace      ^U Paste Text   ^T To Spell     ^_ Go To Line   M-E Redo        M-6 Copy Text   ^Q Where Was
~~~
{: .output}

We need to include the column headers by typing at the first row the titles <<name>> followed by 'tab' and then <<contigs_db_path>> 
To exit from nano, press **Ctrl+X**, respond **Yes** (Y) to save the changes and press **Enter** to use the same file name (external-genomes.txt).

When we explore this file, we must see the new headers of the columns

~~~
head external-genomes.txt
~~~
{: .source}

~~~
name	contigs_db_path
Mtb_N0004_L3	Mtb_N0004_L3-contigs.db
Mtb_N0031_L2	Mtb_N0031_L2-contigs.db
Mtb_N0052_L2	Mtb_N0052_L2-contigs.db
Mtb_N0054_L3	Mtb_N0054_L3-contigs.db
Mtb_N0069_L1	Mtb_N0069_L1-contigs.db
Mtb_N0072_L1	Mtb_N0072_L1-contigs.db
Mtb_N0091_L6	Mtb_N0091_L6-contigs.db
Mtb_N0136_L4	Mtb_N0136_L4-contigs.db
Mtb_N0145_L2	Mtb_N0145_L2-contigs.db
~~~
{: .output}  

Step 6
===============================================
**Rename the .db files**
  
~~~
rename s'/.fa.fasta-contigs.db/.db/' *db
ls *.db
~~~
{: .source}

~~~
Mtb_N0004_L3-contigs.db  Mtb_N0069_L1-contigs.db  Mtb_N0145_L2-contigs.db  Mtb_N1201_L6-contigs.db  Mtb_N1272_L5-contigs.db
Mtb_N0031_L2-contigs.db  Mtb_N0072_L1-contigs.db  Mtb_N0155_L2-contigs.db  Mtb_N1202_L6-contigs.db  Mtb_N1274_L3-contigs.db
Mtb_N0052_L2-contigs.db  Mtb_N0091_L6-contigs.db  Mtb_N0157_L1-contigs.db  Mtb_N1216_L4-contigs.db  Mtb_N1283_L4-contigs.db
Mtb_N0054_L3-contigs.db  Mtb_N0136_L4-contigs.db  Mtb_N1176_L5-contigs.db  Mtb_N1268_L5-contigs.db  Mtb_N3913_L7-contigs.db
~~~
{: .output}

Step 7
===============================================
**Create the pangenome database**
  
~~~
anvi-gen-genomes-storage -e external-genomes.txt -o MTBC_GENOMES.db
~~~
{: .source}

~~~
WARNING
===============================================
None of your genomes seem to have any functional annotation. No biggie. Things
will continue to work. But then your genomes have no functional annotation. SAD.

Internal genomes .............................: 0 have been initialized.                                                                          
External genomes .............................: 20 found.                                                                                         
                                                                                                                                                  
WARNING
===============================================
The contigs databases you are using for this analysis are missing HMMs for
single-copy core genes. Maybe you haven't run `anvi-run-hmms` on your contigs
database, or they didn't contain any hits. It is perfectly legal to have anvi'o
contigs databases without HMMs or SCGs for things to work, but we wanted to give
you heads up so you can have your 'aha' moment if you see funny things in the
interface.

                                                                                                                                                  
* Mtb_N0004_L3 is stored with 3,890 genes (108 of which were partial)
* Mtb_N0031_L2 is stored with 3,964 genes (136 of which were partial)                                                                             
* Mtb_N0052_L2 is stored with 3,935 genes (116 of which were partial)                                                                             
* Mtb_N0054_L3 is stored with 3,995 genes (113 of which were partial)                                                                             
* Mtb_N0069_L1 is stored with 4,009 genes (137 of which were partial)                                                                             
* Mtb_N0072_L1 is stored with 3,991 genes (114 of which were partial)                                                                             
* Mtb_N0091_L6 is stored with 3,975 genes (135 of which were partial)                                                                             
* Mtb_N0136_L4 is stored with 3,976 genes (207 of which were partial)                                                                             
* Mtb_N0145_L2 is stored with 3,975 genes (126 of which were partial)                                                                             
* Mtb_N0155_L2 is stored with 3,971 genes (134 of which were partial)                                                                             
* Mtb_N0157_L1 is stored with 3,980 genes (143 of which were partial)                                                                             
* Mtb_N1176_L5 is stored with 4,034 genes (158 of which were partial)                                                                             
* Mtb_N1201_L6 is stored with 3,968 genes (130 of which were partial)                                                                             
* Mtb_N1202_L6 is stored with 3,975 genes (139 of which were partial)                                                                             
* Mtb_N1216_L4 is stored with 3,962 genes (147 of which were partial)                                                                             
* Mtb_N1268_L5 is stored with 3,992 genes (155 of which were partial)                                                                             
* Mtb_N1272_L5 is stored with 3,978 genes (155 of which were partial)                                                                             
* Mtb_N1274_L3 is stored with 3,998 genes (115 of which were partial)                                                                             
* Mtb_N1283_L4 is stored with 3,994 genes (215 of which were partial)                                                                             
* Mtb_N3913_L7 is stored with 3,988 genes (138 of which were partial)                                                                             

The new genomes storage ......................: MTBC_GENOMES.db (v7, signature: hashb220d200)
Number of genomes ............................: 20 (internal: 0, external: 20)
Number of gene calls .........................: 79,550
Number of partial gene calls .................: 2,821
~~~
{: .output}

Now you must obtained a new *.db file containing all the information about your genomes of interest

~~~
ls *.db
~~~
{: .source}

~~~
MTBC_GENOMES.db          Mtb_N0069_L1-contigs.db  Mtb_N0155_L2-contigs.db  Mtb_N1216_L4-contigs.db  Mtb_N3913_L7-contigs.db
Mtb_N0004_L3-contigs.db  Mtb_N0072_L1-contigs.db  Mtb_N0157_L1-contigs.db  Mtb_N1268_L5-contigs.db
Mtb_N0031_L2-contigs.db  Mtb_N0091_L6-contigs.db  Mtb_N1176_L5-contigs.db  Mtb_N1272_L5-contigs.db
Mtb_N0052_L2-contigs.db  Mtb_N0136_L4-contigs.db  Mtb_N1201_L6-contigs.db  Mtb_N1274_L3-contigs.db
Mtb_N0054_L3-contigs.db  Mtb_N0145_L2-contigs.db  Mtb_N1202_L6-contigs.db  Mtb_N1283_L4-contigs.db
~~~
{: .output}

  
Step 8
===============================================
**Construct the pangenome with the database created above**
~~~
anvi-pan-genome -g MTBC_GENOMES.db \
                --project-name "Pangenome_MTBC" \
                --output-dir Pangenome-MTBC \
                --num-threads 6 \
                --minbit 0.5 \
                --mcl-inflation 10 \
                --use-ncbi-blast
~~~
{: .source}

~~~
WARNING
===============================================
If you publish results from this workflow, please do not forget to cite DIAMOND
(doi:10.1038/nmeth.3176), unless you use it with --use-ncbi-blast flag, and MCL
(http://micans.org/mcl/ and doi:10.1007/978-1-61779-361-5_15)

Functions found ..............................:                                                                                                   
Genomes storage ..............................: Initialized (storage hash: hashb220d200)                                                          
Num genomes in storage .......................: 20
Num genomes will be used .....................: 20
Pan database .................................: A new database, /home/betterlab/Pangenomics/MTBC/Pangenome-MTBC/Pangenome_MTBC-PAN.db, has been   
                                                created.
Exclude partial gene calls ...................: False

AA sequences FASTA ...........................: /home/betterlab/Pangenomics/MTBC/Pangenome-MTBC/combined-aas.fa                                   

Num AA sequences reported ....................: 79,550
Num excluded gene calls ......................: 0
Unique AA sequences FASTA ....................: /home/betterlab/Pangenomics/MTBC/Pangenome-MTBC/combined-aas.fa.unique                            

WARNING
===============================================
You elected to use NCBI's `blastp` for amino acid sequence search. Running
blastp will be significantly slower than DIAMOND, but in some cases, slightly
more sensitive. We are unsure about whether the slight increase in sensitivity
may justify significant increase in run time, but you are the boss.


NCBI BLAST MAKEDB
===============================================
BLAST search db ..............................: /home/betterlab/Pangenomics/MTBC/Pangenome-MTBC/combined-aas.fa.unique                            

NCBI BLAST SEARCH
===============================================
BLAST results ................................: /home/betterlab/Pangenomics/MTBC/Pangenome-MTBC/blast-search-results.txt                          

MCL INPUT
===============================================
Min percent identity .........................: 0.0                                                                                               
Minbit .......................................: 0.5
Filtered search results ......................: 1,863,310 edges stored                                                                            
MCL input ....................................: /home/betterlab/Pangenomics/MTBC/Pangenome-MTBC/mcl-input.txt

MCL
===============================================
MCL inflation ................................: 10.0
MCL output ...................................: /home/betterlab/Pangenomics/MTBC/Pangenome-MTBC/mcl-clusters.txt                                  
Number of MCL clusters .......................: 4,192
                                                                                                                                                  
CITATION
===============================================
The workflow you are using will likely use 'muscle' by Edgar,
doi:10.1093/nar/gkh340 (http://www.drive5.com/muscle) to align your sequences.
If you publish your findings, please do not forget to properly credit this tool.

* Your pangenome is ready with a total of 4,192 gene clusters across 20 genomes 🎉 
~~~
{: .output}

You can find the results into the Pangenome-MTBC directory created 

~~~
ls Pangenome-MTBC/
~~~
{: .source}

~~~
Pangenome_MTBC-PAN.db            combined-aas.fa.unique        combined-aas.fa.unique.pin  combined-aas.fa.unique.pto
blast-search-results.txt         combined-aas.fa.unique.names  combined-aas.fa.unique.pot  log.txt
blast-search-results.txt.unique  combined-aas.fa.unique.pdb    combined-aas.fa.unique.psq  mcl-clusters.txt
combined-aas.fa                  combined-aas.fa.unique.phr    combined-aas.fa.unique.ptf  mcl-input.txt
~~~
{: .output}

Step 9
===============================================
**Create the imagen of the results**

~~~
anvi-display-pan -g MTBC_GENOMES.db \
    -p Pangenome-MTBC/Pangenome_MTBC-PAN.db
~~~
{: .source}

~~~
Interactive mode .............................: pan                                                                                               
Functions found .............................................:                                                                                    
Genomes storage .............................................: Initialized (storage hash: hashb220d200)                                           
Num genomes in storage ......................................: 20
Num genomes will be used ....................................: 20
Pan DB ......................................................: Initialized: Pangenome-MTBC/Pangenome_MTBC-PAN.db (v. 15)
Gene cluster homogeneity estimates ..........................: Functional: [YES]; Geometric: [YES]; Combined: [YES]
                                                                                                                                                  
* Gene clusters are initialized for all 4192 gene clusters in the database.

                                                                                                                                                  
WARNING
==============================================================
Genomes storage does not have any info about gene functions. Certain parts of
the pangenomic workflow will not be accessible.

                                                                                                                                                  
WARNING
==============================================================
Genomes storage does not have any info about gene functions. Certain parts of
the pangenomic workflow will not be accessible.


WARNING
==============================================================
Someone asked anvi'o to initialize a gene cluster functions summary dict, but it
seems there are no gene cluster functions even after initializing functions for
the pangenome. So we move on without any summary dict for functions and/or drama
about it to let the downstream analyses decide how to punish the unlucky.

                                                                                                                                                  
* The server is up and running 🎉

WARNING
===============================================
If you are using OSX and if the server terminates prematurely before you can see
anything in your browser, try running the same command by putting 'sudo ' at the
beginning of it (you will be prompted to enter your password if sudo requires
super user credentials on your system). If your browser does not show up, try
manually entering the URL shown below into the address bar of your favorite
browser. *cough* CHROME *cough*.


Server address ...............................: http://0.0.0.0:8080

* When you are ready, press CTRL+C once to terminate the server and go back to the
command line.
~~~
{: .output}

Step 10
===============================================
**Visualize the results by cliking the link of the server address**

~~~
http://132.248.196.38:8080
~~~
{: .output}

> ## Exercise 2. Splitting the pangenome.
> 1. Read about [anvi-split](https://anvio.org/help/main/programs/anvi-split/) 
> 2. With this program split your pangenome in independent pangenomes that:
> > * Contains only singletons.
> > * Contains only core gene clusters.
>
> Tip: [anvi-display-pan](https://anvio.org/help/main/programs/anvi-display-pan/) can be usefull
{: .challenge}

{% include links.md %}
