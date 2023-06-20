---
title: "Downloading Genomic Data"
teaching: 30
exercises: 15
questions:
- "How to download public genomes by using the command line?"
objectives:
- "Explore `ncbi-genome-download` as a tool for fetching genomic data from the NCBI."
keypoints:
- "The `ncbi-genome-download` package is a set of scripts designed to download genomes from the NCBI."
---

## Getting Genomic Data from the NCBI

We already have the genomes of strains 18RS21 and H36B in our pan_workshop/data directory. However, the remaining six GBS strains will be downloaded in this episode. We will obtain these genomic sequences from the National Center for Biotechnology Information (NCBI) database, the primary source of publicly available genomes. To automate the downloading process, we are going to utilize the specialized `ncbi-genome-download` [package](https://github.com/kblin/ncbi-genome-download), which includes convenient shell commands that allows users to download genomes directly from the NCBI. This package offers great flexibility by enabling users to specify their desired queries. It simplifies the process of retrieving the data and ensures that is conveniently saved into the working directory. 

The `ncbi-genome-download` package can be installed with conda. In our case, we have already installed it into the environment under the same name of the package *ncbi-genome-download*. Thus, in order to use the package, we just have to activate the *ncbi-genome-download* conda environment. 

> ## Know more
> If you want to know more about what is *conda* and its *environments* visit this [link](https://docs.conda.io/en/latest/).
{: .callout}


Let's activate the *ncbi-genome-download* conda environment to begin.  
~~~
$ conda activate ncbi-genome-download
~~~
{: .language-bash}

~~~
(ncbi-genome-download) $
~~~
{: .output}

For practicality, the prompt will be written only as `$` instead of `(ncbi-genome-download) $`.

Now, you are able to run the package `ncbi-genome-download`. 
Exploring the range of options available in the package is highly reccomended in order to choose well and get what you really need. To access the complete list of parameters to incorporate in your downloads, simply type the following commmand: 

~~~
$ ncbi-genome-download --help
~~~
{: .language-bash}

~~~
usage:  
 ncbi-genome-download [-h] [-s {refseq,genbank}] [-F FILE_FORMATS]  
                      	[-l ASSEMBLY_LEVELS] [-g GENERA] [--genus GENERA]  
                      	[--fuzzy-genus] [-S STRAINS] [-T SPECIES_TAXIDS]  
                      	[-t TAXIDS] [-A ASSEMBLY_ACCESSIONS]  
                      	[-R REFSEQ_CATEGORIES]  
                      	[--refseq-category REFSEQ_CATEGORIES] [-o OUTPUT]  
                      	[--flat-output] [-H] [-P] [-u URI] [-p N] [-r N]  
                      	[-m METADATA_TABLE] [-n] [-N] [-v] [-d] [-V]  
                      	[-M TYPE_MATERIALS]
                      	groups
	-F FILE_FORMATS, --formats FILE_FORMATS  
                    	Which formats to download (default: genbank).A comma-
                    	separated list of formats is also possible. For
                    	example: "fasta,assembly-report". Choose from:
                    	['genbank', 'fasta', 'rm', 'features', 'gff',
                    	'protein-fasta', 'genpept', 'wgs', 'cds-fasta', 'rna-
                    	fna', 'rna-fasta', 'assembly-report', 'assembly-
                    	stats', 'all']
	-g GENERA, --genera GENERA  
                    	Only download sequences of the provided genera. A
                    	comma-separated list of genera is also possible. For
                    	example: "Streptomyces coelicolor,Escherichia coli".
                    	(default: [])  
	-S STRAINS, --strains STRAINS   
                    	Only download sequences of the given strain(s). A
                    	comma-separated list of strain names is possible, as
                    	well as a path to a filename containing one name per
                    	line.
	-A ASSEMBLY_ACCESSIONS, --assembly-accessions ASSEMBLY_ACCESSIONS  
                    	Only download sequences matching the provided NCBI
                    	assembly accession(s). A comma-separated list of
                    	accessions is possible, as well as a path to a
                    	filename containing one accession per line.
	-o OUTPUT, --output-folder OUTPUT   
                    	Create output hierarchy in specified folder (default:
                    	/home/betterlab)
	-n, --dry-run     	Only check which files to download, don't download
                    	genome files.  
~~~
{: .output}

If you type `ncbi-genome-download` and you get the error `command-not-found`,
it could be because you are in the `base` and not inside the *ncbi-genome-download*
conda environment. Come back to the previous instruction and let's move on. 

Now, we have to move into our `data/` directory

~~~
$ cd ~/pan_workshop/data
~~~
{: .language-bash}

If you list the contents of this directory (using the `ls` command),
you'll see several directories, each of which contains the raw data
of different strains of *Streptococcus agalactiae* used
in [Tettelin *et al*., (2005)](https://www.pnas.org/doi/10.1073/pnas.0506758102)
in `.gbk` and `.fasta` formats.

~~~
$ ls
~~~
{: .language-bash}

~~~
agalactiae_18RS21  agalactiae_H36B annotated_mini
~~~
{: .output}

Prior to downloading anything from the NCBI, we recommend to verify if the genome or genomes you are interested are available in the database. The package `ncbi-genome-download` includes the `--dry-run` or `-n` flag, which means that the algorithm only will check the genomes you specify to download, but without downloading the files. Other useful flags are `--formats`, which serve to specify the desired format; `--genera` to specify the genus (or species), and `-S` for the strains. Importantly, the `ncbi-genome-download` command must always end with the name of the group of organisms where the search will be performed, which in our case will be `bacteria`. 

So, first, let's check if one of the genomes we are interested to download, "*Streptococcus agalactiae* 515" is available in NCBI. We will use the flags mentioned above.

~~~
$ ncbi-genome-download -n --formats fasta --genera "Streptococcus agalactiae" -S 515 bacteria 
~~~
{: .language-bash}
~~~
Considering the following 1 assemblies for download:
GCF_012593885.1 Streptococcus agalactiae 515    515
~~~
{: .output}

Great! The genome is available! 

Now, we can proceed to download it. To better organize our data, we can save this file into an specific directory for this strain. We can indicate this instruction with the `--output-folder` or `-o` flag followed by the name we choose. In this case, will be `-o agalactie_515`. Notice that now we no longer need the flag the `-n`.  

~~~
$ ncbi-genome-download --formats fasta --genera "Streptococcus agalactiae" -S 515 -o agalactiae_515 bacteria 
~~~
{: .language-bash}


Once the prompt `$` appears again, use the command `tree` to show the contents of the recently downloaded directory `agalactiae_515`.

~~~
$ tree agalactiae_515
~~~
{: .language-bash}
~~~
agalactiae_515
└── refseq
    └── bacteria
        └── GCF_012593885.1
            ├── GCF_012593885.1_ASM1259388v1_genomic.fna.gz
            └── MD5SUMS

3 directories, 2 files
~~~
{: .output}

The genome file `GCF_012593885.1_ASM1259388v1_genomic.fna.gz` 
is a compressed file located inside the directory
`agalactiae_515/refseq/bacteria/GCF_012593885.1/`. Let's
decompress the file with `gunzip` and visualize with `tree`
to corroborate the file status.

~~~
$ gunzip agalactiae_515/refseq/bacteria/GCF_012593885.1/GCF_012593885.1_ASM1259388v1_genomic.fna.gz
$ tree agalactiae_515/
~~~
{: .language-bash}
~~~
agalactiae_515/
└── refseq
    └── bacteria
        └── GCF_012593885.1
            ├── GCF_012593885.1_ASM1259388v1_genomic.fna
            └── MD5SUMS

3 directories, 2 files
~~~
{: .output}

`GCF_012593885.1_ASM1259388v1_genomic.fna` is now with `fna` extension
which means is in a nucleotide `fasta` format. Let's move the file to the
`agalactiae_515/` directory and remove the extra content that we will not 
use again in this lesson.

> ## MD5SUMS file
> Apart from the fasta file that we wanted, a file called `MD5SUMS` was also downloaded. This file has a unique code that identifies the contents of the 
> files of interest, so you can use it to check the integrity of your downloaded copy. We will not cover that step in the lesson but you can check 
> this [article](https://www.geeksforgeeks.org/md5sum-linux-command/) to see how you can use it.
{: .callout}

~~~
$ mv agalactiae_515/refseq/bacteria/GCF_012593885.1/GCF_012593885.1_ASM1259388v1_genomic.fna agalactiae_515/.
$ rm -r agalactiae_515/refseq
$ ls agalactiae_515/
~~~
{: .language-bash}   
~~~ 
GCF_012593885.1_ASM1259388v1_genomic.fna  
~~~
{: .output}

## Download multiple genomes

Right now, you have the _S. agalactiae 505_ genomic `fasta` file in a directory.
However, five more strains were included in Tettelin's first pangenome. It is the moment
to practice what you've learned about loops in the shell lesson. Instead
of downloading the genomes one by one, we will write a `while` loop.  

With the `nano` editor, create a file to add the other four strains. The missing strains are A909, 
COH1, CJB111, NEM316 and 2603V/R. Write one strain per line in the file and name it
"TettelinList.txt".

~~~
$ nano TettelinList.txt  
~~~
{: .language-bash}

Visualize "Tettelin.txt" contents with the `cat` command. 

~~~
$ cat TettelinList.txt  
~~~
{: .language-bash}
~~~
A909  
COH1  
CJB111 
NEM316
2603V/R
~~~
{: .output}

First, let's read the lines of the file file ina loop, and print
them in the terminal with the `echo strain $line` command.  
`strain` is just a word that we will print, and `$line` will 
store the value of each of the lines of the `Tettelin.txt` file.

~~~
$ cat TettelinList.txt | while read line 
> do 
> echo strain $line
> done
~~~
{: .language-bash}
~~~
strain A909  
strain COH1  
strain CJB111 
strain NEM316
strain 2603V/R
~~~
{: .output}

We can now check if these strains are available in NCBI (remember to use
the `-n` flag so genome files aren't downloaded).

~~~
$ cat TettelinList.txt | while read line; 
> do
> echo strain $line
> ncbi-genome-download --formats fasta --genera "Streptococcus agalactiae" -S $line -n bacteria
> done
~~~
{: .language-bash}
~~~ 
strain A909  
Considering the following 1 assemblies for download:  
GCF_000012705.1 Streptococcus agalactiae A909   A909  
strain COH1  
Considering the following 1 assemblies for download:  
GCF_000689235.1 Streptococcus agalactiae COH1   COH1  
strain CJB111  
Considering the following 2 assemblies for download:  
GCF_000167755.1 Streptococcus agalactiae CJB111 CJB111  
GCF_015221735.2 Streptococcus agalactiae CJB111 CJB111  
strain NEM316
Considering the following 1 assemblies for download:
GCF_000196055.1 Streptococcus agalactiae NEM316 NEM316
strain 2603V/R
Considering the following 1 assemblies for download:
GCF_000007265.1 Streptococcus agalactiae 2603V/R        2603V/R
~~~
{: .output}

The tool has successfully found the five strain. Notice that
the strain CJB111 contains two versions.

We can now proceed to download these strains to their corresponding
output directories by adding the `-o` flag followed by the directory
name and removing the `-n` flag).

~~~
$ cat TettelinList.txt | while read line 
> do
> echo downloading strain $line
> ncbi-genome-download --formats fasta --genera "Streptococcus agalactiae" -S $line -o agalactiae_$line bacteria
> done
~~~
{: .language-bash}
~~~
downloading strain A909
downloading strain COH1
downloading strain CJB111
downloading strain NEM316
downloading strain 2603V/R
~~~
{: .output}

Just as before, we should decompress the downloaded genome files using `gunzip`.
To do so, we can use the `*` wildcard, which means "anything", instead of unzipping
one by one.

~~~
$ gunzip agalactiae_*/refseq/bacteria/*/*gz
 ~~~
{: .language-bash}

Let's visualize the structure of the results
~~~
$ tree agalactiae_*
~~~
{: .language-bash}

~~~
agalactiae_2603V
└── R
    └── refseq
        └── bacteria
            └── GCF_000007265.1
                ├── GCF_000007265.1_ASM726v1_genomic.fna.gz
                └── MD5SUMS
agalactiae_515
└── GCF_012593885.1_ASM1259388v1_genomic.fna
agalactiae_A909
└── refseq
    └── bacteria
        └── GCF_000012705.1
            ├── GCF_000012705.1_ASM1270v1_genomic.fna
            └── MD5SUMS
agalactiae_CJB111
└── refseq
    └── bacteria
        ├── GCF_000167755.1
        │   ├── GCF_000167755.1_ASM16775v1_genomic.fna
        │   └── MD5SUMS
        └── GCF_015221735.2
            ├── GCF_015221735.2_ASM1522173v2_genomic.fna
            └── MD5SUMS
agalactiae_COH1
└── refseq
    └── bacteria
        └── GCF_000689235.1
            ├── GCF_000689235.1_GBCO_p1_genomic.fna
            └── MD5SUMS
agalactiae_H36B
├── Streptococcus_agalactiae_H36B.fna
└── Streptococcus_agalactiae_H36B.gbk
agalactiae_NEM316
└── refseq
    └── bacteria
        └── GCF_000196055.1
            ├── GCF_000196055.1_ASM19605v1_genomic.fna
            └── MD5SUMS

3 directories, 2 files
~~~
{: .output}

We noticed that all fasta files but  `GCF_000007265.1_ASM726v1_genomic.fna.gz` have been decompressed.
That decompression failure was because the 2603V/R strain has a different directory structure. This structure
is consequence of the name of the strain, because the characters "/R" are part of the name,
a directory named `R` has been added to the output, changing the directory structure.
Differences like this are expected to occur in big datasets, and must be manually
curated after the general cases has been treated with scripts. In this case the `tree`
command has helped us to identify that the error. Let's decompress the file 
`GCF_000007265.1_ASM726v1_genomic.fna.gz` and move it to the `agalactiae_2603V/` directory, although it doesn't have the real strain name.

~~~
$  gunzip agalactiae_2603V/R/refseq/bacteria/*/*gz
$  mv  agalactiae_2603V/R/refseq/bacteria/GCF_000007265.1/GCF_000007265.1_ASM726v1_genomic.fna agalactiae_2603V/
$  rm -r agalactiae_2603V/R/
$  ls agalactiae_2603V
~~~
{: .language-bash}

~~~
GCF_000007265.1_ASM726v1_genomic.fna
~~~
{: .output}

Finally, we need to move the other genome files to their corresponding locations and
get rid of unnecessary directories. To do so, we'll use a `while` cycle as follows.  
Beware of the typos! Take it slowly and make sure you are sending the files to the correct location. 

~~~
 $ cat TettelinList.txt | while read line
 > do 
 > echo moving fasta file of strain $line
 > mv agalactiae_$line/refseq/bacteria/*/*.fna agalactiae_$line/. 
 > done
 ~~~
{: .language-bash}
~~~
moving fasta file of strain A909
moving fasta file of strain COH1
moving fasta file of strain CJB111
moving fasta file of strain NEM316
moving fasta file of strain 2603V/R
mv: cannot stat 'agalactiae_2603V/R/refseq/bacteria/*/*.fna': No such file or directory
~~~
{: .output}

Thats ok, it is just telling us that the `agalactiae_2603V/R/` does not have an `fna` file, which is what we wanted.

Use the `tree` command to make sure that everything is in its right place.

Now let's remove the `refseq/` directories completely:

~~~
 $ cat TettelinList.txt | while read line
 > do 
 > echo removing refseq directory of strain $line
 > rm -r agalactiae_$line/refseq
 > done
 ~~~
{: .language-bash}

~~~
removing refseq directory of strain A909
removing refseq directory of strain COH1
removing refseq directory of strain CJB111
removing refseq directory of strain NEM316
removing refseq directory of strain 2603V/R
rm: cannot remove 'agalactiae_2603V/R/refseq': No such file or directory
~~~
{: .output}

At this point, you should have eight directories starting with `agalactiae_` containing
the following:

~~~
$ tree agalactiae_*
~~~
{: .language-bash}

~~~
agalactiae_18RS21
├── Streptococcus_agalactiae_18RS21.fna
└── Streptococcus_agalactiae_18RS21.gbk
agalactiae_2603V
└── GCF_000007265.1_ASM726v1_genomic.fna
agalactiae_515
└── GCF_012593885.1_ASM1259388v1_genomic.fna
agalactiae_A909
└── GCF_000012705.1_ASM1270v1_genomic.fna
agalactiae_CJB111
├── GCF_000167755.1_ASM16775v1_genomic.fna
└── GCF_015221735.2_ASM1522173v2_genomic.fna
agalactiae_COH1
└── GCF_000689235.1_GBCO_p1_genomic.fna
agalactiae_H36B
├── Streptococcus_agalactiae_H36B.fna
└── Streptococcus_agalactiae_H36B.gbk
agalactiae_NEM316
└── GCF_000196055.1_ASM19605v1_genomic.fna

0 directories, 1 file
~~~
{: .output}
We can see that the strain `CJB111` has two files, since we will only need one, let's remove the second one:
~~~
$ rm agalactiae_CJB111/GCF_015221735.2_ASM1522173v2_genomic.fna
~~~
{: .language-bash}

> ## Downloading specific formats
> In this example, we have downloaded the genome in `fasta` format. However, we can
>  use the `--format` or `-F` flags to get any other format of interest.
>  For example, the `gbk` format files (which contain information about the
>  coding sequences, their locus, the name of the protein and the full
>  nucleotide sequence of the assembly, and are useful for annotation double-checking)
>  can be downloaded by specifying our queries with `--format genbank`.
{: .callout}

> ## Exercise 1: Searching for desired strains
>  Until now we have downloaded only specific strains that we were looking for. Write a command that would tell you which genomes are 
>  available for all the *Streptococcus* genus.
>
> **Bonus**: Make a file with the output of your search.
> 
> > ## Solution
> > Use the `-n` flag to make it a dry run. 
> > Search only for the genus *Streptococcus* without using the strain flag.
> > ~~~
> > $  ncbi-genome-download -F fasta --genera "Streptococcus" -n bacteria
> > ~~~
> > {: .language-bash}
> > ~~~
> > Considering the following 18331 assemblies for download:
> > GCF_000959925.1	Streptococcus gordonii	G9B
> > GCF_000959965.1	Streptococcus gordonii	UB10712
> > GCF_000963345.1	Streptococcus gordonii	I141
> > GCF_000970665.2	Streptococcus gordonii	IE35
> > .
> > .
> > .
> > ~~~
> > {: .output}
> > **Bonus**: Redirect your command output to a file with the `>` command.
> > ~~~
> > $  ncbi-genome-download -F fasta --genera "Streptococcus" -n bacteria > ~/pan_workshop/data/streptococcus_available_genomes.txt
> > ~~~
> > {: .language-bash}
> > 
> {: .solution}
{: .challenge}



