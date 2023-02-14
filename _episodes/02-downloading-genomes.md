---
title: "Downloading Genomic Data"
teaching: 30
exercises: 15
questions:
- "How to download public genomic data from the command line?"
objectives:
- "Explore `ncbi-genome-download` as a tool for genomic data fetching from the NCBI."
keypoints:
- "The `ncbi-genome-download` package is a set of scripts to download genomes from the NCBI."
---

## Getting Genomic Data from the NCBI

The NCBI Genome Downloading Scripts provide a shell command that
allows users to download genomes from the NCBI. This [package](https://github.com/kblin/ncbi-genome-download)
is highly useful, it enables to specify queries as much as it is desired. It simplifies
the process of getting the data directly into the working directory. Firstly, activate
the `ncbi-genome-download` Conda environment.  

~~~
$ conda activate ncbi-genome-download
~~~
{: .language-bash}

~~~
(ncbi-genome-download) $
~~~
{: .output}
For practicality, the prompt is written only as `$` instead of `(ncbi-genome-download) $`.

The full list of parameters to incorporate in your downloads can be obtained by typing:
~~~
$ ncbi-genome-download --help
~~~
{: .language-bash}

This command outputs a long user manual; some of the most important parameters are:

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
it could be because you are in the `base` and not inside the `ncbi-genome-download`
Conda environment. Once inside the environment we are ready to use the package.
We'll now change to our data directory created in the previous section.

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
18RS21  antismash_db.csv  H36B
~~~
{: .output}

To construct the first pangenome, several _Streptococcus_ strains were considered.
Let's search the genome of _Streptoccocus agalactie 515_, one of the original strains.

Prior to downloading anything from the NCBI, it is advisable to verify if the information we seek
is available in the database and its content. Include
the `-n` flag within your command to do a "dry-run", which means to show the available accessions instead of downloading them.
Use the flag `--formats` to specify the desired format, `--genera` to specify the genus (or species), and `-S` for the strain. The `ncbi-genome-download` command must always end with the name of the group of organisms where the search will be performed, which in our case is `bacteria`.

~~~
$ ncbi-genome-download --formats fasta --genera "Streptococcus agalactiae" -S 515 -n bacteria 
~~~
{: .language-bash}
~~~
Considering the following 1 assemblies for download:
GCF_012593885.1 Streptococcus agalactiae 515    515
~~~
{: .output}

Once we know that the genome is available in NCBI,
let us download its corresponding fasta file. The flag 
`-o agalactie_515` specifies `agalactie_515` as the output directory. 
Notice that the `-n` flag is not included in this command. This is because
now we will downloading the genome instead of finding if it is available
in NCBI.

~~~
$ ncbi-genome-download --formats fasta --genera "Streptococcus agalactiae" -S 515 -o agalactiae_515 bacteria 
~~~
{: .language-bash}


Once the prompt `$` is available again, use the command `tree` 
to show the contents of the recently downloaded directory `agalactiae_515`.

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
which means is in a nucleotide `fasta` format. Let us move the file to the
`agalactiae_515` directory and remove the extra content that we will not 
use again in this lesson.

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

Right now, you have the _Streptococcus agalactiae 505_ genomic `fasta` file in a directory.
However, five more strains were included in Tettelin's first pangenome. It is the moment
to practice what you've learned about cycles in the shell lesson. Instead
of downloading the genomes one by one, we will write a `while` cycle. 
With the nano editor, create a file to add the other four strains that 
Tettelin included in the first pangenome. The missing strains are A909, 
COH1, CJB111, NEM316 and 2603V/R. Write one strain per line in the file and named it
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

First, let us read the lines of Tettelin file, and print
them in the terminal with the `echo strain $line` command.  
`strain` is just a word that we will print, and `$line` will 
store the value of each of the lines of the `Tettelin.txt` file.

~~~
$ cat TettelinList.txt | while read line 
> do echo strain $line
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
strain 2603V/R
~~~
{: .output}

The tool has successfully found the three strain. Notice that
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

Finally, we need to move the genome files to their corresponding locations and
get rid of unnecessary directories. To do so, we'll use a `while` cycle as follows:

~~~
 $ cat TettelinList.txt | while read line
 > do 
 > echo removing refseq directory of strain $line
 > mv agalactiae_$line/refseq/bacteria/*/*.fna agalactiae_$line/. # Move file to current directory 
 > rm -r agalactiae_$line/refseq  # Remove refseq directory
 > done
 ~~~
{: .language-bash}
~~~
removing refseq directory of strain A909
removing refseq directory of strain COH1
removing refseq directory of strain CJB111
removing refseq directory of strain NEM316
removing refseq directory of strain 2603V/R
~~~
{: .output}

At this point, you should have four directories starting with `agalactiae_` containing
the following:

~~~
$ tree agalactiae_*
~~~
{: .language-bash}

~~~
agalactiae_515
└── GCF_012593885.1_ASM1259388v1_genomic.fna
agalactiae_A909
└── GCF_000012705.1_ASM1270v1_genomic.fna
agalactiae_CJB111
├── GCF_000167755.1_ASM16775v1_genomic.fna
└── GCF_015221735.2_ASM1522173v2_genomic.fna
agalactiae_COH1
└── GCF_000689235.1_GBCO_p1_genomic.fna
~~~
{: .output}


> ## Downloading specific formats
> In this example, we have downloaded the genome in `fasta` format. However, we can
>  use the `--format` or `-F` flags to get any other format of interest.
>  For example, the `gbk` format files (which contain information about the
>  coding sequences, their locus, the name of the protein and the full
>  nucleotide sequence of the assembly, and are useful for annotation double-checking)
>  can be downloaded by specifying our queries with `--format genbank`.
{: .callout}

> ## Exercise 1: Downloading data from NCBI using the command line FIXME
>
> To download the genome of another strain of *S. thermophilus* we need to perform the following steps:
>
> 1. Download the genome of *Streptococcus thermophilus* with the NCBI assembly number `GCF_000011825.1` in a `fasta` format and save it to an output directory named `thermophilusLMG18311`.
> 2. Change to the directory where the `fasta` file is located and unzip it.
> 3. Move the `fasta` file all the way back to the `thermophilusLMG18311` directory.
> 4. Change to the `thermophilusLMG18311` directory and delete the `refseq` subdirectory created by the downloading tool.
>
> Complete the following set of commands to carry out the previous steps:
>
> Step 1.
>
> ~~~
> $ ncbi-genome-download -F __________ --genera __________ -A __________ -o __________ bacteria
> ~~~
> {: .language-bash}
>
> Step 2.
>
> ~~~
> $ cd __________/refseq/bacteria/GCF_000011825.1/
> $ __________ GCF_000011825.1_ASM1182v1_genomic.fna.gz
> ~~~
> {: .language-bash}
>
> Step 3.
>
> ~~~
> $ mv GCF_000011825.1_ASM1182v1_genomic.fna __________
> ~~~
> {: .language-bash}
>
> Step 4.
>
> ~~~
> $ cd __________
> $ rm -rf refseq
> ~~~
> {: .language-bash}
>
> > ## Solution
> >
> > Step 1. Using the provided information, we complete each blank space with the corresponding word:
> >
> > ~~~
> > $ ncbi-genome-download -F fasta --genera "Streptococcus thermophilus" -A GCF_000011825.1 -o thermophilusLMG18311 bacteria
> > ~~~
> > {: .language-bash}
> >
> > Step 2. The previous command creates a `thermophilusLMG18311` subdirectory. To get to the `fasta` file we must go through a sequence of subdirectories and then apply the `gunzip` command to unzip the `fasta` file:
> >
> > ~~~
> > $ cd thermophilusLMG18311/refseq/bacteria/GCF_000011825.1/
> > $ gunzip GCF_000011825.1_ASM1182v1_genomic.fna.gz
> > ~~~
> > {: .language-bash}
> >
> > Step 3. We have now an unzipped `fasta` file. The parent directory `thermophilusLMG18311` is located three directories above the current one. To get to the first parent directory, you would type `..`; if you want to get to the second parent directory, you would use `../..`. Thus, to move the `fasta` file to the `thermophilusLMG18311` directory, we need to type:
> >
> > ~~~
> > $ mv GCF_000011825.1_ASM1182v1_genomic.fna ../../..
> > ~~~
> > {: .language-bash}
> >
> > Step 4. Finally, we move back to the `thermophilusLMG18311` directory (in a similar manner as in the previous step) and delete the `refseq` directory.
> >
> > ~~~
> > $ cd ../../..
> > $ rm -r refseq
> > ~~~
> > {: .language-bash}
> >
> {: .solution}
{: .challenge}


Finally, it is a good practice
to keep your raw data untouched, taking this into account you can remove the writing
permission of the data directory.  

~~~
$ cd ~/pan_workshop/data/
$ chmod -w ~/pan_workshop/data/*
$ ls -lh ~/pan_workshop/data/
~~~
{: .language-bash}

~~~
dr-xr-xr-x 2 alumno17 alumno17 4.0K feb  9 18:27 18RS21
dr-xr-xr-x 2 alumno17 alumno17 4.0K feb  9 18:27 agalactiae_515
dr-xr-xr-x 2 alumno17 alumno17 4.0K feb  9 18:27 agalactiae_A909
dr-xr-xr-x 2 alumno17 alumno17 4.0K feb  9 18:27 agalactiae_CJB111
dr-xr-xr-x 2 alumno17 alumno17 4.0K feb  9 18:27 agalactiae_COH1
-r--r--r-- 1 alumno17 alumno17 1.4M jun 12  2022 antismash_db.csv
dr-xr-xr-x 2 alumno17 alumno17 4.0K feb  9 18:27 H36B
-r--r--r-- 1 alumno17 alumno17   17 feb  9 18:03 TettelinList.txt
dr-xr-xr-x 2 alumno17 alumno17 4.0K feb  9 18:27 thermophilusLMG18311  
~~~
{: .output}

> ## Exercise 2: Unknown strains
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
> > $  ncbi-genome-download -F fasta --genera "Streptococcus" -n bacteria > streptococcus_available_genomes.txt
> > ~~~
> > {: .language-bash}
> > 
> {: .solution}
{: .challenge}



