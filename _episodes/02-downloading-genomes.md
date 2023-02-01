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
For practicality the prompt is written only as `$` instead of `(ncbi-genome-download) $`.

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
it could be because you are in `base` and not inside the `ncbi-genome-download`
Conda environment. Once inside the environment we are ready to use the package.
Next, we need to go to our data directory.
~~~
$ cd ~/gm_workshop/data
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
18RS21/  515/  A909/  COH1/  CJB111/  H36B/
~~~
{: .output}

Prior to downloading anything from the NCBI, it is advisable to verify if the information we seek
is available in the database and its content. Include
the `-n` flag within your command to avoid this information to be downloaded. For instance, if we wish to check the availability of the genome of the
LMD-9 strain of the *Streptococcus thermophilus* bacterium in FASTA format, we would type the following:

~~~
$ ncbi-genome-download --formats fasta --genera "Streptococcus thermophilus" -S LMD-9 -n bacteria
~~~
{: .language-bash}

~~~
Considering the following 1 assemblies for download:
GCF_000014485.1 Streptococcus thermophilus LMD-9    	LMD-9
~~~
{: .output}

As you can see, there is one assembly available assigned to the number `GCF_000014485.1`.
We proceed to download it to an output directory named `thermophilusLMD9`:

~~~
$ ncbi-genome-download --formats fasta --genera "Streptococcus thermophilus" -S LMD-9 -o thermophilusLMD9 bacteria
~~~
{: .language-bash}

This script downloads a compressed FASTA file into a specific set of directories:

~~~
$ tree thermophilusLMD9/
~~~
{: .language-bash}

~~~
thermophilusLMD9/
└── refseq
	└── bacteria
    	└── GCF_000014485.1
        	├── GCF_000014485.1_ASM1448v1_genomic.fna.gz
        	└── MD5SUMS
~~~
{: .output}

In order to view it, we must decompress it using `gunzip`:

~~~
$ gunzip thermophilusLMD9/refseq/bacteria/GCF_000014485.1/GCF_000014485.1_ASM1448v1_genomic.fna.gz
~~~
{: .language-bash}

Now we can explore the file and move it to the main directory `thermophilusLMD9`, and delete the `refseq` directory as it is not longer needed:

```
$ mv thermophilusLMD9/refseq/bacteria/GCF_000014485.1/GCF_000014485.1_ASM1448v1_genomic.fna thermophilusLMD9/
$ rm -r thermophilusLMD9/refseq
```
{: .language-bash}


> ## Downloading specific formats
> In this example, we have downloaded the genome in FASTA format. However, we can
>  use the `--format` or `-F` flags to get any other format of interest.
>  For example, the `gbk` format files (which contain information about the
>  coding sequences, their locus, the name of the protein and the full
>  nucleotide sequence of the assembly, and are useful for annotation double-checking)
>  can be downloaded by specifying our queries with `--format genbank`.
{: .callout}

> ## Exercise 1: Downloading data from NCBI using the command line
>
> To download the genome of another strain of *S. thermophilus* we need to perform the following steps:
>
> 1. Download the genome of *Streptococcus thermophilus* with the NCBI assembly number `GCF_000011825.1` in a FASTA format and save it to an output directory named `thermophilusLMG18311`.
> 2. Change to the directory where the FASTA file is located and unzip it.
> 3. Move the FASTA file all the way back to the `thermophilusLMG18311` directory.
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
> > Step 2. The previous command creates a `thermophilusLMG18311` subdirectory. To get to the FASTA file we must go through a sequence of subdirectories and then apply the `gunzip` command to unzip the FASTA file:
> >
> > ~~~
> > $ cd thermophilusLMG18311/refseq/bacteria/GCF_000011825.1/
> > $ gunzip GCF_000011825.1_ASM1182v1_genomic.fna.gz
> > ~~~
> > {: .language-bash}
> >
> > Step 3. We have now an unzipped FASTA file. The parent directory `thermophilusLMG18311` is located three directories above the current one. To get to the first parent directory, you would type `..`; if you want to get to the second parent directory, you would use `../..`. Thus, to move the FASTA file to the `thermophilusLMG18311` directory, we need to type:
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

Make sure you have downloaded both strains of *S. thermophilus*, LMD-9 and LMG 18311,
as they will be required in upcoming episodes.   

Finally, it is a good practice
to keep your raw data untouched, taking this into account you can remove the writing
permission of the data directory.  

~~~
$ cd ~/gm_workshop/data/
$ chmod -w ~/gm_workshop/data/
$ ls -lh ~/gm_workshop/data/
~~~
{: .language-bash}

~~~
drwxr-xr-x 2 alumno17 alumno17	4096 Jun  6 13:20 18RS21  
drwxr-xr-x 2 alumno17 alumno17	4096 Jun  6 13:21 515  
drwxr-xr-x 2 alumno17 alumno17	4096 Jun  6 13:21 A909  
-rw-r--r-- 1 alumno17 alumno17 1460118 Jun 12 12:20 antismash_db.csv  
drwxr-xr-x 2 alumno17 alumno17	4096 Jun  6 13:21 CJB111  
drwxr-xr-x 2 alumno17 alumno17	4096 Jun  6 13:21 COH1  
drwxr-xr-x 2 alumno17 alumno17	4096 Jun  6 13:21 H36B  
drwxr-xr-x 2 alumno17 alumno17	4096 Jun 13 15:25 thermophilusLMD9  
drwxr-xr-x 2 alumno17 alumno17	4096 Jun 13 15:26 thermophilusLMG18311  
~~~
{: .output}


> ## Discussion
>  Until now we have downloaded only one file at a time. How would you look for several
>  genomes at the same time?
> > ## Solution
> > One idea would be to look for all the genomes in one genera such as _Streptococcus_
> > ~~~
> > $  ncbi-genome-download -F fasta --genera "Streptococcus " -n bacteria
> > ~~~
> > {: .language-bash}
> >
> {: .solution}
{: .challenge}



