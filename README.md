# Oatk: an organelle genome assembly toolkit [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7631376.svg)](https://doi.org/10.5281/zenodo.7631376)

## Overview

Oatk is designed for *de novo* assembly of complex plant organelle genomes using PacBio HiFi data. It can also be used to assemble other simple organelle genomes such as animal mitochondria. The toolkit consists of four major tools. `syncasm` is a *de novo* HiFi read assembler using a sparse de Bruijn graph constructed from closed syncmers ([Edgar, R. 2021](https://peerj.com/articles/10805/)). `hmm_annotation` is a [HMMER](http://hmmer.org/) wrapper for convenient annotation of organelle sequences using a pre-built HMM profile database which is available at [OatkDB](https://github.com/c-zhou/OatkDB.git). `pathfinder` is a tool used for parsing and circularising organelle genomes from the assembled sequences combining the HMM annotations and assembly graph structure. `oatk` is a wrapper for running `syncasm`, `hmm_annotation` and `pathfinder` collectively. There is also an auxiliary tool `path_to_fasta` used to extract FASTA sequences from a [GFA format](https://github.com/GFA-spec/GFA-spec) file given a path.

## Installation

You need to have a C compiler, GNU make and zlib development files installed. Download the source code from this repo or with `git clone https://github.com/c-zhou/oatk.git`. Then type `make` in the source code directory to compile.

## Dependencies

- You need the [nhmmscan](http://hmmer.org/) executable to run `oatk` and `hmm_annotation`. It could be made available to the program through the environmental path or specified as a program parameter.
- You also need the organelle gene profile database to run `oatk` and `hmm_annotation`. It could be downloaded from [OatkDB](https://github.com/c-zhou/OatkDB.git). There are many pre-built databases ready to use. You can also build your own database following the instructions on that page. For land plant organelles, while there are databases for several lower-rank taxonomies, the embryophyta is the recommended one to start with.

## Run Oatk

The easiest way is to use the `oatk` wrapper. Internally, it consecutively runs three modules: HiFi read assembly with `syncasm`, HMM annotation with `hmm_annotation` and organelle genome extraction with `pathfinder`. You can also run three components separately as detailed below.

The `oatk` program also allows you to extract organelle genomes from your own genome assembly graph built from such as [MBG](https://github.com/maickrau/MBG), [Spades](https://github.com/ablab/spades) or [Hifiasm](https://github.com/chhylp123/hifiasm). To do this, you need to include the `-G` option amd the input sequence file(s) needs to be replaced by a GFA file. It should be noted that `oatk` (and `pathfinder`) highly relies on the sequence and arc coverage (especially the sequence coverage) to solve the genome structure, so the overlap/string graphs constructed from many long-read assemblers (e.g., Hifiasm) are not ideal for complete organelle genome construction. Also, different assemblers use different tags for sequence and arc coverage in the GFA file. By default, `oatk` (and `pathfinder`) uses `KC:i`, `SC:i` and `EC:i` as the tag for kmer, sequence and arc coverage respectively, where the value of kmer coverage is approximately equal to the product of the sequence coverage and the sequence length. This is the default output format for `syncasm`. When you use GFA inputs from other assemblers, you may need to explicitly pass this information to `oatk` (and `pathfinder`) with `--kmer-c-tag` (or `--seq-c-tag`) and `--edge-c-tag` options. For example, with MBG output, you need to specify `--kmer-c-tag FC:f` and `--edge-c-tag ec:i`, otherwise, all coverage values will be considered as 1, and will probably lead to incorrect genome structures.

Below are examples of running these programs using the Arabidopsis thaliana organelle genome data download from this [Zenodo data repository](https://zenodo.org/records/10367917).

### Use oatk wrapper

Here is an example to run `oatk`,

    oatk -k 1001 -c 30 -t 8 --nhmmscan /bin/nhmmscan -m embryophyta_mito.fam -p embryophyta_pltd.fam -o ddAraThal4 ddAraThal4_organelle.hifi.fa.gz

Seven optional and one positional parameters are set in this example:

`-k` specifies the syncmer size (same as the default value).

`-c` specifies the synmcer coverage threshold. Syncmers with coverages below the threshold will be excluded from the assembly graph. This is a important parameter for organelle genome assembly. A proper choice of this threshold will help filter out the nuclear genome and numts. **As an empirical instruction, this parameter can be set as 5-10 times the value of nuclear sequence coverage**.

`-t` specifies the thread number.

`--nhmmscan` specifies the path to the [nhmmscan](http://hmmer.org/) executable. If not specified, the program will assume it is in the environmental path.

`-m` specifies the mitochondrion (MT) gene profile database. With this parameter, the program will attempt to parse the mitochondrion genome.

`-p` specifies the chloroplast (PT) gene profile database. With this parameter, the program will attempt to parse the chloroplast genome.

`-o` specifies the prefix of the output files.

The positional parameter specifies the input PacBio HiFi data file. The program recognises both FASTA and FASTQ format, plain and gzipped. Multiple input files are allowed.

Upon a successful run, the command will generate these major files:
~~~
ddAraThal4.utg.final.gfa     the GFA file for the final genome assembly             | syncasm
ddAraThal4.annot_mito.txt    the MT gene annotation file for assembled sequences    | hmm_annotation
ddAraThal4.annot_pltd.txt    the PT gene annotation file for assembled sequences    | hmm_annotation
ddAraThal4.mito.gfa          the subgraph for the MT genome                         | pathfinder
ddAraThal4.mito.bed          the gene annotation for the MT sequences               | pathfinder
ddAraThal4.mito.ctg.fasta    the structure-solved MT contigs                        | pathfinder
ddAraThal4.mito.ctg.bed      the genome annotation for MT contigs                   | pathfinder
ddAraThal4.pltd.gfa          the subgraph for the PT genome                         | pathfinder
ddAraThal4.pltd.bed          the gene annotation for the PT sequences               | pathfinder
ddAraThal4.pltd.ctg.fasta    the structure-solved PT contigs                        | pathfinder
ddAraThal4.pltd.ctg.bed      the genome annotation for PT contigs                   | pathfinder
~~~

### Use individual programs

These three programs share many parameters with the `oatk` wrapper. Unless specified, the parameters used this section is the same as those in the previous section.

#### 1. HiFi read assembly

Here is an example to run `syncasm`,

    syncasm -k 1001 -c 30 -t 8 -o ddAraThal4 ddAraThal4_organelle.hifi.fa.gz

The major file generated is `ddAraThal4.utg.final.gfa` for the GFA file of the final genome assembly. 

#### 2. HMM annotation

Here is an example to run `hmm_annotation`,

    hmm_annotation -t 8 --nhmmscan /bin/nhmmscan -o ddAraThal4.annot_mito.txt embryophyta_mito.fam ddAraThal4.utg.final.gfa
    hmm_annotation -t 8 --nhmmscan /bin/nhmmscan -o ddAraThal4.annot_pltd.txt embryophyta_pltd.fam ddAraThal4.utg.final.gfa
 
Three optional and two positional parameters are set in this example:

The `-t` and `--nhmmscan` options are the same as the `oatk`. The `-o` option set the output file name. 

The first positional parameter specifies the HMM profile database. Proper database file for the target species should be selected. In this example, the databases for embryophyta were used. Here, 'mito' for mitochondrion and 'pltd' for chloroplast. 

The second positional parameter specifies the file of sequences for annotation. Here we use the GFA file generated by `syncasm`. The program also recognises FASTA and FASTQ format, both plain and gzipped.

#### 3. Organelle genome extraction

Here is an example to run `pathfinder`,

    pathfinder -m ddAraThal4.annot_mito.txt -p ddAraThal4.annot_pltd.txt -o ddAraThal4 ddAraThal4.utg.final.gfa
    
Three optional and one positional parameters are set in this example:

`-m` specifies the mitochondrion annotation file. With this parameter, the program will attempt to parse the mitochondrion genome.

`-p` specifies the chloroplast annotation file. With this parameter, the program will attempt to parse the chloroplast genome.

`-o` specifies the prefix of the output files.

The positional parameter specifies the assembly graph file.

The major output files are the organelle subgraph files, contig FASTA files and sequence annotation files as listed under the `oatk` example.

## Other auxiliary tools

***path_to_fasta*** is a tool used to extract FASTA sequences from a GFA file with a path. For example, `path_to_fasta oatk_utg_final.gfa u1+,u2-,u3-,u2+`.

