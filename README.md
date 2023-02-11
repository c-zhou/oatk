# Oatk: an organelle genome assembly toolkit

## Overview
Oatk is designed for *de novo* assembly of complex plant organelle genomes using PacBio HiFi data. The toolkit consists of three major tools. `run_syncasm` is a *de novo* HiFi read assembler using a sparse de Bruijn graph constructed from closed syncmers ([Edgar, R. 2021](https://peerj.com/articles/10805/)). `hmm_annotation` is a HMMER wrapper for convenient annotation of organelle sequences using a pre-built HMM profile database which is available at [OatkDB](https://github.com/c-zhou/OatkDB.git). `path_finder` is a tool used for parsing and circularising organelle genomes from the assembled sequences combining the HMM annotations and assembly graph structure.

## Installation
You need to have a C compiler, GNU make and zlib development files installed. Download the source code from this repo or with `git clone https://github.com/c-zhou/oatk.git`. Then type `make` in the source code directory to compile.

## Run Oatk

A typical run consists of three steps: HiFi read assembly with `run_syncasm`, HMM annotation with `hmm_annotation` and organelle genome extraction with `path_finder`.

### HiFi read assembly

Here is an example to run `run_syncasm`,

    run_syncasm -k 1001 -c 150 -t 8 -o oatk $input_hifi

Four optional and one positional parameters are set in this example:

`-k` specifies the syncmer size.

`-c` specifies the synmcer coverage threshold. Syncmers with a coverage below the threshold will be excluded from the assembly graph. This is a crucial parameter for organelle genome assembly. A proper choice of this threshold will help filter out the nuclear genome and numts. **As an empirical instruction, this parameter can be set as 5-10 times the value of nuclear sequence coverage**.

`-t` specifies the thread number.

`-o` specifies the prefix of the output files.

The positional parameter specifies the input PacBio HiFi data file. The program recognises both FASTA and FASTQ format, plain and gzipped.

The command will generate several files. The `oatk_utg_final.gfa` is the file for the final genome assembly in [GFA format](https://github.com/GFA-spec/GFA-spec). 

### HMM annotation

Here is an example to run `hmm_annotation`,

    hmm_annotation -t 4 --nhmmscan /usr/bin/nhmmscan -o oatk_utg_final.mito.txt $hmm_db_dir/angiosperm_mito.fam oatk_utg_final.gfa
    hmm_annotation -t 4 --nhmmscan /usr/bin/nhmmscan -o oatk_utg_final.pltd.txt $hmm_db_dir/angiosperm_pltd.fam oatk_utg_final.gfa
 
Three optional and two positional parameters are set in this example:

`-t` specifies the thread number.

`-o` specifies the output file.

`--nhmmscan` specifies the path to the [nhmmscan](http://hmmer.org/) executable. If not specified, the program will assume it is in the environmental path.

The first positional parameter specifies the HMM profile database. Proper database file for the target species should be selected. In this example, the database for angiosperms was used. Here, 'mito' for mitochondrial and 'pltd' for chloroplast. 

The second positional parameter specifies the file of sequences for annotation. The program recognises FASTA, FASTQ, and GFA format, both plain and gzipped.

### Organelle genome extraction

Here is an example to run `path_finder`,

    path_finder -m oatk_utg_final.mito.txt -p oatk_utg_final.pltd.txt -o oatk oatk_utg_final.gfa
    
Three optional and one positional parameters are set in this example:

`-m` specifies the mitochondrial annotation file. With this parameter, the program will attempt to parse the mitochondrial genome.

`-p` specifies the chloroplast annotation file. With this parameter, the program will attempt to parse the chloroplast genome.

`-o` specifies the prefix of the output files.

The positional parameter specifies the assembly graph file.

The major output files in this example include: `oatk.mito.fasta` and `oatk.mito.gfa` for mitochondrial genome assembly; and `oatk.pltd.fasta` and `oatk.pltd.gfa` for chloroplast genome assembly. If the assembly is too complicated to solve, there will also be `oatk.mito.unassembled.fasta` and `oatk.pltd.unassembled.fasta` for the unitigs.

## Other auxiliary tools

* ***path_to_fasta*** is a tool used to extract FASTA sequences from GFA file with a path. For example, `path_to_fasta oatk_utg_final.gfa u1+,u2-,u3-,u2+`.

