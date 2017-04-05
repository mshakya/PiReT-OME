
#PiReT

Pipeline for Reference based Transcriptomics - for VIOME.

[![Build Status](https://travis-ci.com/mshakya/PiReT.svg?token=xwcWcg2wroskmENQQapz&branch=master)](https://travis-ci.com/mshakya/PiReT)
[![codecov](https://codecov.io/gh/mshakya/PiReT/branch/master/graph/badge.svg?token=B0PzdxRcdj)](https://codecov.io/gh/mshakya/PiReT)

#Overview

## Requirements
Before installing PiReT-OME, user must have following installed. These normally come pre-installed as part of a unix machine. 
- curl/wget
- make
- git

## Installing PiReT-OME
Use `git clone` from command line.

```
git clone https://github.com/mshakya/PiReT.git
```

`cd` into the `PiReT` directory

```
cd PiReT
./INSTALL-PiReTOME.sh
```

PiReT-OME is a wrapper of RNA seq tools, many of which are available in [bioconda](https://bioconda.github.io). For installing `PiReT` we have provided a script `INSTALL-PiReTOME.sh` that first checks for required dependencies (including their versions) are installed and in your `PATH`. If not found, dependencies will be installed in `thirdParty` directories within `PiReT`. `sudo` privileges are not needed for installation. A log of all installation can be found in `install.log`

##Test
We have provided test data set to check if the installation was successful or not. `fastq` files can be found in `test_data/fastqs` and corresponding reference fasta files are found in `test_data/data`. To run the test, from within `PiReT` directory:

```
cd test_data

# if you are in a LINUX system:
sh ./test_pipeline_linux.sh

# if you are in Mac OS X:
sh ./test_pipeline_MacOSX.sh
```
These shell script automatically creates `experimental_design.txt` and runs the pipeline. All, but index files are written in `working directory`. See section below for details on output.

Pipeline run status can be checked in either `process.log` or `error.log`.


## Installed Dependencies
PiReT-OME requires following dependencies, all of which should be installed and added to PATH. All of the dependencies should be installed by `INSTALL-PiReTOME.sh`. However, some manual work may require depending on your system configuration.

### Programming/Scripting languages
- [Python >=v2.7](https://www.python.org/downloads/release/python-2712/)
    - The pipeline is not compatible with Python v3.0 or higher.
- [Perl >=v5.22.0](https://www.perl.org/get.html)
    - The pipeline has only been tested in v5.16.3 and v5.22.0
- [R >=v3.3.1](https://www.r-project.org)


### Required dependencies
This is the core list of dependencies. However, there are secondary dependencies for many of the listed tools, which will also be installed by `bioconda`.
- [conda v4.3.11](http://conda.pydata.org/docs/index.html)
    If [anaconda](https://www.continuum.io/downloads) is not installed, `INSTALL-PiReTOME.sh` will download and install [miniconda](http://conda.pydata.org/miniconda.html), a "mini" version of `conda` that only installs handful of packages compared to [anaconda](https://docs.continuum.io/anaconda/pkg-docs). However, we strongly recommend that users install anaconda as it comes with many packages that one will eventually use for other bioinformatic analyses.

### Third party softwares/packages
- [jellyfish (v2.2.6)](http://www.genome.umd.edu/jellyfish.html)
    - It is used by our fastQC program.
- [samtools (v1.3.1)](http://www.htslib.org)
    - It is used to parse (arrange, fetch, convert, etc.) the mapped `SAM/BAM` files.
- [HiSat2 (v2.0.5)](https://ccb.jhu.edu/software/hisat/index.shtml)
    - It is used to align reads to the reference. Before alignment, reference fasta is used for creating indices. 
- [bedtools (v2.26.0)](http://bedtools.readthedocs.io/en/latest/index.html)
    - It is used to detect allelic variants.
- [gffread (v0.9.6)](http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread_dl)
    - It is used to parse gff and gtf files.
- [StringTie (v1.3.3)](https://ccb.jhu.edu/software/stringtie/index.shtml)
    - It is used to quantify mapped reads, calculate TPM/FPKM, and count coverage.

### Perl modules
- [Parallel::ForkManager (v1.17)](http://search.cpan.org/~yanick/Parallel-ForkManager-1.19/lib/Parallel/ForkManager.pm)
- [String::Approx (v3.27)](http://search.cpan.org/dist/String-Approx/Approx.pm)

## Running PiReT

The pipeline can be run in a multiprocessor server with the ability to submit jobs in a queue system through `qsub` or in a single processor system where all jobs are run sequentially. The former requires a `qsub` system. Also, the current state of pipeline only works for a single processor.


```
    perl runPiReT-OME -exp experimental_design.txt \
    -d working_dir \
    -eukarya_fasta data/eukarya_test.fa -index_ref_ht2 RefIndex \
    -gff_eukarya data/eukarya_test.gff3 -gff_prokaryote data/test_prok.gff \
    -test_method both -gene_coverage_fasta data/test_prok.fa
```

`-d`: working directory where all output files/directories will be written, users must have write permission.

`-eukarya_fasta` : comma-separated list of reference genome (eukarya) fasta files. [optional]

`-gff_eukarya`: comma-separated list of gff files for corresponding reference genome fasta files (contig names must match reference sequence header). [optional]

`-index_ref_ht2`: HISAT2 mapping index file, if the file exists, pipeline skips this step. [optional]

<!-- `-cpu`: number of CPU to be used (default 1) -->

`-exp`: A tab delimited file that contains at least 3 columns with following header `ID`, `Rawread_files`, and  `group`. `Rawread_files` must have an absolute path.

## What is in the working directory (-d)?

Here are the list of directories that you should expect in `working directory` from test run
in `test_data`.

```
ls -R | grep ":$" | sed -e 's/:$//' -e 's/[^-][^\/]*\//--/g' -e 's/^/   /' -e 's/-/|/'

   |---eukarya
   |-----eukarya_test3
   |-logdir
   |-merged_gtfs
   |---samp1
   |---samp2
   |---samp3
   |---samp4
   |---samp5
   |---samp6
   |-samp1
   |---mapping_results
   |---trimming_results
   |-samp2
   |---mapping_results
   |---trimming_results
   |-samp3
   |---mapping_results
   |---trimming_results
   |-samp4
   |---mapping_results
   |---trimming_results
   |-samp5
   |---mapping_results
   |---trimming_results
   |-samp6
   |---mapping_results
   |---trimming_results
   |-sum_gene_count
   |---read_count
   |-----eukarya
   |-------eukarya_test3
   |---tmp_count
   |-----eukarya
   |-------eukarya_test3
   |---process.log
   |---error.log
   |---splice_sites_gff.txt
   |---eukarya.fa.fai
   |---eukarya.gtf
   |---

```

- `splice_sites_gff.txt`: contains known splice sites, generated using `scripts/extract_splice_sites.py`, a python script provided with *HISAT*.

`process.log`: report of all the commands/scripts/ that were ran as part of the pipeline.

`error.log`: any error are reported here.

`samp1-samp6`: The name of this directory corresponds to sample name (as given in experimental design file). Within this folder there are two sub-folders:

- `trimming_results`
            This folder contains results of quality trimming or filtering. This folder was generated using the same script that filters reads in [EDGE](https://bioedge.lanl.gov/edge_ui/) pipeline. It contains following files:
            
            - `fastqCount.txt`: summary (read length, count and nt content) of input fastq.
            - `samp1.1.trimmed.fastq`: QCed read1s of a pair.
            - `samp1.2.trimmed.fastq`: QCed read2s of a pair.
            - `samp1_qc_report.pdf: summary of QC run with figures.
            - `samp1.stats.txt: summary of QC run.
            - `samp1.unpaired.trimmed.fastq : unpaired reads from the QC.

- `mapping_results`
    This folder contains reads mapped using *HISAT2* in following formats. I
        - `mapped.sam`: output of *HISAT2*.
        - `samp1.alns.sorted.bam`: ordered `mapped.sam` file.
        - `mapped.log`: Alignment summary file from `HISAT2`
        - `samp1.summary.gtf`: A `GTF` style file with information of coverage (`cov`) , FPKM (`FPKM`), TPM (`TPM`), gene id (`ref_gene_id`), etc. as `attributes`. The file is the output (-o) of `StringTie` quantifying reads.
        - `samp1_full_coverage.gtf`: A `GTF` file that contains all transcripts that are fully covered by reads. See (`-C`) option of `StringTie` for details.
        - `samp1_gene_abundance.tab`: A tab delimited file with information on gene coverage, FPKM, and TPM. See (-A) opetion of `StringTie` for details.


- `eukarya.fa`: A copy of input eukarya `fna` files.

- `eukarya.fa.fai`: Indexed reference sequence from `eukarya.fa` using `samtools faidx`. A four column table with NAME, LENGTH, OFFSET, LINEBASES, and LINEWIDTH 

- `eukarya.gtf`: A `GTF` file from input `GFF` files.

- `eukarya.genedesc.txt`: A tab delimited file with gene id and function.

- `eukarya.gff`: A `GFF` file from input GFF files.

## Removing PiReT

For removal, since all dependencies that are not in your system are installed in `PiReT`, delete (`rm -rf`) `PiReT` folder is sufficient to uninstall the package. **Before removing check if your project files are within `PiReT` directory**.


##Contributions
- Migun Shakya
- Shihai Feng

## Citations:
If you use PiReT please cite following papers:

- **samtools**: Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. [PMID: 19505943]
- **HISAT2**: Kim, D., Langmead, B., & Salzberg, S. L. (2015). HISAT: a fast spliced aligner with low memory requirements. Nature methods, 12(4), 357-360. [PMID: 25751142]
- **BEDTools**: Quinlan AR and Hall IM, 2010. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 26, 6, pp. 841–842. [PMID: 20110278]
- *StringTie*: Pertea M, Pertea GM, Antonescu CM, Chang TC, Mendell JT  & Salzberg SL. StringTie enables improved reconstruction of a transcriptome from RNA-seq reads Nature Biotechnology 2015, doi:10.1038/nbt.3122 [PMID: 25690850]
