
#PiReT

Pipeline for Reference based Transcriptomics.

[![Build Status](https://travis-ci.com/mshakya/PiReT.svg?token=xwcWcg2wroskmENQQapz&branch=master)](https://travis-ci.com/mshakya/PiReT)

## Installing PiReT

Download the latest version of PiReT from [github](https://github.com/mshakya/PiReT.git) into a
Linux server.

```
git clone https://github.com/mshakya/PiReT.git
```

`cd` into the `PiReT` directory and

```
cd PiReT
```

There are two ways to install the required dependencies. We recommened using the `bioconda_INSTALL.sh`.

### Manual installation of dependencies.

```
./INSTALL.sh

```

### Bioconda installation of dependecies

PiReT uses many of the biofinformatic tools that are already available in [bioconda](https://bioconda.github.io). Many of these dependencies (see below) are likely already available in your UNIX system. The following script checks if the required dependencies (see below) are in your path and installs (download binaries within the PiReT directory and adds the path to your `~/.bashrc` or `~/.bash_profile`) it.

```
./bioconda_INSTALL.SH
```





##Dependencies
PiReT run require fowllowing dependencies and in your path.

### Programming/Scritpting languages
- [Python 2.7](https://www.python.org/downloads/release/python-2712/)
    - The pipeline does not work in Python 3.0 or higher yet due
- [Perl v5.16.3](https://www.perl.org/get.html)
    - The pipeline has only been tested in v5.16.3
- [R v3.3.1](https://www.r-project.org)


### Unix
- sed
- awk
- find
- curl/wget

### Installation dependecies
- [conda](http://conda.pydata.org/docs/index.html)
    If conda is not installed, `bioconda_INSTALL.sh` will download and install [miniconda](http://conda.pydata.org/miniconda.html), a "mini" version of `conda` that only install handful of packages compared to [anaconda](https://docs.continuum.io/anaconda/pkg-docs)

### Third party softwares/packages
- [jellyfish](http://www.genome.umd.edu/jellyfish.html)
- [samtools](http://www.htslib.org)
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [bwa](http://bio-bwa.sourceforge.net)
- [R](https://www.r-project.org)
- [HTseq](http://www-huber.embl.de/HTSeq/doc/overview.html)
- [HiSat2](https://ccb.jhu.edu/software/hisat/index.shtml)
- [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)

### R packages
- [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
- [DEseq](http://bioconductor.org/packages/release/bioc/html/DESeq.html)

### Python packages
- [numpy](http://www.numpy.org)
- [matplotlib](http://matplotlib.org)

### Perl modules
- [Parallel::ForkManager](http://search.cpan.org/~yanick/Parallel-ForkManager-1.19/lib/Parallel/ForkManager.pm)
- [String::Approx](http://search.cpan.org/dist/String-Approx/Approx.pm)

## How to run the pipeline?


```
    perl runPipeline_rRNA_noqsub_commandline.pl [options] -exp exp_descriptfile.txt -d workdir -prokaryote_fasta indexprokaryote.fa -eukarya_fasta indexeukarya.fa -index_ref_bt2 indexfile -gff_prokaryote prokaryote.gff -gene_coverage_ref gene_coverage_reference.fa
```

`-d`: working directory where all output files will be written, the user must have write permission.

`-prokaryote_fasta`: comma-separated list of referecnce genome (prokarya) fasta files (for making bowtie2 mapping index file). [optional]

`-gff_prokaryote`: comma-separated list of gff files for corresponding referecnce genome fasta files (contig names must match reference sequence header). [optional]

`-eukarya_fasta` : comma-separated list of referecnce genome (eukarya) fasta files (for making bowtie2 mapping index file). [optional]

`-gff_eukarya`: comma-separated list of gff files for corresponding referecnce genome fasta files (contig names must match reference sequence header). [optional]

`-index_ref_bt2`: bowtie2 mapping index file,  if the file exists, pipeline skips this step. [optional]

`-gene_coverage_fasta`: fasta file  (for directional coverage analysis, sequnce  must be part of prokaryote mapping reference sequence). [optional]

`-test_kingdom`: desired differential gene expression kingdom (both (for both eukarya and prokaryote), prokaryote, or eukarya (default:`prokaryote`));

`-test_method`: dessired differential gene expression method (both (for both EdgeR and Deseq2 method), EdgeR, or Deseq (default: `both`)); must have have at least 3 duplicates if using Deseq2.

`-cpu`: number of cpu to be used (default 1)

`-BAM_ready`: if mapping file are provided for samples by users (`yes` or `no`). default: `no`

`-significant_pvalue`: floating number cutoff to define significant differentially express genes, (default =0.001)

`-exp`: A tab delimited file that contains at least 3 columns with following header `ID`, `Rawread_files`, and  `group`. `Rawread_files` must have an absolute path.

`-pair_comparison`: tab delimited txt file descripting pairwise comparison. If the file is not specified, all possible pairwise analysis will be condected.


## Citations:
If you use PiReT please cite following papers:

- Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. [PMID: 19505943]
- Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature methods, 9(4), 357-359. [PMID: 22388286]
- Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168]
- Anders S and Huber W (2010). Differential expression analysis for sequence count data. Genome Biology, 11, pp. R106. [PMID: 20979621]
- McCarthy, J. D, Chen, Yunshun, Smyth and K. G (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. Nucleic Acids Research, 40(10), pp. -9. [PMID: 22287627]
- Anders, S., Pyl, P. T., & Huber, W. (2014). HTSeq–a Python framework to work with high-throughput sequencing data. Bioinformatics. [PMID: 25260700]
- Kim, D., Langmead, B., & Salzberg, S. L. (2015). HISAT: a fast spliced aligner with low memory requirements. Nature methods, 12(4), 357-360. [PMID: 25751142]
- Quinlan AR and Hall IM, 2010. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 26, 6, pp. 841–842. [PMID: 20110278]





