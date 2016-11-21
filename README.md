
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
./INSTALL.sh

```


##Dependencies
PiReT run require fowllowing dependencies and in your path.

### Unix
- sed
- awk
- find


### Third party softwares/packages must be installed and in the path
- [samtools](http://www.htslib.org)
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [bwa](http://bio-bwa.sourceforge.net)
- [R](https://www.r-project.org)
- [HTseq](http://www-huber.embl.de/HTSeq/doc/overview.html)
- [HiSat](https://ccb.jhu.edu/software/hisat/index.shtml)
- [Python 2.5 or higher](https://www.python.org/download/releases/2.7/)

### R packages
- [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
- [DEseq](http://bioconductor.org/packages/release/bioc/html/DESeq.html)

### Python packages
- [numpy](http://www.numpy.org)
- [matplotlib](http://matplotlib.org)


## pipeline options

If you have


```
    perl runPipeline_rRNA_noqsub_commandline.pl [options] -exp exp_descriptfile.txt -d workdir -prokaryote_fasta indexprokaryote.fa -eukarya_fasta indexeukarya.fa -index_ref_bt2 indexfile -gff_prokaryote prokaryote.gff -gene_coverage_ref gene_coverage_reference.fa
```

`-d`: string, working directory where all output files will be written, the user must have write permission





## Citations:
If you use PiReT please cite following papers:

- Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. [PMID: 19505943]
- Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature methods, 9(4), 357-359. [PMID: 22388286]
- Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168]
- Anders S and Huber W (2010). Differential expression analysis for sequence count data. Genome Biology, 11, pp. R106. [PMID: 20979621]
- McCarthy, J. D, Chen, Yunshun, Smyth and K. G (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. Nucleic Acids Research, 40(10), pp. -9. [PMID: 22287627]
- Anders, S., Pyl, P. T., & Huber, W. (2014). HTSeq–a Python framework to work with high-throughput sequencing data. Bioinformatics. [PMID: 25260700]
- Kim, D., Langmead, B., & Salzberg, S. L. (2015). HISAT: a fast spliced aligner with low memory requirements. Nature methods, 12(4), 357-360. [PMID: 25751142]





