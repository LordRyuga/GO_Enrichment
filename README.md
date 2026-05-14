# Understanding biological functions from a list of genes

Many curiosities or analysis boils down to a potential list of candidates. For
example, in cases of genome sequence or genome-wide assays, a set of analysis
leads to a list of genes to investigate. A natural next step is to first look
at in what set of biological functions they are involved. To do this, all we
need to do is to feed that list to existing R or Python packages and get a
visualization.

In this tiny project, we are going to simply do that. We will get a list of
genes after searching for a sequence of interest in promoter region of all the
genes and then look at their biological functions.

You will need the following datasets to do the project.

- [Human genome](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/), search for `hg38.fa.gz` and download

- Dreg from [EMBOSS](https://anaconda.org/bioconda/emboss)

- [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) Bioconductor Package

- [bedtools](https://anaconda.org/channels/bioconda/packages/bedtools/overview) 

## Environment setup

If you have a Windows system, install Windows Subsystem for Linux, WSL. It will give you direct access to Linux system. Students with either a Linux or MacOS system will have easy time setting up the environment to execute the project. 

Setup [Miniforge](https://github.com/conda-forge/miniforge), and create a virtual environment with the following command:

```
mamba create -n go_enrichment python=3.12 
```

Once the above run successfully, you can do the following. 

```
mamba activate go_enrichment
```

```
mamba install bioconda::emboss
```

## Target sequences to search in the promoter region

As explained in the class, there about 1700 Transcription Factors in humans categorized in different families. You can know more about them [here](https://www.sciencedirect.com/science/article/pii/S0092867418301065) if you are interested. For the assignment, we will have one example sequence "GCGC..GCGC" (the dot represents any nucleotide). It is a recognition sequence for NRF1 Transcription Factor, a well known factor for regulating genes for Mitochondria Biogenesis, and a host of other biological functions. We are going to explore TSS upstream sequences and perform gene enrichment analysis. 

I hope you all have downloaded the human genome, and also have a working version of WSL/Linux/MacOS. Make sure the data is accessible in WSL (not in powershell). 

### Step 1

Convert human gene annotation file (`human_gene_annotation.tsv.gz`) into a bed file. Here, I am sharing a few lines of the bed file.

```
chr1    11869   11870   chr1@11869-11870|DDX11L2        .       +
chr1    12010   12011   chr1@12010-12011|DDX11L1        .       +
chr1    17436   17437   chr1@17436-17437|MIR6859-1      .       -
chr1    24886   24887   chr1@24886-24887|WASH7P .       -
```

Column description are below:

- Col1: Chromosome

- Col2: Transcription Start Site (TSS) 

- Col3: Col2 + 1

- Col4: Col1@Col2-Col3|<gene name>

- Col5: "." (place holder for score)

- Col6: Strand (+/-) 


Once you have the bed file - your goal is to extend the coordinates (Col2 or Col3) by 500 bases in a strand aware manner. Install `bedtools` and look at the `bedtools slop` command


In the `go_enrichment` environment, install `bedtools`. 

