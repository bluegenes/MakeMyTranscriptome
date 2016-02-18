###MakeMyTranscriptome (MMT): an automated pipeline for *de novo* transcriptomics

[![Join the chat at https://gitter.im/bluegenes/MakeMyTranscriptome](https://badges.gitter.im/bluegenes/MakeMyTranscriptome.svg)](https://gitter.im/bluegenes/MakeMyTranscriptome?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

For a full guide, check out the [Wiki](http://github.com/bluegenes/makeMyTranscriptome/wiki)

MMT is an attempt to automate a current best-practices transcriptome pipeline from start (read quality control) to finish (annotated transcriptome). If you provide metadata in a [csv file](https://github.com/bluegenes/MakeMyTranscriptome/wiki/CSV_Input), we will also do some basic pairwise differential expression with DESeq2. Since DE is not usually one-size fits-all, we recommend you look into these results and rerun as necessary. 

MMT is composed of four principle tools: 
- **assembly** - assemble input reads into a transcriptome
- **quality** - assess the quality of a transcriptome using common length-based and annotation-based metrics
- **annotation** - annotate a transcriptome using ORF-prediction and BLAST-style mappings
- **expression**  - identify differentially expressed transcripts and genes 

A fifth tool is currently being developed:
- **compare** - compare MMT assemblies constructed with different parameters

Each tool can be run individually, or all can be executed as a single virtual tool, **full**. All tools make use of two additional convenience modules, **databases** and **tools** which check for required databases and tools, respectively. If you're working on a Linux machine, MMT can download the tools for you; otherwise, it will assess the tools you have available and print recommended installation instructions.

One important note: if you already have a transcriptome, the **full** module can still be used to perform all downstream analyses. If you provide your transcriptome with the ```--assembly``` option, MMT will simply skip the assembly steps (currently this means skipping read QC as well, but this can be addressed if there's interest). Similarly, if you do not provide a csv file, MMT will skip the expression steps.  


### Quick Start Guide:

MMT is a convenience tool that runs a number of bioinformatics tools in an automated fashion.

First, download MMT via git (below) or click the "Download ZIP" button.
```
git clone https://github.com/bluegenes/MakeMyTranscriptome.git
```
cd into the MakeMyTranscriptome directory.
```
cd MakeMyTranscriptome
```

To install tools and databases, run:

``` 
mmt setup --install --email youremail
``` 
If you just want to check the installed tools, omit the ```--install``` option. In release v0.1, databases will be downloaded in setup, but we'll update in a future release to allow you to specify a database path if you've previously installed these databases. Since downloading databases takes a while, we recommend using our ```--email``` option: MMT will email you with updates and when all databases have been installed.


We've included test data to check your installation:
```
mmt full -test --email youremail@
```
A full test takes about 40 minutes on 12 cores (30 mins of which is the assembly with Trinity).




