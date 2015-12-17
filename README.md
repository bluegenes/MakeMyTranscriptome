###MakeMyTranscriptome (MMT) is an automated pipeline for *de novo* transcriptomics.

For a full guide, check out the [Wiki](http://github.com/bluegenes/makeMyTranscriptome/wiki)



MMT is composed of four principle tools: 
- **assembly** - assemble input reads into a transcriptome
- **quality** - assess the quality of a transcriptome using common length-based and annotation-based metrics
- **annotation** - annotate a transcriptome using ORF-prediction and BLAST-style mappings
- **expression**  - identify differentially expressed transcripts and genes. 

Each tool can be run individually, or all four can be executed as a single virtual tool, **full**.


### Quick Start Guide:

MMT is a convenience tool that runs a number of bioinformatics tools in an automated fashion. In order to use it, you'll have to download these other tools, and put these into your PATH. If you're working on an Ubuntu machine with root access, MMT can download the tools for you (coming soon). Instructions [here](https://github.com/bluegenes/makeMyTranscriptome/wiki/Install).

Once you have installed the required tools, cd into the MakeMyTranscriptome directory and run:
```
mmt [TOOL_SELECTOR] [ARGUMENTS]  
``` 
for example,
```
mmt full -test
```
will run a full test on the sample data.
