###MakeMyTranscriptome (MMT) is an automated pipeline for *de novo* transcriptomics.

For a full guide, check out the [Wiki](http://github.com/bluegenes/makeMyTranscriptome/wiki)



MMT is composed of four principle tools: 
- **assembly** - assemble input reads into a transcriptome
- **quality** - assess the quality of a transcriptome using common length-based and annotation-based metrics
- **annotation** - annotate a transcriptome using ORF-prediction and BLAST-style mappings
- **expression**  - identify differentially expressed transcripts and genes 

A fifth tool is currently being developed:
- **compare** - compare MMT assemblies constructed with different parameters

Each tool can be run individually, or all four can be executed as a single virtual tool, **full**. All tools make use of two additional convenience modules, **databases** and **tools** which download and install databases and tools, respectively. For non-linux machines, the tools module just checks available tools and provides installation instructions.


### Quick Start Guide:

MMT is a convenience tool that runs a number of bioinformatics tools in an automated fashion. In order to use it, you'll have to download these other tools, and put these into your PATH. If you're working on a Linux machine, MMT can download the tools for you; otherwise, it will assess the tools you have available and print recommended installation instructions.

First, download MMT via git (below) or click the "Download ZIP" button.
```
git clone https://github.com/bluegenes/MakeMyTranscriptome.git
```
cd into the MakeMyTranscriptome directory.
```
cd MakeMyTranscriptome
```

In order to run MMT, the general structure is:
```
mmt [TOOL_SELECTOR] [ARGUMENTS]  
``` 

For example, to run a full test of the pipeline on sample data:
```
mmt full -test
```



