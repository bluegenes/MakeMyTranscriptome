## Input
The Assembly tool is a pipeline that generates an assembly from a set of fastq reads. To execute the assembly tool, call
```
python pipeline.py assembly [INPUT] [OPTIONS] [FLAGS]


Input - At least one of the following must be provided.
  -u UNPAIRED, --unpaired UNPAIRED
                        A comma separated list of unpaired fastq files.
  -1 FASTQ1, --fastq1 FASTQ1
                        A comma separated list of fastq files. Each file
                        should be paired with the same indexed file in fastq2.
  -2 FASTQ2, --fastq2 FASTQ2
                        A comma seperated list of fastq files. Each file
                        should be paired with the same indexed file in fastq1.

Options
  --subsample_size SUBSAMPLE_SIZE
                        If greater than this number of reads (in millions) is
                        provided, sub sample down to this number. Use 0 to
                        signal that no subsampling should be performed. The
                        default value is 50 meaning that 50 million randomly 
                        selected paired reads will be provided to the assembler.
  --subsample_seed SUBSAMPLE_SEED
                        A seed used to initialize the random number generator
                        used during random sampling. You can set this option in 
                        order to recreate an earlier generated sample. The default 
                        value is the current system time.
  --busco_ref BUSCO_REF
                        A lineage-specific profile library available with busco. 
                        You can think of it as a reference that busco will use to 
                        evaluate the quality of your assembly. Acceptable values are
                        "Metazoa", "Eukaryota", "Arthropoda", "Vertebrata", "Fungi",
                        or "Bacteria". The default value is set arbitrarily to 
                        "Metazoa". It is strongly recommended that you set this to 
                        be an approrpriatte value for the organism whose assembly 
                        you are trying to generate.
  --cpu CPU             Sets the process cap for execution. Default is 12. Use
                        0 to indicate no process cap should be used.
  --email EMAIL         Pipeline will send emails informing you of current runstate.
  --out_dir OUT_DIR     Path to the output location. Defaults to assemblies
                        directory inside pipeline
  -o OUT_NAME, --out_name OUT_NAME
                        The name of the output directory to be made in
                        out_dir. If unused, name will be inherited from input
                        file names
Flags
  -rnaspades            Use this flag to specify that assembly should be
                        performed by rnaSPAdes rather than the default
                        Trinity.
  -no_rmdup             Use this flag to disable the removing duplicates
                        portion of the pre-assembly read cleaning.
  -no_trim              Use this flag to disable all trimming portions of pre-
                        assembly read cleaning. Duplicate and low quality
                        reads will not be removed. Subsampling will still be
                        executed. 
  -cegma                Use this flag to run cegma on the created assembly. 
                        Cegma is a legacy  tool that analyzes the quality 
                        of an assembly. It is a no longer supported tool that 
                        is being replaced by Busco. Using this flag only 
                        signals the pipeline to run cegma. Busco will still 
                        be performed. 
  -no_log               Pipeline will delete log files.
  -force                Use this flag to perform a fresh run of the pipeline.
                        All steps will be executed regardless of what has
                        already been performed.
```
## Output
The output of the assembly tool is an assembly directory, identical in basic structure as every other tool. The principle output is located in the root of the assembly directory and called "myassembly.fasta". It is the complete assembly constructed given your specified parameters. What follows is a full description of the assembly directory constructed by running the assembly tool. In the following tree, "--" mean that the object is a file while ">>" means the object is a directory.

```
<out_name>
    --<out_name>.fasta 
    --run.log 
    >>annotation_files 
    >>assembly_files
        >>fastqc_pre_trimming 
        >>fastqc_post_trimming_paired
        >>fastqc_post_trimming_unpaired
        >>fastqc_trinity_input_paired 
        --prinseq_output_<K>_<input_name>
        --assembly_stats.json
        >>rna_spades_out_dir
        >>run_busco_<busco_ref>
        >>trinity
        --final_reads_1.fastq
        --final_reads_1.fastq.readcount 
        --final_reads_2.fastq
        --final_reads_2.fastq.readcount
    >>expression_files
    >>logs
```

### Description of output files and directories
**&lt;out\_name>** - This directory is the output directory generated by running OCT. It is taken directly from the "--out\_name"  options. If the option is not specified, OCT will construct a name from the input fastq file names. It is strongly recommended that you use the "--out\_name" option.

**&lt;out\_name>.fasta** - This file contains the assembly generated by running the assembly tool. It inherits its name from the name of the output directory

**run.log** - This file contains a log of all of the jobs performed by OCT, and when they were performed.

**annotation\_files** - This directory contains all files generated by running the annotation tool. It is not populated by the assembly tool. to populate this directory, execute the annotation tool.

**assembly\_files** - This directory contains all files generated by running the assembly tool. 

**fastqc\_pre\_trimming** - This directory contains the output from running fastqc on the unaltered input fastq files. For details on the contents of this directory, consult the fastqc documentation found [at this link](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

**fastqc\_post\_trimming\_paired** - This directory contains the output from running fastqc on the cleaned paired fastq files generated by running prinseq on the provided paired files. This directory will not be present if the "-no\_trim" flag is used while running the assembly tool or if no paired reads were provided. For details on the contents of this directory, consult the fastqc documentation found [at this link](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

**fastqc\_post\_trimming\_unpaired** - This directory contains the output from running fastqc on the cleaned unpaired fastq files generated by running prinseq on the provided unpaired files. This directory will not be present if the "-no\_trim" flag is used while running the assembly tool or if no unpaired reads were provided. For details on the contents of this directory, consult the fastqc documentation found [at this link](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

**fastqc\_final\_reads\_paired** - This directory contains the output from running fastqc on the subsampled paired fastq. This directory will not be generated if no paired reads were provided. For details on the contents of this directory, consult the fastqc documentation found [at this link](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

**prinseq\_output\_&lt;K>\_&lt;input\_name>** - A file of this format will be generated for every fastq file provided as input. The K value is an index specifying the order the input files were provided in starting with paired fastq files then moving to unpaired files. The input name is inhereted directly from the input files. For example, if you were to provide "-1 reads_1.fastq -2 reads_2.fastq -u sinlges.fastq" then you would see three files named prinseq\_output\_0\_1\_reads\_1.fastq, prinseq\_output\_0\_2\_reads\_2.fastq, and prinseq\_output\_1\_singles.fastq inside the assembly_files directory.

**assembly_stats.json** - This file contains a few basic stats about the assembly such as median transcript length, mean transcript length, n50, total length across all transcripts, gc content, and the number of transcripts.

**rna\_spades\_out\_dir** - This is the output directory generated by running rnaSPAdes on the finalized reads. This directory will only be present if the '-rnaspades' flag is used. For details on the contents of the directory consult the rnaSPAdes documentation found [at this link](http://spades.bioinf.spbau.ru/rnaspades0.1.1/rnaspades_manual.html).

**run\_busco\_&lt;busco\_ref>** - This directory contains the output of running Busco on the created assembly. It will inherit part of its name from the busco lineage specific profile library specified by "--busco_ref". For details on the contents of this directory, consult the busco documentation found [at this link](http://busco.ezlab.org/).

**trinity** - This is the output directory generated by running Trinity on the finalized reads. If the '-rnaspades' flag is used, then rnaSPAdes will be used in place of trinity and the rna\_spades\_out\_dir directory will be generated instead of the trinity directory. For details on the contents of the directory consult the Trinity documentation found [at this link](https://github.com/trinityrnaseq/trinityrnaseq/wiki).

**final\_reads\_1.fastq** - Contains all left (1) reads that will be passed to the appropriate assembler for assembly. This file will only be generated if paired reads were provided as input.

**final\_reads\_1.fastq.readcount** - The total number of reads in final\_reads\_1.fastq. This file will only be generated if paired reads were provided as input.

**final\_reads\_2.fastq** - Contains all right (2) reads that will be passed to the appropriate assembler for assembly. This file will only be generated if paired reads were provided as input.

**final\_reads\_2.fastq.readcount** - The total number of reads in final\_reads\_2.fastq. This file will only be generated if paired reads were provided as input.

**expression\_files** - This directory contains all files generated by running the expression tool. It is not populated by the assembly tool. To populate this directory, execute the expression tool.

**logs** -  This directory contains logs for every job performed by OCT. It is partially populated by every tool.
