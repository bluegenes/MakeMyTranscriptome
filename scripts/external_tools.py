# author: bluegenes

from tools_class import ext_tool as tc 
from os.path import join,dirname,abspath
import platform

PATH_ROOT = dirname(dirname(abspath(__file__)))
#PATH_SCRIPTS = join(PATH_ROOT, 'scripts')
#PATH_DATABASES = join(PATH_ROOT, 'databases')
#PATH_ASSEMBLIES = join(PATH_ROOT, 'assemblies')
PATH_TOOLS = join(PATH_ROOT, 'external_tools')

TOOL_LIST = []
#TOOL_LIST = [trinity_tool, trimmomatic_tool, prinseq_tool, transdecoder_tool, transrate_tool,busco_tool, hmmer_tool, diamond_tool, salmon_tool, busco_plant_tool, fastqc_tool]

### Trinity ###
trinity_source_url = 'https://github.com/trinityrnaseq/trinityrnaseq/archive/v2.1.1.tar.gz'
trinity_source_target = join(PATH_TOOLS, 'trinityrnaseq-2.1.1')
trinity_exe = ['Trinity', 'util/support_scripts/get_Trinity_gene_to_trans_map.pl']
trinity_instructions = "trinity instructions here"
trinity_cmd = 'make'
#trinity_dependencies= ['bowtie-1'] ### not doing anything with this yet... need to add to tools class?

trinity_tool = tc('trinity', trinity_source_url, trinity_source_target, trinity_exe, trinity_instructions)
trinity_tool.set_install(trinity_cmd)
TOOL_LIST.append(trinity_tool)

### Trimmomatic ###
trimmomatic_binary_url = 'http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.35.zip'
trimmomatic_binary_target = join(PATH_TOOLS, 'Trimmomatic-0.35')
trimmomatic_exe = ['trimmomatic-0.35.jar']
trimmomatic_instructions = "trimmomatic instructions here"
trimmomatic_urltype = 'zip'

trimmomatic_tool = tc('trimmomatic', trimmomatic_binary_url, trimmomatic_binary_target, trimmomatic_exe, trimmomatic_instructions, urltype=trimmomatic_urltype)
TOOL_LIST.append(trimmomatic_tool)

### Prinseq ###
prinseq_url = 'http://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz'
prinseq_target = join(PATH_TOOLS, 'prinseq-lite-0.20.4')
prinseq_exe = ['prinseq-lite.pl']
prinseq_instructions = "prinseq instructions here"

prinseq_tool = tc('prinseq', prinseq_url, prinseq_target, prinseq_exe, prinseq_instructions) 
TOOL_LIST.append(prinseq_tool)

### Transdecoder ###
transdecoder_url = 'https://github.com/TransDecoder/TransDecoder/archive/2.0.1.tar.gz'
transdecoder_target = join(PATH_TOOLS, 'TransDecoder-2.0.1')
transdecoder_exe = ['TransDecoder.Predict','TransDecoder.LongOrfs']
transdecoder_instructions = "transdecoder instructions here"
transdecoder_cmd= 'make'

transdecoder_tool = tc('transdecoder', transdecoder_url, transdecoder_target, transdecoder_exe, transdecoder_instructions)
transdecoder_tool.set_install(transdecoder_cmd)
TOOL_LIST.append(transdecoder_tool)

### Transrate ###
transrate_linux_url = 'https://bintray.com/artifact/download/blahah/generic/transrate-1.0.1-linux-x86_64.tar.gz'
transrate_linux_target = join(PATH_TOOLS, 'transrate-1.0.1-linux-x86_64')
transrate_exe = ['transrate']
transrate_instructions = "transrate instructions here"
transrate_cmd = 'transrate --install-deps=ref'

transrate_tool = tc('transrate', transrate_linux_url, transrate_linux_target, transrate_exe, transrate_instructions)
transrate_tool.set_install(transrate_cmd)
TOOL_LIST.append(transrate_tool)

### BUSCO ###
busco_url = 'http://busco.ezlab.org/files/BUSCO_v1.1b1.tar.gz'
busco_target = join(PATH_TOOLS, 'BUSCO_v1.1b1')
busco_exe = ['BUSCO_v1.1b1.py']
busco_instructions = "busco instructions here"

busco_tool = tc('busco', busco_url, busco_target, busco_exe, busco_instructions)
TOOL_LIST.append(busco_tool)

#plant_busco
busco_plant_url = 'http://buscos.ezlab.org/files/plant_early_release.tar.gz'
busco_plant_target = join(PATH_TOOLS, 'plant_early_release')
busco_plant_exe = ['BUSCO_plants.py']
busco_plant_instructions = 'busco plant instructions here'

busco_plant_tool = tc('busco_plant', busco_plant_url, busco_plant_target, busco_plant_exe, busco_plant_instructions)
TOOL_LIST.append(busco_plant_tool)

### HMMER ###
hmmer_linux_url = 'http://selab.janelia.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz'
hmmer_linux_target = join(PATH_TOOLS, 'hmmer-3.1b2-linux-intel-x86_64')
hmmer_exe = ['binaries/hmmscan','binaries/hmmpress']
hmmer_instructions = "hmmer instructions here"

hmmer_tool = tc('hmmer',hmmer_linux_url, hmmer_linux_target, hmmer_exe, hmmer_instructions)
TOOL_LIST.append(hmmer_tool)

### Diamond ###
diamond_linux_url = 'http://github.com/bbuchfink/diamond/releases/download/v0.7.10/diamond-linux64.tar.gz'
diamond_target = join(PATH_TOOLS, 'diamond')
diamond_exe = ['diamond']
diamond_instructions = "diamond instructions here"

diamond_tool = tc('diamond', diamond_linux_url, diamond_target, diamond_exe, diamond_instructions)
diamond_tool.change_exe_fullpath(PATH_TOOLS) # bc reg is PATH_TOOLS/target/exe
TOOL_LIST.append(diamond_tool)

### Salmon ###
salmon_linux_url = 'https://github.com/COMBINE-lab/salmon/releases/download/v0.6.0/SalmonBeta-0.6.0_DebianSqueeze.tar.gz'
salmon_linux_target = join(PATH_TOOLS, 'SalmonBeta-0.6.1_DebianSqueeze')
salmon_exe = ['bin/salmon']
salmon_instructions = "salmon instructions here"

salmon_tool = tc('salmon', salmon_linux_url, salmon_linux_target, salmon_exe, salmon_instructions)
TOOL_LIST.append(salmon_tool)

### FastQC ###
fastqc_linux_url = 'http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.4.zip'
fastqc_linux_target = 'fastqc_v0.11.4'
fastqc_exe = ['fastqc']
fastqc_instructions = 'fastqc instructions here'
fastqc_folder_name = 'FastQC' #unzips into FastQC dir
fastqc_urltype = 'zip'

fastqc_tool = tc('fastqc', fastqc_linux_url, fastqc_linux_target, fastqc_exe, fastqc_instructions, urltype='zip', folder_name=fastqc_folder_name)
fastqc_tool.change_exe_fullpath(join(PATH_TOOLS,'FastQC')) #unzips into FastQC dir
fastqc_cmd = 'chmod 755 '+ fastqc_tool.full_exe[0]
fastqc_tool.set_install(fastqc_cmd)
TOOL_LIST.append(fastqc_tool)

cegma_url =  'http://korflab.ucdavis.edu/datasets/cegma/CEGMA_v2.5.tar.gz'
cegma_target = 'CEGMA_v2.5'
cegma_exe = ['cegma']
cegma_instructions = 'MMT cannot install cegma. Please see http://korflab.ucdavis.edu/datasets/cegma/#SCT3 for installation instructions.'

cegma_tool = tc('cegma', cegma_url, cegma_target, cegma_exe, cegma_instructions)
cegma_tool.change_exe_fullpath('') # they need to put cegma into their $path
TOOL_LIST.append(cegma_tool)


bedtools_url = 'https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz'
bedtools_target = 'bedtools-2.25.0'
bedtools_cmd = 'make'
bedtools_instructions = 'see installation instructions here: http://bedtools.readthedocs.org/en/latest/content/installation.html'
bedtools_exe = ['intersectBed']
bedtools_folder_name = 'bedtools2' #unpacks to 'bedtools2'

bedtools_tool = tc('bedtools',bedtools_url, bedtools_target, bedtools_exe, bedtools_instructions, folder_name=bedtools_folder_name)
bedtools_tool.change_exe_fullpath(join(PATH_TOOLS,'bedtools2/bin'))
TOOL_LIST.append(bedtools_tool)



# set up tools dictionary for use in all mmt
#if(platform.system().lower() == 'linux')
TOOLS_DICT = {tool.name: tool for tool in TOOL_LIST}


""" TOOLS/INFO TO ADD !!!!!!
PATH_TRANSDECODER = 'TransDecoder'
PATH_TRIMMOMATIC_ADAPTERS_SINGLE = join(PATH_TOOLS,'Trimmomatic-0.35/adapters/TruSeq3-SE.fa')
PATH_TRIMMOMATIC_ADAPTERS_PAIRED = join(PATH_TOOLS,'Trimmomatic-0.35/adapters/TruSeq3-PE.fa')

''' set these up as executables that we just look for in the $path '''
PATH_BEDTOOLS = 'bedtools'
PATH_BLASTP = 'blastp'
PATH_BLASTX = 'blastx'
PATH_BOWTIE2 = ''
PATH_CEGMA = 'cegma'
PATH_EXPRESS = 'express'
PATH_KALLISTO = 'kallisto'

'''optional tools -- set up the same way as above? just don't check for them with mmt setup; only install instructions in tools'''
PATH_RNAMMER = '/matta1/biotools/redhat/rnammer-1.2/rnammer'
PATH_RNAMMER_PL = 'RnammerTranscriptome.pl'
PATH_SIGNALP = 'signalp'
PATH_RNASPADES = 'rnaspades.py'
PATH_TMHMM = 'tmhmm'


"""
