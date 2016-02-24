# author: bluegenes

from tools_class import ext_tool as tc 
from functions_general import PATH_TOOLS
from os.path import join,dirname,abspath
import platform

#samtools --> need to either add here, or make sure user has downloaded somehow..
TOOL_LIST = []

### Trinity ###
trinity_source_url = 'https://github.com/trinityrnaseq/trinityrnaseq/archive/v2.1.1.tar.gz'
trinity_source_target = join(PATH_TOOLS, 'trinityrnaseq-2.1.1')
trinity_exe = ['Trinity', 'util/support_scripts/get_Trinity_gene_to_trans_map.pl']
trinity_instructions = '\n\t After downloading the binary: \n\n\t\t tar zxf trinityrnaseq-2.1.1.tar.gz \n\t\t ' \
                        'cd trinityrnaseq-1.1.1 \n\t\t make \n ' \
                        '\n\tThen: soft link /path/to/trinityrnaseq-2.1.1 into the MMT "external_tools" folder\n' \
                        '\t  Or: 1. add /path/to/trinityrnaseq-2.1.1 to your $path variable\n' \
                        '\t      2. add /path/to/trinityrnaseq-2.1.1/util/support_scripts/ to your $path variable\n'
trinity_cmd = 'make'
trinity_website = 'https://github.com/trinityrnaseq/trinityrnaseq/wiki'
#trinity_dependencies= ['bowtie-1'] ### not doing anything with this yet... need to add to tools class?

trinity_tool = tc('trinity', trinity_source_url, trinity_source_target, trinity_exe, trinity_instructions, web=trinity_website)
trinity_tool.set_install(trinity_cmd)
TOOL_LIST.append(trinity_tool)
#rcorrector
rcorrector_url = 'https://github.com/mourisl/Rcorrector/archive/master.zip'
rcorrector_urltype = 'zip'
rcorrector_target = join(PATH_TOOLS, 'master')
rcorrector_folder_name = 'Rcorrector-master'
rcorrector_exe = ['rcorrector'] 
rcorrector_instructions = '\n\t After downloading: \n\n\t\t unzip master.zip \n\t\t ' \
                        'cd Rcorrector-master \n\t\t make \n ' \
                        '\n\tThen: soft link /path/to/Rcorrector-master into the MMT "external_tools" folder\n' \
                        '\t  Or:  add /path/to/Rcorrector-master to your $path variable\n'
rcorrector_cmd = 'make'
rcorrector_website = 'https://github.com/mourisl/Rcorrector'

rcorrector_tool = tc('rcorrector', rcorrector_url, rcorrector_target, rcorrector_exe, rcorrector_instructions, web=rcorrector_website, urltype=rcorrector_urltype, folder_name=rcorrector_folder_name)
rcorrector_tool.change_exe_fullpath(join(PATH_TOOLS,rcorrector_folder_name))
rcorrector_tool.set_install(rcorrector_cmd)
TOOL_LIST.append(rcorrector_tool)

### Trimmomatic ###
trimmomatic_binary_url = 'http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.35.zip'
trimmomatic_binary_target = join(PATH_TOOLS, 'Trimmomatic-0.35')
trimmomatic_exe = ['trimmomatic-0.35.jar', 'adapters/TruSeq3-PE.fa','adapters/TruSeq3-SE.fa']
trimmomatic_instructions = '\n\t After downloading the binary: \n\n\t\t unzip Trimmomatic-0.35.zip \n ' \
                        '\n\tThen: soft link /path/to/Trimmomatic-0.35 into the MMT "external_tools" folder\n' \
                        '\t  Or: add /path/to/Trimmomatic-0.35 to your $path variable\n'

trimmomatic_urltype = 'zip'
trimmomatic_website = 'http://www.usadellab.org/cms/index.php?page=trimmomatic'

trimmomatic_tool = tc('trimmomatic', trimmomatic_binary_url, trimmomatic_binary_target, trimmomatic_exe, trimmomatic_instructions, urltype=trimmomatic_urltype, web=trimmomatic_website )
TOOL_LIST.append(trimmomatic_tool)

### Prinseq ###
prinseq_url = 'http://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz'
prinseq_target = join(PATH_TOOLS, 'prinseq-lite-0.20.4')
prinseq_exe = ['prinseq-lite.pl']
prinseq_instructions = '\n\t After downloading the binary: \n\n\t\t tar zxf prinseq-lite-0.20.4.tar.gz \n\t\t ' \
                        '\n\tThen: soft link /path/to/prinseq-lite-0.20.4 into the MMT "external_tools" folder\n' \
                        '\t  Or: add /path/to/prinseq-lite-0.20.4 to your $path variable\n' 
prinseq_website = 'http://prinseq.sourceforge.net'

prinseq_tool = tc('prinseq', prinseq_url, prinseq_target, prinseq_exe, prinseq_instructions, web = prinseq_website) 
TOOL_LIST.append(prinseq_tool)

### Transdecoder ###
PATH_TRANSDECODER = 'TransDecoder' #which exe to use!???

transdecoder_url = 'https://github.com/TransDecoder/TransDecoder/archive/2.0.1.tar.gz'
transdecoder_target = join(PATH_TOOLS, 'TransDecoder-2.0.1')
transdecoder_exe = ['TransDecoder.LongOrfs', 'TransDecoder.Predict']
transdecoder_instructions ='\n\t After downloading the binary: \n\n\t\t tar zxf 2.0.1.tar.gz \n\t\t ' \
                        'cd TransDecoder-2.0.1 \n\t\t make \n ' \
                        '\n\tThen: soft link /path/to/TransDecoder-2.0.1 into the MMT "external_tools" folder\n' \
                        '\t  Or: add /path/to/TransDecoder-2.0.1 to your $path variable\n'  
transdecoder_cmd= 'make'
transdecoder_website = 'http://transdecoder.github.io'

transdecoder_tool = tc('transdecoder', transdecoder_url, transdecoder_target, transdecoder_exe, transdecoder_instructions, web=transdecoder_website)
transdecoder_tool.set_install(transdecoder_cmd)
TOOL_LIST.append(transdecoder_tool)

### Transrate ###
transrate_linux_url = 'https://bintray.com/artifact/download/blahah/generic/transrate-1.0.1-linux-x86_64.tar.gz'
transrate_linux_target = join(PATH_TOOLS, 'transrate-1.0.1-linux-x86_64')
transrate_exe = ['transrate']
transrate_instructions = '\n\t After downloading the binary: \n\n\t\t tar zxf transrate-1.0.1-linux-x86_64.tar.gz \n\t\t ' \
                        'cd transrate-1.0.1-linux-x86_64 \n\t\t transrate --install-deps=ref \n ' \
                        '\n\tThen: soft link /path/to/transrate-1.0.1-linux-x86_64 into the MMT "external_tools" folder\n' \
                        '\t  Or: add /path/to/transrate-1.0.1-linux-x86_64 to your $path variable\n' 
transrate_cmd = 'transrate --install-deps=ref'
transrate_website = 'http://hibberdlab.com/transrate/'

transrate_tool = tc('transrate', transrate_linux_url, transrate_linux_target, transrate_exe, transrate_instructions, web=transrate_website)
transrate_tool.set_install(transrate_cmd)
TOOL_LIST.append(transrate_tool)

### BUSCO ###
busco_url = 'http://busco.ezlab.org/files/BUSCO_v1.1b1.tar.gz'
busco_target = join(PATH_TOOLS, 'BUSCO_v1.1b1')
busco_exe = ['BUSCO_v1.1b1.py']
busco_instructions = '\n\t After downloading the binary: \n\n\t\t tar zxf BUSCO_v1.1b1.tar.gz \n\n ' \
                        '\n\t\t cd BUSCO_v1.1b1\n' \
                        '\n\t\t chmod u+x BUSCO_v1.1b1.py\n' \
                        '\tThen: soft link /path/to/BUSCO_v1.1b1 into the MMT "external_tools" folder\n' \
                        '\t  Or: add /path/to/BUSCO_v1.1b1 to your $path variable\n'
busco_website='http://busco.ezlab.org'

busco_tool = tc('busco', busco_url, busco_target, busco_exe, busco_instructions, web=busco_website)
TOOL_LIST.append(busco_tool)

#plant_busco
busco_plant_url = 'http://buscos.ezlab.org/files/plant_early_release.tar.gz'
busco_plant_target = join(PATH_TOOLS, 'plant_early_release')
busco_plant_exe = ['BUSCO_plants.py']
busco_plant_instructions = '\n\t After downloading the binary: \n\n\t\t tar zxf plant_early_release.tar.gz \n\n ' \
                        '\n\t\t cd plant_early_release\n' \
                        '\n\t\t chmod u+x BUSCO_plants.py\n' \
                        '\tThen: soft link /path/to/plant_early_release into the MMT "external_tools" folder\n' \
                        '\t  Or: add /path/to/plant_early_release to your $path variable\n'

busco_plant_tool = tc('busco_plant', busco_plant_url, busco_plant_target, busco_plant_exe, busco_plant_instructions, web=busco_website)
TOOL_LIST.append(busco_plant_tool)

### HMMER ###
hmmer_linux_url = 'http://selab.janelia.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz'
hmmer_linux_target = join(PATH_TOOLS, 'hmmer-3.1b2-linux-intel-x86_64')
hmmer_exe = ['binaries/hmmscan','binaries/hmmpress']
hmmer_instructions = '\n\t After downloading the binary: \n\n\t\t tar zxf hmmer-3.1b2-linux-intel-x86_64.tar.gz \n\n ' \
                        '\tThen: soft link /path/to/hmmer-3.1b2-linux-intel-x86_64 into the MMT "external_tools" folder\n' \
                        '\t  Or: add /path/to/hmmer-3.1b2-linux-intel-x86_64/binaries to your $path variable\n'
hmmer_website = 'http://hmmer.org'

hmmer_tool = tc('hmmer',hmmer_linux_url, hmmer_linux_target, hmmer_exe, hmmer_instructions, web=hmmer_website)
TOOL_LIST.append(hmmer_tool)

### Diamond ###
diamond_linux_url = 'http://github.com/bbuchfink/diamond/releases/download/v0.7.10/diamond-linux64.tar.gz'
diamond_target = join(PATH_TOOLS, 'diamond')
diamond_exe = ['diamond']
diamond_instructions = '\n\t After downloading the binary: \n\n\t\t tar zxf diamond-linux64.tar.gz \n\n ' \
                        '\n\tThen: soft link path/to/diamond into the MMT "external_tools" folder\n' \
                        '\t  Or: add /path/to/diamond to your $path variable\n'
 
diamond_website = 'https://github.com/bbuchfink/diamond'

diamond_tool = tc('diamond', diamond_linux_url, diamond_target, diamond_exe, diamond_instructions,web=diamond_website)
diamond_tool.change_exe_fullpath(PATH_TOOLS) # bc reg is PATH_TOOLS/target/exe
TOOL_LIST.append(diamond_tool)

### Salmon ###
salmon_linux_url = 'https://github.com/COMBINE-lab/salmon/releases/download/v0.6.0/SalmonBeta-0.6.0_DebianSqueeze.tar.gz'
salmon_linux_target = join(PATH_TOOLS, 'SalmonBeta-0.6.1_DebianSqueeze')
salmon_exe = ['bin/salmon']
salmon_instructions = '\n\t After downloading the binary: \n\n\t\t tar zxf SalmonBeta-0.6.0_DebianSqueeze.tar.gz \n\n ' \
                         '\tThen: soft link SalmonBeta-0.6.1_DebianSqueeze into the MMT "external_tools" folder\n' \
                        '\t  Or: 1. add /path/to/SalmonBeta-0.6.1_DebianSqueeze/bin to your $path variable\n' \
                        '\t      2. add /path/to/SalmonBeta-0.6.1_DebianSqueeze/lib to your shared libraries variable\n'
 
salmon_website = 'http://combine-lab.github.io/salmon'

salmon_tool = tc('salmon', salmon_linux_url, salmon_linux_target, salmon_exe, salmon_instructions, web=salmon_website)
TOOL_LIST.append(salmon_tool)

### FastQC ###
fastqc_linux_url = 'http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.4.zip'
fastqc_linux_target = join(PATH_TOOLS,'fastqc_v0.11.4')
fastqc_exe = ['fastqc']
fastqc_instructions = '\n\t After downloading the binary: \n\n\t\t unzip fastqc_v0.11.4.zip \n\n ' \
                       '\n\tThen: soft link /path/to/fastqc_v0.11.4 into the MMT "external_tools" folder\n' \
                       '\t  Or: add /path/to/fastqc_v0.11.4 to your $path variable \n'
fastqc_folder_name = 'FastQC' #unzips into FastQC dir
fastqc_urltype = 'zip'
fastqc_website = 'http://www.bioinformatics.babraham.ac.uk/projects/fastqc'

fastqc_tool = tc('fastqc', fastqc_linux_url, fastqc_linux_target, fastqc_exe, fastqc_instructions, urltype='zip', folder_name=fastqc_folder_name, web=fastqc_website)
fastqc_tool.change_exe_fullpath(join(PATH_TOOLS,'FastQC')) #unzips into FastQC dir
fastqc_cmd = 'chmod 755 '+ fastqc_tool.full_exe[0]
fastqc_tool.set_install(fastqc_cmd)
TOOL_LIST.append(fastqc_tool)

bedtools_url = 'https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz'
bedtools_target = join(PATH_TOOLS,'bedtools-2.25.0')
bedtools_cmd = 'make'
bedtools_instructions = '\n\t After downloading the binary: \n\n\t\t tar zxf bedtools-2.25.0 \n ' \
                        '\t\t cd bedtools2 \n' \
                        '\t\t make \n\n' \
                        '\tThen: soft link /path/to/bedtools2 into the MMT "external_tools" folder\n' \
                        '\t  Or: add /path/to/bedtools2/bin to your $path variable \n\n' \
                        '\tFor more help, see: http://bedtools.readthedocs.org/en/latest/content/installation.html\n'
bedtools_exe = ['bin/intersectBed']
bedtools_folder_name = 'bedtools2' #unpacks to 'bedtools2'
bedtools_website = 'https://github.com/arq5x/bedtools2'

bedtools_tool = tc('bedtools',bedtools_url, bedtools_target, bedtools_exe, bedtools_instructions, folder_name=bedtools_folder_name, web=bedtools_website)
#bedtools_tool.change_exe_fullpath(join(PATH_TOOLS,folder_name))
TOOL_LIST.append(bedtools_tool)

#rnaspades
rnaspades_url = 'http://spades.bioinf.spbau.ru/rnaspades0.1.1/rnaSPAdes-0.1.1-Linux.tar.gz'
rnaspades_target = join(PATH_TOOLS, 'rnaSPAdes-0.1.1-Linux')
rnaspades_exe = ['bin/rnaspades.py']
rnaspades_instructions = 'rnaspades is currently not fully supported in MMT. Installation instructions will go here once it is fully supported.'

rnaspades_tool = tc('rnaspades', rnaspades_url, rnaspades_target, rnaspades_exe, rnaspades_instructions)
TOOL_LIST.append(rnaspades_tool)

blastplus_url = 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.3.0+-x64-linux.tar.gz'
blastplus_target = 'ncbi-blast-2.3.0+-x64-linux'
blastplus_folder_name = 'ncbi-blast-2.3.0+'
blastplus_instructions = '\n\t After downloading the binary: \n\n\t\t tar zxf ncbi-blast-2.3.0+-x64-linux.tar.gz \n ' \
                        '\n\t cd ncbi-blast-2.3.0+ \n' \
                        '\n\tThen: soft link /path/to/ncbi-blast-2.3.0+ into the MMT "external_tools" folder\n' \
                        '\t  Or: add /path/to/ncbi-blast-2.3.0+/bin to your $path variable \n\n'
blastplus_exe = ['bin/makeblastdb','bin/blastx', 'bin/blastp']
blastplus_website = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download'

blastplus_tool = tc('blast', blastplus_url, blastplus_target,blastplus_exe,blastplus_instructions,folder_name=blastplus_folder_name, web=blastplus_website) 
#blastplus_tool.change_exe_fullpath(join(PATH_TOOLS,folder_name))
TOOL_LIST.append(blastplus_tool)

bowtie2_url = 'http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.6/bowtie2-2.2.6-linux-x86_64.zip' 
bowtie2_urltype = 'zip'
bowtie2_target = join(PATH_TOOLS,'bowtie2-2.2.6-linux-x86_64')
bowtie2_instructions = '\n\t After downloading the binary: \n\n\t\t unzip bowtie2-2.2.6-linux-x86_64.zip \n ' \
                        '\n\t cd bowtie2-2.2.6 \n' \
                        '\n\tThen: soft link /path/to/bowtie2-2.2.6 into the MMT "external_tools" folder\n' \
                        '\t  Or: add /path/to/bowtie2-2.2.6 to your $path variable \n\n'
bowtie2_folder_name = 'bowtie2-2.2.6'
bowtie2_exe = ['bowtie2-build','bowtie2']
bowtie2_website = 'http://bowtie-bio.sourceforge.net/bowtie2/index.shtml'

bowtie2_tool = tc('bowtie2', bowtie2_url, bowtie2_target, bowtie2_exe, bowtie2_instructions, urltype=bowtie2_urltype, folder_name=bowtie2_folder_name, web=bowtie2_website)
TOOL_LIST.append(bowtie2_tool)

express_url = 'http://bio.math.berkeley.edu/eXpress/downloads/express-1.5.1/express-1.5.1-linux_x86_64.tgz'
express_urltype = '.tgz'
express_target = join(PATH_TOOLS, 'express-1.5.1-linux_x86_64.tgz')
express_instructions = '\n\t After downloading the binary: \n\n\t\t tar zxf express-1.5.1-linux_x86_64.tgz \n ' \
                        '\n\tThen: soft link /path/to/express-1.5.1-linux_x86_64 into the MMT "external_tools" folder\n' \
                        '\t  Or: add /path/to/express-1.5.1-linux_x86_64 to your $path variable \n\n'
express_exe = ['express']
express_website = 'http://bio.math.berkeley.edu/eXpress/index.html'

express_tool = tc('express',express_url, express_target, express_exe, express_instructions, urltype=express_urltype, web=express_website)
TOOL_LIST.append(express_tool)

##### tools we can't install ####
cegma_url =  'http://korflab.ucdavis.edu/datasets/cegma/CEGMA_v2.5.tar.gz'
cegma_target = join(PATH_TOOLS,'CEGMA_v2.5')
cegma_exe = ['cegma']
cegma_instructions = '\n\t MMT cannot install cegma. Please see http://korflab.ucdavis.edu/datasets/cegma/#SCT3 for installation instructions.\n'
cegma_website = 'http://korflab.ucdavis.edu/datasets/cegma/'
cegma_tool = tc('cegma', cegma_url, cegma_target, cegma_exe, cegma_instructions, install=False, web=cegma_website)
cegma_tool.change_exe_fullpath('') # they need to put cegma into their $path
TOOL_LIST.append(cegma_tool)

kallisto_url = 'https://github.com/pachterlab/kallisto/releases/download/v0.42.4/kallisto_linux-v0.42.4.tar.gz'
kallisto_target = join(PATH_TOOLS, 'kallisto_linux-v0.42.4')
kallisto_instructions = "\n\tAs kallisto is distributed under a non-commercial license, MMT cannot download kallisto for you. Please see https://pachterlab.github.io/kallisto/about.html for information about kallisto. To use, install kallisto yourself and place tool in your $path variable\n"
kallisto_exe = ['kallisto']
kallisto_website = 'https://pachterlab.github.io/kallisto/'
kallisto_tool = tc("kallisto", kallisto_url, kallisto_target, kallisto_exe, kallisto_instructions, install=False, web=kallisto_website)
kallisto_tool.change_exe_fullpath('') # look for exe in $path
TOOL_LIST.append(kallisto_tool)

#PATH_RNAMMER = '/matta1/biotools/redhat/rnammer-1.2/rnammer'
rnammer_url = ''
rnammer_target = '' #join(PATH_TOOLS,'rnammer-1.2')
rnammer_exe = ['RnammerTranscriptome.pl']
rnammer_instructions = '\n\tRNAMMER is freely available for academic use only. See http://www.cbs.dtu.dk/services/RNAmmer/ for download and installation instructions. RNAMMER is currently supported as an optional tool, but this support may be removed at any time in favor of openly licensed tools.\n'
rnammer_website = 'http://www.cbs.dtu.dk/cgi-bin/sw_request?rnammer'
rnammer_tool = tc('rnammer', rnammer_url, rnammer_target, rnammer_exe, rnammer_instructions, install=False, web=rnammer_website)
#rnammer_tool.change_exe_fullpath('') # look for exe in $path
TOOL_LIST.append(rnammer_tool)

tmhmm_url = ''
tmhmm_target = ''
tmhmm_exe = ['tmhmm']
tmhmm_instructions = '\n\tTMHMM is freely available for academic use only. See http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm for download and installation instructions. TMHMM is currently supported as an optional tool, but this support may be removed at any time in favor of openly licensed tools.\n'
tmhmm_website = 'http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm'
tmhmm_tool = tc('tmhmm', tmhmm_url, tmhmm_target, tmhmm_exe, tmhmm_instructions, install=False, web=tmhmm_website)
tmhmm_tool.change_exe_fullpath('') # look for exe in $path
TOOL_LIST.append(tmhmm_tool)

signalp_url = ''
signalp_target = ''
signalp_exe = ['signalp']
signalp_instructions = '\n\tSignalp is freely available for academic use only. See http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp for download and installation instructions. SignalP is currently supported as an optional tool, but this support may be removed at any time in favor of openly licensed tools.\n'
signalp_website= 'http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp'
signalp_tool = tc('signalp', signalp_url, signalp_target, signalp_exe, signalp_instructions, install=False, web=signalp_website)
signalp_tool.change_exe_fullpath('') # look for exe in $path
TOOL_LIST.append(signalp_tool)

TOOLS_DICT = {tool.name: tool for tool in TOOL_LIST}
