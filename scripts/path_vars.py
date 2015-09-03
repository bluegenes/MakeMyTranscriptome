import os


''' name variables '''
NAME_ASSEMBLY = 'myassembly'
NAME_OUT_DIR = 'pipeline_output'

''' static path variables '''
PATH_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PATH_SCRIPTS = os.path.join(PATH_ROOT,'scripts')
PATH_DATABASES = os.path.join(PATH_ROOT,'databases')
PATH_ASSEMBLIES = os.path.join(PATH_ROOT,'assemblies')
PATH_BEDTOOLS = '/matta1/biotools/redhat/bedtools2-master/bin/bedtools'
PATH_BLASTP = ' /matta1/biotools/redhat/ncbi-blast-2.2.30+/bin/blastp'
PATH_BLASTX = '/matta1/biotools/redhat/ncbi-blast-2.2.30+/bin/blastx'
PATH_BOWTIE2 = '/matta1/biotools/redhat/bowtie2-2.1.0'
PATH_BUSCO = '/matta1/biotools/redhat/BUSCO_v1.1b/BUSCO_v1.1b.py'
PATH_BUSCO_REFERENCE = '/matta1/hitsdata/reference_files/BUSCO'
PATH_CEGMA = 'cegma'
PATH_EXPRESS = '/matta1/biotools/redhat/express-1.5.1-linux_x86_64/express'
PATH_FASTQC = '/matta1/biotools/redhat/FastQC/fastqc'
PATH_GENE_TRANS_MAP = PATH_SCRIPTS+'/get_Trinity_gene_to_trans_map.pl'
PATH_KALLISTO = '/matta1/biotools/redhat/kallisto_linux-v0.42.1/kallisto'
PATH_PFAM = 'hmmscan'
PATH_PFAM_DB = '/matta1/hitsdata/reference_files/for_trinotate/for_transDecoder/Pfam-AB.hmm.bin'
PATH_PRINSEQ = 	'/matta1/biotools/redhat/prinseq-lite-0.20.4/prinseq-lite.pl'
PATH_RNAMMER = '/matta1/biotools/redhat/rnammer-1.2/rnammer'
PATH_RNAMMER_PL = '/matta1/biotools/redhat/Trinotate_r20140708/util/rnammer_support/RnammerTranscriptome.pl'
PATH_SALMON = '/matta1/biotools/redhat/SalmonBeta-0.4.2_DebianSqueeze/bin/salmon'
PATH_SIGNALP = 'signalp'
PATH_RNASPADES = '/matta1/biotools/redhat/rnaSPAdes-0.1.0-Linux/bin/rnaspades.py'
PATH_TMHMM = 'tmhmm'
PATH_TRANSDECODER = '/matta1/biotools/redhat/trinityrnaseq_r20140717/trinity-plugins/transdecoder/TransDecoder'
PATH_TRINITY = '/matta1/biotools/redhat/trinityrnaseq_r20140717/Trinity'
PATH_SWISS_PROT = os.path.join(PATH_DATABASES,'uniprot_sprot','uniprot_sprot.fasta')
PATH_UNIREF90 = os.path.join(PATH_DATABASES,'uniref90','uniref90.fasta')
PATH_NOG_CATEGORIES = os.path.join(PATH_DATABASES,'nog_categories')

''' dynamic path variable functions'''
GEN_PATH_DIR = lambda : os.path.join(PATH_ASSEMBLIES,NAME_OUT_DIR)
GEN_PATH_ASSEMBLY_FILES = lambda : os.path.join(GEN_PATH_DIR(),'assembly_files')
GEN_PATH_ANNOTATION_FILES = lambda : os.path.join(GEN_PATH_DIR(),'annotation_files')
GEN_PATH_EXPRESSION_FILES = lambda : os.path.join(GEN_PATH_DIR(),'expression_files')
GEN_PATH_LOGS = lambda : os.path.join(GEN_PATH_DIR(),'logs')
GEN_PATH_ASSEMBLY = lambda : os.path.join(GEN_PATH_DIR(),NAME_ASSEMBLY+'.fasta')

GEN_LOGS = lambda x : (os.path.join(GEN_PATH_LOGS(),x+'.out_log'),os.path.join(GEN_PATH_LOGS(),x+'.err_log'))


if(not os.path.isdir(PATH_DATABASES)):
	os.makedirs(PATH_DATABASES)
if(not os.path.isdir(PATH_ASSEMBLIES)):
	os.makedirs(PATH_ASSEMBLIES)

