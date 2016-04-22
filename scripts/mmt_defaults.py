from os.path import join, dirname, abspath
from os import makedirs


''' dir variables'''
PATH_ROOT = dirname(dirname(abspath(__file__)))
PATH_TEST = join(PATH_ROOT, 'test_data')
PATH_SCRIPTS = join(PATH_ROOT, 'scripts')
PATH_DATABASES = join(PATH_ROOT, 'databases')
PATH_DATABASE_LOGS = join(PATH_DATABASES, 'log_files')
PATH_ASSEMBLIES = join(PATH_ROOT, 'assemblies')
PATH_TOOLS = join(PATH_ROOT, 'external_tools')
PATH_UTIL = join(PATH_SCRIPTS, 'util')

''' db variables '''
PATH_UNIPROT_SPROT_DIR = join(PATH_DATABASES, 'uniprot_sprot')
PATH_UNIPROT_SPROT = join(PATH_UNIPROT_SPROT_DIR, 'uniprot_sprot.fasta')
PATH_UNIREF90_DIR = join(PATH_DATABASES, 'uniref90')
PATH_UNIREF90 = join(PATH_UNIREF90_DIR, 'uniref90.fasta')
PATH_NR_DIR = join(PATH_DATABASES, 'nr')
PATH_NR = join(PATH_NR_DIR, 'nr.fasta')
PATH_BUSCO_REFERENCE_DIR = join(PATH_DATABASES, 'busco')
PATH_BUSCO_METAZOA = join(PATH_BUSCO_REFERENCE_DIR, 'metazoa_buscos')
PATH_BUSCO_ANTHROPODA = join(PATH_BUSCO_REFERENCE_DIR, 'arthropoda_buscos')
PATH_BUSCO_VERTEBRATA = join(PATH_BUSCO_REFERENCE_DIR, 'vertebrata_buscos')
PATH_BUSCO_EUKARYOTA = join(PATH_BUSCO_REFERENCE_DIR, 'eukaryota_buscos')
PATH_BUSCO_FUNGI = join(PATH_BUSCO_REFERENCE_DIR, 'fungi_buscos')
PATH_BUSCO_BACTERIA = join(PATH_BUSCO_REFERENCE_DIR, 'bacteria_buscos')
PATH_BUSCO_PLANT = join(PATH_BUSCO_REFERENCE_DIR, 'plantae_buscos')
PATH_PFAM_DIR = join(PATH_DATABASES, 'pfam')
PATH_PFAM_DATABASE = join(PATH_PFAM_DIR, 'Pfam-A.hmm')
PATH_NOG_CATEGORIES = join(PATH_DATABASES, 'nog_categories')
PATH_GO_PATHWAY = join(PATH_DATABASES, 'go_pathway.txt')
PATH_SWISS_ENZYME = join(PATH_DATABASES, 'swiss_enzyme.list')
PATH_PFAM_ENZYME = join(PATH_DATABASES, 'pfam_enzyme.list')
PATH_ID_MAPPING_DIR = join(PATH_DATABASES, 'id_mapping')
PATH_ID_MAPPING = join(PATH_ID_MAPPING_DIR, 'idmapping.dat')
PATH_ID_MAPPING_BIOCYC = PATH_ID_MAPPING+'.biocyc'
PATH_ID_MAPPING_EGGNOG = PATH_ID_MAPPING+'eggNOG'
PATH_ID_MAPPING_KO = PATH_ID_MAPPING+'.KO'
PATH_ID_MAPPING_ORTHODB = PATH_ID_MAPPING+'orthodb'
PATH_UNIPROT_SPROT_MAP = join(PATH_DATABASES, 'uniprot_sprot.dat')
PATH_KOG_FUNCTIONAL = join(PATH_DATABASES, 'allKOG_functional_info.txt')
PATH_SLIM_GENERIC = join(PATH_DATABASES, 'goslim_generic.obo')
PATH_ID_MAPPING_SELECTED = join(PATH_DATABASES, 'idmapping_selected.tab')
PATH_DATABASE_LOG = join(PATH_DATABASES, '.database_supervisor_log ')
PATH_DB_CONFIG_FILE = join(PATH_ROOT, 'db_config.json')
PATH_SWISS_ENZYME = join(PATH_DATABASES, 'swiss_enzyme.list')
PATH_ENZYME_PATHWAY = join(PATH_DATABASES, 'enzyme_pathway.list')
PATH_ORTHOLOGY_PATHWAY = join(PATH_DATABASES, 'orthology_pathway.list')


''' url variables '''
URL_BUSCO_BACTERIA = 'http://busco.ezlab.org/files/bacteria_buscos.tar.gz'
URL_BUSCO_FUNGI = 'http://busco.ezlab.org/files/fungi_buscos.tar.gz'
URL_BUSCO_EUKARYOTA = 'http://busco.ezlab.org/files/eukaryota_buscos.tar.gz'
URL_BUSCO_VERTEBRATA = 'http://busco.ezlab.org/files/vertebrata_buscos.tar.gz'
URL_BUSCO_ANTHROPODA = 'http://busco.ezlab.org/files/arthropoda_buscos.tar.gz'
URL_BUSCO_METAZOA = 'http://busco.ezlab.org/files/metazoa_buscos.tar.gz'
URL_BUSCO_PLANT = 'http://buscos.ezlab.org/files/plant_early_release.tar.gz'
URL_UNIPROT_SPROT = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz'
URL_UNIREF90 = 'ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz'
URL_NR = 'ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz'
URL_GO_PATHWAY = 'http://rest.genome.jp/link/go/pathway'
URL_PFAM_ENZYME = 'http://rest.genome.jp/link/enzyme/pfam'
URL_ID_MAPPING = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz'
URL_UNIPROT_SPROT_MAP = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz'
URL_KOG_FUNCTIONAL =  'http://eggnogdb.embl.de/download/eggnog_4.1/data/NOG/NOG.annotations.tsv.gz'
URL_SLIM_GENERIC = 'http://www.geneontology.org/ontology/subsets/goslim_generic.obo'
URL_PFAM_DATABASE = 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam//releases/Pfam28.0/Pfam-A.hmm.gz'
URL_ID_MAPPING_SELECTED = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz'
URL_SWISS_ENZYME = 'http://rest.genome.jp/link/enzyme/swissprot'
URL_ORTHOLOGY_PATHWAY = None
URL_ENZYME_PATHWAY = None
URL_NOG_CATEGORIES = None


class Output_Path_Vars:

    def __init__(self, basename, out_dir=PATH_ASSEMBLIES):
        ''' A simple class that simplifies generating names during creation of
            tasks. The only two input are "basename", which should be the name
            of the output dir you wish to create and have paths for, along with
            an optional "out_dir", which species the directory in which the
            out_dir should be created. The paths will be assesible as member
            variables. The variables are listed below...

            self.assembly_name          : basename for output
            self.path_dir               : path to output dir
            self.path_assembly_files    : path to assembly files
            self.path_quality_files     : path to quality files
            self.path_annotation_files  : path to annotation files
            self.path_expression_files  : path to expression files
            self.path_filter_files      : path to filtered_assemblies dir
            self.path_logs              : path to log file containing dir
            self.path_assembly          : path to the fasta assembly
            self.path_gene_trans_map    : path to .gene_trans_map output
            self.path_transdecoder_dir  : path to transdecoder dir
            self.path_transrate_dir     : path to transrate dir
            self.path_pep               : path to transdecoder pep file
            self.path_annot_table       : path to annotation table
            self.path_history           : path to history.json
        '''
        self.assembly_name = basename
        self.path_dir = join(out_dir, basename)
        self.path_assembly_files = join(self.path_dir, 'assembly_files')
        self.path_quality_files = join(self.path_dir, 'quality_files')
        self.path_annotation_files = join(self.path_dir, 'annotation_files')
        self.path_expression_files = join(self.path_dir, 'expression_files')
        self.path_filter_files = join(self.path_dir, 'filtered_assemblies')
        self.path_logs = join(self.path_dir, 'log_files')
        self.path_assembly = join(self.path_dir, basename + '.fasta')
        self.path_gene_trans_map = join(self.path_assembly_files, basename + '.gene_trans_map')
        self.path_transdecoder_dir = join(self.path_annotation_files, 'transdecoder')
        self.path_transrate_dir = join(self.path_quality_files, 'transrate')
        self.path_pep = join(self.path_transdecoder_dir, basename + '.fasta.transdecoder.pep')
        self.path_annot_table = join(self.path_dir, basename + '_annotation.txt')
        self.path_history = join(self.path_logs, 'history.json')

    def build(self):
        dirs = [self.path_dir, self.path_assembly_files, self.path_quality_files,
                self.path_annotation_files, self.path_expression_files, self.path_filter_files,
                self.path_logs, self.path_transdecoder_dir, self.path_transrate_dir]
        for d in dirs:
            try:
                makedirs(d)
            except:
                pass
