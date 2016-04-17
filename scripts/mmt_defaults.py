from os.path import join, dirname, abspath


''' dir variables'''
PATH_ROOT = dirname(dirname(abspath(__file__)))
PATH_TEST = join(PATH_ROOT, 'test_data')
PATH_SCRIPTS = join(PATH_ROOT, 'scripts')
PATH_DATABASES = join(PATH_ROOT, 'databases')
PATH_DATABASE_LOGS = join(PATH_DATABASES, 'log_files')
PATH_ASSEMBLIES = join(PATH_ROOT, 'assemblies')
PATH_TOOLS = join(PATH_ROOT, 'external_tools')

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
PATH_ID_MAPPING = join(PATH_DATABASES, 'idmapping.dat')
PATH_UNIPROT_SPROT_MAP = join(PATH_DATABASES, 'uniprot_sprot.dat')
PATH_KOG_FUNCTIONAL = join(PATH_DATABASES, 'allKOG_functional_info.txt')
PATH_SLIM_GENERIC = join(PATH_DATABASES, 'goslim_generic.obo')
PATH_IDMAPPING_SELECTED = join(PATH_DATABASES, 'idmapping_selected.tab')
PATH_DATABASE_LOG = join(PATH_DATABASES, '.database_supervisor_log ')

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
URL_IDMAPPING_SELECTED = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz'


''' misc variables'''
PATH_DB_CONFIG_FILE = join(PATH_ROOT, 'db_config.json')
