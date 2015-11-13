###############################################################################
###############################################################################
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
###############################################################################
###############################################################################


def conversion(conversion_file):
	"""
	   Description:	Generic dictionary creator for database conversions
	   Input:		Conversion File
	   Return:		Dictionary
	"""
    conv_dt = {}
    conversion = open(conversion_file)
    for line in conversion:
        line = line.strip().split()
        first = line[0].split(":")[1]
        second = line[1].split(":")[1]
        conv_dt[first] = second
    return conv_dt


def conversion_nogf(conversion_file):
	"""
	   Description:	eggNog -> eggNog function dictionary creator for 
	   				database conversions
	   Input:		Conversion File
	   Return:		Dictionary
	"""
    conv_dt = {}
    conversion = open(conversion_file)
    for line in conversion:
        line = line.strip().split("\t")
        first = line[0][3:]
        second = line[1]
        conv_dt[first] = second
    return conv_dt


def conversion_contig_blastnr(conversion_file):
	"""
	   Description:	Contig -> BlastNR dictionary creator for database 
	   				conversions
	   Input:		Conversion File
	   Return:		Dictionary
	"""
    conv_dt = {}
    conversion = open(conversion_file)
    for line in conversion:
        line = line.strip().split("\t")
        first = line[0]
        second = line[1]
        third = line[2]
        fourth = line[3]
        conv_dt[first] = (second, third, fourth)
    return conv_dt


"""Question for Tessa: Need to account for duplicate contigs?"""
def conversion_contig_closest(conversion_file):
	"""
	   Description:	Contig -> Closest Hit dictionary creator for 
	   				database conversions
	   Input:		Conversion File
	   Return:		Dictionary
	"""
    conv_dt = {}
    conversion = open(conversion_file)
    prev = ""
    for line in conversion:
        line = line.strip().split("\t")
        if prev == line[0]:
        	continue
        first = line[0]
        second = line[12]
        third = line[14]
        conv_dt[first] = (second, third)
        prev = first
    return conv_dt


def conversion_ez_path(conversion_file):
	"""
	   Description:	Enzyme -> Kegg Pathway dictionary creator for
	   				database conversions
	   Input:		Conversion File
	   Return:		Dictionary
	"""
    conv_dt = {}
    conversion = open(conversion_file)
    for line in conversion:
        line = line.strip().split("\t")
        if line[2] == 'original':
            first = line[0].split(":")[1]
            second = line[1].split(":")[1]
            conv_dt[first] = second.replace(",", ";")
    return conv_dt


def conversion_goslim(conversion_file):
	"""
	   Description:	GO -> GO Slim dictionary creator for database 
	   				conversions
	   Input:		Conversion File
	   Return:		Dictionary
	"""
    conv_dt = {}
    conversion = open(conversion_file)
    for line in conversion:
        if line[0:7] == "id: GO:":
            conv_dt[line[4:].strip()] = line[4:].strip()
    return conv_dt


def conversion_idmap(conversion_file):
	"""
	   Description:	Generic dictionary creator for IDMAP database 
	   				conversions
	   Input:		Conversion File
	   Return:		Dictionary
	"""
    conv_dt = {}
    conversion = open(conversion_file)
    for line in conversion:
        line = line.strip("\n").split("\t")
        conv_dt[line[0]] = line[2].replace(",", ";")
    return conv_dt


def conversion_go_entrez(conversion_file):
	"""
	   Description:	SwissProt -> GO & Entrez dictionary creator for
	   				database conversions
	   Input:		Conversion File
	   Return:		Dictionary
	"""
    conv_dt = {}
    conversion = open(conversion_file)
    for line in conversion:
        line = line.split("\t")
        temp = line[6].replace(" ", "")
        conv_dt[line[0]] = (temp, line[2])
    return conv_dt

