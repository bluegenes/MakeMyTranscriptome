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

import re

def conversion(conversion_file):
    """
       Description:    Generic dictionary creator for database conversions
       Input:        Conversion File
       Return:        Dictionary
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
       Description:    eggNog -> eggNog function dictionary creator for 
                       database conversions
       Input:        Conversion File
       Return:        Dictionary
    """
    conv_dt = {}
    conversion = open(conversion_file)
    for line in conversion:
        line = line.strip().split("\t")
        first = line[0][3:]
        second = line[1]
        conv_dt[first] = second
    return conv_dt

def conversion_contig_blastnr(nr_file, lengths):
    """
       Description: Generates NR BlastX Best Words conversion file and
                    the BlastX Top Hit conversion file
       Input:       BlastX NR File, Contig Length
       Return:      Dictionary
    """
    if lengths == 0:
        return {}
    file = open(nr_file)
    prev = ''
    conv_dt = {}
    contig_blasts = []
    for line in file:
        row = line.strip('\n').split('\t')
        if prev == row[0]:
            contig_blasts.append(line)
        else:
            if contig_blasts != []:
                best = find_best_words(contig_blasts, lengths[prev])
                conv_dt[best[0]] += [best[1]]
            prev = row[0]
            contig_blasts = [line]
            conv_dt[row[0]] = [row[1], row[3]]
    best = find_best_words(contig_blasts, lengths[prev])
    conv_dt[best[0]] += [best[1]]
    file.close()
    return conv_dt


def find_best_words(contig_blasts, length):
    """
       Description: Finds the Best Words match per contig
       Input:       Set of hits for contig, Contig Length
       Return:      Best Words Contig
    """
    badwords = ['and', 'or', 'to' 'similar', 'protein', 'predicted', 'hypothetical', 'quality']
    #stuff
    weights = {}
    words_per_contig = []
    for line in contig_blasts:
        line = line.strip('\n').split('\t')
        percent_identity = float(line[2])
        query_coverage = float(line[3]) / length
        words = []
        for word in line[12].split():
            if '[' in word:
                break
            word = word.lower().replace(',','').replace(':','').replace('-like','')
            if word in badwords:
                continue
            if words:
                if (re.compile('^[a-zA-z]{1}$').match(word) or re.compile('^[0-9]+$').match(word)) or 'isoform' in words[-1] and 'x' in word:
                    words[-1] = words[-1] + ' ' + word
            else:
                words.append(word)
        for word in words:
            try:
                weights[word] += percent_identity * query_coverage
            except:
                weights[word] = percent_identity * query_coverage
        words_per_contig.append(words)
    best = 0
    best_index = 0
    for index, line in enumerate(words_per_contig):
        score = sum([weights[word] for word in line])
        if best < score:
            best = score
            best_index = index
    return contig_blasts[best_index].strip('\n').split('\t')


"""Question for Tessa: Need to account for duplicate contigs?"""
def conversion_contig_closest(conversion_file):
    """
       Description:    Contig -> Closest Hit dictionary creator for 
                       database conversions
       Input:        Conversion File
       Return:        Dictionary
    """
    conv_dt = {}
    conversion = open(conversion_file)
    prev = ""
    for line in conversion:
        line = line.strip().split("\t")
        if prev == line[0]:
            continue
        first = line[0]
        second = line[1]
        ###
        #Blast Length
        #third = line[14]
        third = '.'
        ###
        conv_dt[first] = (second, third)
        prev = first
    return conv_dt


def conversion_ez_path(conversion_file):
    """
       Description:    Enzyme -> Kegg Pathway dictionary creator for
                       database conversions
       Input:        Conversion File
       Return:        Dictionary
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
       Description:    GO -> GO Slim dictionary creator for database 
                       conversions
       Input:        Conversion File
       Return:        Dictionary
    """
    conv_dt = {}
    conversion = open(conversion_file)
    for line in conversion:
        if line[0:7] == "id: GO:":
            conv_dt[line[4:].strip()] = line[4:].strip()
    return conv_dt


def conversion_idmap(conversion_file):
    """
       Description:    Generic dictionary creator for IDMAP database 
                       conversions
       Input:        Conversion File
       Return:        Dictionary
    """
    conv_dt = {}
    conversion = open(conversion_file)
    for line in conversion:
        line = line.strip("\n").split("\t")
        conv_dt[line[0]] = line[2].replace(",", ";")
    return conv_dt


def conversion_go_entrez(conversion_file):
    """
       Description:    SwissProt -> GO & Entrez dictionary creator for
                       database conversions
       Input:        Conversion File
       Return:        Dictionary
    """
    conv_dt = {}
    conversion = open(conversion_file)
    for line in conversion:
        line = line.split("\t")
        temp = line[6].replace(" ", "")
        conv_dt[line[0]] = (temp, line[2])
    return conv_dt

