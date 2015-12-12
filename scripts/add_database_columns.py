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

def get_conversion(key, conversion_dt):
    ''' Takes in a conversoin dictionary to convert from one database to another'''
    try:
        return conversion_dt[key]
    except:
        return '.'


def add_ko_to_path(line, ko_path_dt, index):
    ''' Finds KO Pathway using KO number and adds it to Column 20
        Kegg_Orthology:18 -> Kegg_Pathway:20
    '''
    line[index['Kegg_Pathway']] = get_conversion(line[index['Kegg_Orthology']], ko_path_dt)
    return line


def add_sp_to_enzyme_to_path(line, sp_ez_dt, ez_path_dt, index):
    '''Finds Kegg Enzyme using SP number and adds it to Column 19
       Sprot_BLASTX:3 -> Kegg_Enzyme:19
       Finds KO Pathway using Kegg Enzyme if no Pathway present
       Kegg_Enzyme:19 -> Kegg_Pathway:20
    '''
    line[index['Kegg_Enzyme']] = get_conversion(line[index['Sprot_BLASTX']], sp_ez_dt)
    if line[index['Kegg_Enzyme']] != '.' and line[index['Kegg_Pathway']] == '.':
        line[index['Kegg_Pathway']] = get_conversion(line[index['Kegg_Enzyme']], ez_path_dt)
    return line


def add_pfam_to_enzyme_to_path(line, pf_ez_dt, ez_path_dt, index):
    ''' Finds Kegg Enzyme using PFAM and adds it to Column 19
        PFAM:8 -> Kegg_Enzyme:19
        Finds KO Pathway using Kegg Enzyme if no Pathway present
        Kegg_Enzyme:19 -> Kegg_Pathway:20
    '''
    if line[index['Kegg_Enzyme']] == '.':
        pfam = re.compile('^PF(\S*)\^')
        pfam_match = pfam.match(line[index['PFAM']])
        if pfam_match is not None:
            pfs = pfam_match.groups()
            for pf in pfs:
                line[index['Kegg_Enzyme']] = get_conversion(pf, pf_ez_dt)
                if line[index['Kegg_Enzyme']] != '.':
                    break
            if line[index['Kegg_Enzyme']] != '.' and line[index['Kegg_Pathway']] == '.':
                line[index['Kegg_Pathway']] = get_conversion(line[index['Kegg_Enzyme']], ez_path_dt)
    return line


def add_go_to_path(line, go_path_dt, index):
    ''' Finds Kegg Pathway using GO Terms and adds it to Column 20
        GO_Terms:16 -> Kegg_Pathway:20
    '''
    go_match = re.compile('(GO:\d*)')
    go_terms = re.findall(go_match, line[index['GO_Terms']])
    if go_terms != []:
        if line[index['Kegg_Pathway']] == '.':
            for go in go_terms:
                if line[index['Kegg_Pathway']] == '.':
                    line[index['Kegg_Pathway']] = get_conversion(go[3:], go_path_dt)
                else:
                    line[index['Kegg_Pathway']] = ';'.join([line[index['Kegg_Pathway']], get_conversion(go[3:], go_path_dt)])
    return line


def add_nog_to_func(line, nog_func_dt, index):
    ''' Finds eggNOG Function using eggNOG and adds it to Column 15
        eggNOG:14 -> eggNOG_Function:15
    '''
    nog_term = line[index['eggNOG']]
    if nog_term != ".":
        line[index['eggNOG_function']] = get_conversion(nog_term, nog_func_dt)
    return line


def add_blastnr(line, contig_blastnr_dt, index):
    ''' Finds BlastNR using Transcript ID and adds relative info to Columns 21,22,23
        Transcript_id:0 -> BLAST_NR:21 + BLAST_NR_Length:22 + BLAST_NR_BestWords:23
    '''
    contig = line[index['Transcript_id']]
    info = get_conversion(contig, contig_blastnr_dt)
    if info != ".":
        line[index['BLAST_NR']] = info[0]
        line[index['BLAST_NR_Length']] = info[1]
        line[index['BLAST_NR_BestWords']] = info[2]
    return line


def add_closest_hit(line, contig_closest_dt, index):
    ''' Finds Closest BlastX Hit using Transcript ID and adds relative info to Columns 24,25
        Transcript_id:0 -> Closest_BLASTX:24 + Closest_BLASTX_Length:25
    '''
    contig = line[index['Transcript_id']]
    best = get_conversion(contig, contig_closest_dt)
    if best != ".":
        line[index['Closest_BLASTX']] = best[0]
        line[index['Closest_BLASTX_Length']] = best[1]
    return line


def add_goslim(line, go_goslim_dt, index):
    ''' Finds GO Slim Terms using GO Terms and adds it to Column 17
        GO_Terms:16 -> GO_Slim:17
    '''
    go_terms = line[index['GO_Terms']].split(";")
    slim = []
    for g in go_terms:
        check = get_conversion(g, go_goslim_dt)
        if check != ".":
            slim.append(check)
    if slim == []:
        return line
    line[index['GO_Slim']] = ";".join(slim)
    return line


def add_idmap(line, conversion_dt, index, column):
    ''' Used for idmapping.dat dependency Columns
        Sprot_BLASTX:3 -> Kegg_Orthology:18
        Sprot_BLASTX:3 -> eggNOG:14
        Sprot_BLASTX:3 -> orthoDB:26
        Sprot_BLASTX:3 -> BioCyc:27
    '''
    if line[index['Sprot_BLASTX']] != ".":
        swissprot_id = line[index['Sprot_BLASTX']].split("|")[1]
        line[column] = get_conversion(swissprot_id, conversion_dt)
    return line


def add_sp_to_go_and_entrez(line, sp_go_entrez_dt, index):
    ''' Used for idmapping_selected dependency Columns
        Sprot_BLASTX:3 -> GO_Terms:16
        Sprot_BLASTX:3 -> entrezGene:28
    '''
    if line[index['Sprot_BLASTX']] != ".":
        swissprot_id = line[index['Sprot_BLASTX']].split("|")[1]
        temp = get_conversion(swissprot_id, sp_go_entrez_dt)
	if temp != ".":
            line[index['GO_Terms']] = temp[0]
            line[index['entrezGene']] = temp[1]
	else:
            line[index['GO_Terms']] = temp
            line[index['entrezGene']] = temp
    return line
