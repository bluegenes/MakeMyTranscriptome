from finishTable import *


def getConversion(lookup, conversDt):
    ''' Takes in a conversoin dictionary to convert from one database to another'''
    try:
        return conversDt[lookup]
    except:
        return '.'


def addKOtoPath(line, koToPathDt, index):
    ''' Finds KO Pathway using KO number and adds it to column 20
        Kegg_Orthology:18 -> Kegg_Pathway:20
    '''
    line[index['Kegg_Pathway']] = getConversion(line[index['Kegg_Orthology']], koToPathDt)
    return line


def addSPtoEnzymetoPath(line, spToEzDt, ezToPathDt, index):
    '''Finds Kegg Enzyme using SP number and adds it to column 19
       Sprot_BLASTX:3 -> Kegg_Enzyme:19
       Finds KO Pathway using Kegg Enzyme if no Pathway present
       Kegg_Enzyme:19 -> Kegg_Pathway:20
    '''
    line[index['Kegg_Enzyme']] = getConversion(line[index['Sprot_BLASTX']], spToEzDt)
    if line[index['Kegg_Enzyme']] != '.' and line[index['Kegg_Pathway']] == '.':
        line[index['Kegg_Pathway']] = getConversion(line[index['Kegg_Enzyme']], ezToPathDt)
    return line


def addPfamtoEnzymetoPath(line, pfToEzDt, ezToPathDt, index):
    ''' Finds Kegg Enzyme using PFAM and adds it to column 19
        PFAM:8 -> Kegg_Enzyme:19
        Finds KO Pathway using Kegg Enzyme if no Pathway present
        Kegg_Enzyme:19 -> Kegg_Pathway:20
    '''
    if line[index['Kegg_Enzyme']] == '.':
        pfam = re.compile('^PF(\S*)\^')
        hasPFmatch = pfam.match(line[index['PFAM']])
        if hasPFmatch is not None:
            pfs = hasPFmatch.groups()
            for pf in pfs:
                line[index['Kegg_Enzyme']] = getConversion(pf, pfToEzDt)
                if line[index['Kegg_Enzyme']] != '.':
                    break
            if line[index['Kegg_Enzyme']] != '.' and line[index['Kegg_Pathway']] == '.':
                line[index['Kegg_Pathway']] = getConversion(line[index['Kegg_Enzyme']], ezToPathDt)
    return line


def addGOtoPath(line, goToPath, index):
    ''' Finds Kegg Pathway using GO Terms and adds it to column 20
        GO_Terms:16 -> Kegg_Pathway:20
    '''
    go_match = re.compile('(GO:\d*)')
    go_terms = re.findall(go_match, line[index['GO_Terms']])
    if go_terms != []:
        if line[index['Kegg_Pathway']] == '.':
            for go in go_terms:
                if line[index['Kegg_Pathway']] == '.':
                    line[index['Kegg_Pathway']] = getConversion(go[3:], goToPath)
                else:
                    line[index['Kegg_Pathway']] = ';'.join([line[index['Kegg_Pathway']], getConversion(go[3:], goToPath)])
    return line


def addNogtoF(line, nogToF, index):
    ''' Finds eggNOG Function using eggNOG and adds it to column 15
        eggNOG:14 -> eggNOG_Function:15
    '''
    nog_match = re.compile('^(\S*)\^', re.IGNORECASE)
    nog_term = re.search(nog_match, line[index['eggNOG']])
    if nog_term is not None:
        nog_term = nog_term.groups()[0]
        line[index['eggNOG_Function']] = getConversion(nog_term, nogToF)
    return line


def addBlastNR(line, contigToBlastNR, index):
    ''' Finds BlastNR using Transcript ID and adds relative info to columns 21,22,23
        Transcript_id:0 -> BLAST_NR:21 + BLAST_NR_Length:22 + BLAST_NR_BestWords:23
    '''
    contig = line[index['Transcript_id']]
    info = getConversion(contig, contigToBlastNR)
    if info != ".":
        line[index['BLAST_NR']] = info[0]
        line[index['BLAST_NR_Length']] = info[1]
        line[index['BLAST_NR_BestWords']] = info[2]
    return line


def addClosestHit(line, contigToClosest, index):
    ''' Finds Closest BlastX Hit using Transcript ID and adds relative info to columns 24,25
        Transcript_id:0 -> Closest_BLASTX:24 + Closest_BLASTX_Length:25
    '''
    contig = line[index['Transcript_id']]
    best = getConversion(contig, contigToClosest)
    if best != ".":
        line[index['Closest_BLASTX']] = best[0]
        line[index['Closest_BLASTX_Length']] = best[1]
    return line


def addGOSlim(line, goToGOSlim, index):
    ''' Finds GO Slim Terms using GO Terms and adds it to column 17
        GO_Terms:16 -> GO_Slim:17
    '''
    goTerms = line[index['GO_Terms']].split(";")
    slim = []
    for g in goTerms:
        check = getConversion(g, goToGOSlim)
        if check != ".":
            slim.append(check)
    if slim == []:
        return line
    line[index['GO_Slim']] = ";".join(slim)
    return line


def addIDMAP(line, convDt, index, column):
    ''' Used for idmapping.dat dependency columns
        Sprot_BLASTX:3 -> Kegg_Orthology:18
        Sprot_BLASTX:3 -> eggNOG:14
        Sprot_BLASTX:3 -> orthoDB:26
        Sprot_BLASTX:3 -> BioCyc:27
    '''
    if line[index['Sprot_BLASTX']] != ".":
        spID = line[index['Sprot_BLASTX']].split("|")[1]
        line[column] = getConversion(spID, convDt)
    return line


def addSPtoGOandEntrez(line, spToGOandEntrez, index):
    ''' Used for idmapping_selected dependency columns
        Sprot_BLASTX:3 -> GO_Terms:16
        Sprot_BLASTX:3 -> entrezGene:28
    '''
    if line[index['Sprot_BLASTX']] != ".":
        spID = line[index['Sprot_BLASTX']].split("|")[1]
        temp = getConversion(spID, spToGOandEntrez)
	if temp != ".":
            line[index['GO_Terms']] = temp[0]
            line[index['entrezGene']] = temp[1]
	else:
            line[index['GO_Terms']] = temp
            line[index['entrezGene']] = temp
    return line
