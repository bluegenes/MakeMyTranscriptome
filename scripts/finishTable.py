###############################################################################
###############################################################################
#
# Author	-	Andrew Walters
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

import sys, re, os
import urllib2
import operator
from buildDicts import *
from addColumns import *
#from list_of_dictionaries import build_dictionaries, build_weights




#Changed so handles unique annotations
def splitByGene(annotTable):
    geneDt = {}
    for line in annotTable:
        gene = line[1]
        if gene in geneDt:
            for column, item in enumerate(geneDt[gene]):
            	if column > 0:
            	    tableColumn = column + 1
            	else:
            	    tableColumn = column
            	entries = item.split(";;")
        	if line[tableColumn] not in entries and line[tableColumn] != ".":
        	    geneDt[gene][column] = geneDt[gene][column] + ";;" + line[tableColumn]
        else:
            geneDt[gene] = [line[0]] + line[2:] # don't need gene ID in the entry
    return geneDt

"""Returns dictionary comprised of only entries with a specific column filled and the int number of total entries
   with that column"""
def extractAnnotations(annotDict, colToExtract):
    extractedDt = {}
    count = 0
    for gene, info in annotDict.items():
        if len(info[colToExtract]) > 1:
            extractedDt[gene] = info # if has hit in this column
            count += 1
    return extractedDt, count #fiName

"""
def writeTable(table, outfile):

    header = ['Transcript_id', 'Gene_id', 'Transcript_Length', 'Sprot_BLASTX', 'Sprot_BLASTX_Length', 'Prot_id',
          'Prot_coordinates', 'Sprot_BLASTP', 'PFAM', 'SignalP', 'TmHMM', 'RNAMMER', 'Uniref90_BLASTX',
          'Uniref90_BLASTP', 'eggNOG', 'eggNOG_function', 'GOer_Terms', 'GO_slim', 'Kegg_Orthology', 'Kegg_Enzyme',
          'Kegg_Pathway', 'BLAST_NR', 'BLAST_NR_Length', 'BLAST_NR_BestWords', 'Closest_BLASTX',
          'Closest_BLASTX_Length', 'orthoDB', 'BioCyc', 'entrezGene']

    outFile = open(outfile + '_annotation.txt', 'w')
    outFileHit = open(outfile + '_only_hits_annotation.txt', 'w')
    outFileTr = open(outfile + '_annotation_longest_transcript.txt', 'w')
    outFileOrf = open(outfile + '_annotation_longest_orf.txt', 'w')

    outFile.write("\t".join(header) + "\n")
    outFileHit.write("\t".join(header) + "\n")
    outFileTr.write("\t".join(header) + "\n")
    outFileOrf.write("\t".join(header) + "\n")

    contigID = table[0][0]
    longestTr = int(table[0][23])
    topTr = table[0]
    prevTr = []
    longestOrf = 0
    topOrf = []
    prevOrf = []
    for line in table:
        outFile.write('\t'.join(map(str,line)) + '\n')
        #Only Hits
        if lineGood(line):
            outFileHit.write('\t'.join(map(str,line)) + '\n')
        #Longest Transcript
        if contigID != line[0] and topTr != prevTr:
            outFileTr.write('\t'.join(map(str,topTr)) + '\n')
            prevTr = topTr
            longestTr = 0
        if longestTr < int(line[23]):
            topTr = line
            longestTr = int(line[23])
        #Longest Orf
        if contigID != line[0] and topOrf != prevOrf:
            outFileOrf.write('\t'.join(map(str,topOrf)) + '\n')
            prevOrf = topOrf
            longestOrf = 0
        if line[6] != ".":
            length = (line[6].split("[")[0]).split("-")
            length = int(length[1]) - int(length[0])
            if longestOrf < length:
                topOrf = line
                longestOrf = length
        contigID = line[0]
    outFile.close()
"""

def writeTable(table, outfile):
    header = ['Transcript_id', 'Gene_id', 'Transcript_Length', 'Sprot_BLASTX', 'Sprot_BLASTX_Length', 'Prot_id',
          'Prot_coordinates', 'Sprot_BLASTP', 'PFAM', 'SignalP', 'TmHMM', 'RNAMMER', 'Uniref90_BLASTX',
          'Uniref90_BLASTP', 'eggNOG', 'eggNOG_function', 'GO_terms', 'GO_slim', 'Kegg_Orthology', 'Kegg_Enzyme',
          'Kegg_Pathway', 'BLAST_NR', 'BLAST_NR_Length', 'BLAST_NR_BestWords', 'Closest_BLASTX',
          'Closest_BLASTX_Length', 'orthoDB', 'BioCyc', 'entrezGene']

    outFile = open(outfile, 'w')

    outFile.write("\t".join(header) + "\n")

    for line in table:
        outFile.write('\t'.join(map(str,line)) + '\n')
    outFile.close()


def lineGood(line):
    #not checking first two (geneID & trascript ID) because they always exist
    for index, entry in enumerate(line[2:]):
        if entry != '.':
            #use 21 instead of 23 to ignore transcript length because it always exists (21+2 = 23)
            if index != 21:
                return True
    return False

def writeDict(dict, output_file):
    header = ['Gene_id', 'Transcript_id', 'Transcript_Length', 'Sprot_BLASTX', 'Sprot_BLASTX_Length', 'Prot_id',
          'Prot_coordinates', 'Sprot_BLASTP', 'PFAM', 'SignalP', 'TmHMM', 'RNAMMER', 'Uniref90_BLASTX',
          'Uniref90_BLASTP', 'eggNOG', 'eggNOG_function', 'GO_Terms', 'GO_slim', 'Kegg_Orthology', 'Kegg_Enzyme',
          'Kegg_Pathway', 'BLAST_NR', 'BLAST_NR_Length', 'BLAST_NR_BestWords', 'Closest_BLASTX',
          'Closest_BLASTX_Length', 'orthoDB', 'BioCyc', 'entrezGene']

    sortByGene = sorted(dict.items(), key=operator.itemgetter(0))
    outFile = open(output_file, 'w')

    outFile.write("\t".join(header) + "\n")

    for (key, val) in sortByGene:
        outFile.write(key + '\t' + '\t'.join(map(str,val)) + '\n')
    outFile.close()


"""Main method to finish full annotation table"""
def finishAnnotTable(annot_table, args, indexDict):

    # create conversion dicts and add kegg annotations

    #koToPathDt = orthology_pathway.list
    #spToEzDt = swiss_enzyme.list
    #ezToPathDt = enzyme_pathway.list
    #pfToEzDt = pfam_enzyme.list
    #goToPath = go_pathway.txt
    #nogToF = allKOG_functional_info.txt
    #contigToClosest = contig2closest
    #goToGOSlim = goslim_generic.obo
    #contigToBlastNR = contig2blastnr.txt
    #spToKO = idmapping.KO
    #spToNog = idmapping.eggNOG
    #spToOrtho = idmapping.orthodb
    #spToBioC = idmapping.biocyc
    #spToGOandEntrez = idmapping_selected.tab


    """Conversion Files"""
    koToPathDt = createConversion(args.ko2path)
    spToEzDt = createConversion(args.sp2enzyme)
    ezToPathDt = createConversionEzPath(args.enzyme2path)
    pfToEzDt = createConversion(args.pfam2enzyme)
    goToPath = createConversion(args.go2path)
    nogToF = createConversionNogF(args.nog2function)
    contigToClosest = createConversionContigClosest(args.contig2closest)
    goToGOSlim = createConversionGOSlim(args.go2slim)
    contigToBlastNR = createConversionContigBlastNR(args.contig2blastnr)
    spToKO = createConversionIDMAP(args.sp2ko)
    spToNog = createConversionIDMAP(args.sp2nog)
    spToOrtho = createConversionIDMAP(args.sp2ortho)
    spToBioC = createConversionIDMAP(args.sp2bioc)
    spToGOandEntrez = createConversoinGOandEntrez(args.sp2goentrez)


    for index, line in enumerate(annot_table):
        annot_table[index] = line + (['.'] * 15)
        annot_table[index] = addIDMAP(annot_table[index], spToKO, indexDict, indexDict['Kegg_Orthology'])
        annot_table[index] = addIDMAP(annot_table[index], spToNog, indexDict, indexDict['eggNOG'])
        annot_table[index] = addIDMAP(annot_table[index], spToOrtho, indexDict, indexDict['orthoDB'])
        annot_table[index] = addIDMAP(annot_table[index], spToBioC, indexDict, indexDict['BioCyc'])
        annot_table[index] = addSPtoGOandEntrez(annot_table[index], spToGOandEntrez, indexDict)
        annot_table[index] = addKOtoPath(annot_table[index], koToPathDt, indexDict)
        annot_table[index] = addSPtoEnzymetoPath(annot_table[index], spToEzDt, ezToPathDt, indexDict)
        annot_table[index] = addPfamtoEnzymetoPath(annot_table[index], pfToEzDt, ezToPathDt, indexDict)
        annot_table[index] = addGOtoPath(annot_table[index], goToPath, indexDict)
        annot_table[index] = addNogtoF(annot_table[index], nogToF, indexDict)
        annot_table[index] = addBlastNR(annot_table[index], contigToBlastNR, indexDict)
        annot_table[index] = addClosestHit(annot_table[index], contigToClosest, indexDict)
        annot_table[index] = addGOSlim(annot_table[index], goToGOSlim, indexDict)
    return annot_table
