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

import sys, re, os
import operator
from build_conversions import *
from add_database_columns import *
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


def finish_annot_table(annot_table, args, index_dict, lengths):     
    """
       Description: Main method to add database conversions
       Input:       Annotation Table, Input Args, Column Index Dictionary
       Return:      Annotation Table
    """
    # Conversion Files
    ko_to_path = conversion(args.ko2path)
    sp_to_ez = conversion(args.sp2enzyme)
    ez_to_path = conversion_ez_path(args.enzyme2path)
    pf_to_ez = conversion(args.pfam2enzyme)
    go_to_path = conversion(args.go2path)
    nog_to_f = conversion_nogf(args.nog2function)
    if args.closest != None:
        contig_to_closest = conversion_contig_closest(args.closest)
    go_to_goslim = conversion_goslim(args.go2slim)
    if args.nr != None:
        contig_to_blastnr = conversion_contig_blastnr(args.nr, lengths)
    sp_to_ko = conversion_idmap(args.sp2ko)
    sp_to_nog = conversion_idmap(args.sp2nog)
    sp_to_ortho = conversion_idmap(args.sp2ortho)
    sp_to_bioc = conversion_idmap(args.sp2bioc)
    sp_to_go_entrez = conversion_go_entrez(args.sp2goentrez)


    for index, line in enumerate(annot_table):
        annot_table[index] = add_idmap(annot_table[index], sp_to_ko, index_dict, index_dict['Kegg_Orthology'])
        annot_table[index] = add_idmap(annot_table[index], sp_to_nog, index_dict, index_dict['eggNOG'])
        annot_table[index] = add_idmap(annot_table[index], sp_to_ortho, index_dict, index_dict['orthoDB'])
        annot_table[index] = add_idmap(annot_table[index], sp_to_bioc, index_dict, index_dict['BioCyc'])
        annot_table[index] = add_sp_to_go_and_entrez(annot_table[index], sp_to_go_entrez, index_dict)
        annot_table[index] = add_ko_to_path(annot_table[index], ko_to_path, index_dict)
        annot_table[index] = add_sp_to_enzyme_to_path(annot_table[index], sp_to_ez, ez_to_path, index_dict)
        annot_table[index] = add_pfam_to_enzyme_to_path(annot_table[index], pf_to_ez, ez_to_path, index_dict)
        annot_table[index] = add_go_to_path(annot_table[index], go_to_path, index_dict)
        annot_table[index] = add_nog_to_func(annot_table[index], nog_to_f, index_dict)
        if args.nr != None:
            annot_table[index] = add_blastnr(annot_table[index], contig_to_blastnr, index_dict)
        if args.closest != None:
            annot_table[index] = add_closest_hit(annot_table[index], contig_to_closest, index_dict)
        annot_table[index] = add_goslim(annot_table[index], go_to_goslim, index_dict)
    return annot_table
