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

from finish_table import *
from init_table import *
import argparse


#MAIN
def main(args):

    indexDict = {'Transcript_id': 0, 'Gene_id': 1, 'Transcript_Length': 2, 'Sprot_BLASTX': 3, 'Sprot_BLASTX_Length': 4,
                 'Prot_id': 5, 'Prot_coordinates': 6, 'Sprot_BLASTP': 7, 'PFAM': 8, 'SignalP': 9, 'TmHMM': 10, 'RNAMMER': 11,
                 'Uniref90_BLASTX': 12, 'Uniref90_BLASTP': 13, 'eggNOG': 14, 'eggNOG_function': 15, 'GO_Terms': 16,
                 'GO_Slim': 17, 'Kegg_Orthology': 18, 'Kegg_Enzyme': 19, 'Kegg_Pathway': 20, 'BLAST_NR': 21,
                 'BLAST_NR_Length': 22, 'BLAST_NR_BestWords': 23, 'Closest_BLASTX': 24, 'Closest_BLASTX_Length': 25,
                 'orthoDB': 26, 'BioCyc': 27, 'entrezGene': 28}


    annot_table = add_fasta(args.fasta, args.geneTransMap)
    annot_table = add_swissprot_top_blastx(annot_table, args.spX)
    annot_table = add_protid(annot_table, args.transdecoder)
    annot_table = add_swissprot_top_blastp(annot_table, args.spP)
    if args.pfam == 'NONE':
        annot_table = add_space(annot_table)
    else:
        annot_table = add_pfam(annot_table, args.pfam)
    if args.signalP == 'NONE':
        annot_table = add_space(annot_table)
    else:
        annot_table = add_signalp(annot_table, args.signalP)
    if args.tmhmm == 'NONE':
        annot_table = add_space(annot_table)
    else:
        annot_table = add_tmhmm(annot_table, args.tmhmm)
    if args.rnammer == 'NONE':
        annot_table = add_space(annot_table)
    else:
        annot_table = add_rnammer(annot_table, args.rnammer)
    if args.ur90X == 'NONE':
        annot_table = add_space(annot_table)
    else:
        annot_table = add_trembl_top_blastx(annot_table, args.ur90X)
    if args.ur90P == 'NONE':
        annot_table = add_space(annot_table)
    else:
        annot_table = add_trembl_top_blastp(annot_table, args.ur90P)
    #INIT_FINISHED
    annot_table = finish_annot_table(annot_table, args, indexDict)
    #FULL_TABLE_BUILT
    annot_by_gene = splitByGene(annot_table)
    #WRITING_FILES
    writeTable(annot_table, args.outfile + '_annotation.txt')
    writeDict(annot_by_gene, args.outfile + '_annotation_by_gene.txt')
    #writeTable(trinotateAnnot, opts.outBase + '_annotation_only_hits.txt', False)

def read_file_test(fname):
    f = argparse.FileType('r')(fname)
    f.close()
    return fname

if __name__ == '__main__':
    psr = argparse.ArgumentParser(description='this module does something')

    #build table files
    psr.add_argument('--fasta',metavar='fasta',type=read_file_test,help='fasta assembly file')
    psr.add_argument('--geneTransMap',metavar='ident',type=read_file_test,help='what is this')
    psr.add_argument('--spX',metavar='ident',type=read_file_test,help='what is this')
    psr.add_argument('--ur90X',metavar='ident',help='what is this')
    psr.add_argument('--rnammer',metavar='ident',type=read_file_test,help='what is this')
    psr.add_argument('--transdecoder',metavar='ident',type=read_file_test,help='what is this')
    psr.add_argument('--spP',metavar='ident',type=read_file_test,help='what is this')
    psr.add_argument('--ur90P',metavar='ident',help='what is this')
    psr.add_argument('--pfam',metavar='ident',type=read_file_test,help='what is this')
    psr.add_argument('--signalP',metavar='ident',type=read_file_test,help='what is this')
    psr.add_argument('--tmhmm',metavar='ident',type=read_file_test,help='what is this')
    #conversion files
    psr.add_argument('--ko2path',metavar='ko2path',type=read_file_test,help='what is this')
    psr.add_argument('--sp2enzyme',metavar='sp2enzyme',type=read_file_test,help='what is this')
    psr.add_argument('--enzyme2path',metavar='enzyme2path',type=read_file_test,help='what is this')
    psr.add_argument('--pfam2enzyme',metavar='pfam2enzyme',type=read_file_test,help='what is this')
    psr.add_argument('--go2path',metavar='go2path',type=read_file_test,help='what is this')
    psr.add_argument('--nog2function',metavar='nog2function',type=read_file_test,help='what is this')
    psr.add_argument('--contig2closest',metavar='contig2closest',type=read_file_test,help='what is this')
    psr.add_argument('--go2slim',metavar='go2slim',type=read_file_test,help='what is this')
    psr.add_argument('--contig2blastnr',metavar='contig2blastnr',type=read_file_test,help='what is this')
    #IDMAPPING
    psr.add_argument('--sp2ko',metavar='sp2ko',type=read_file_test,help='what is this')
    psr.add_argument('--sp2nog',metavar='sp2nog',type=read_file_test,help='what is this')
    psr.add_argument('--sp2ortho',metavar='sp2ortho',type=read_file_test,help='what is this')
    psr.add_argument('--sp2bioc',metavar='sp2bioc',type=read_file_test,help='what is this')
    #IDMAPPING SELECTED
    psr.add_argument('--sp2goentrez',metavar='sp2goentrez',type=read_file_test,help='what is this')
    #outfile name
    psr.add_argument('--outfile',metavar='outfile',help='what is this')

    args = psr.parse_args()
    main(args)
