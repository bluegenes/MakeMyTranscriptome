import pandas as pd
import argparse

def get_transcript_length(fasta):
    with open(fasta, 'r') as f:
        lengths = {}
        length = 0
        for line in f:
            if line.startswith('>'):
                if length > 0:
                    lengths[transcriptid] = length
                    length = 0
                transcriptid = line.split(' ')[0][1:]
            else:
                length += len(line.strip())
        lengths[transcriptid] = length
    lenSeries = pd.Series(lengths, name="Transcript_Length")
    return(lenSeries)

def get_orf_info(orf_bed):
    #was: 'Prot_id': 5, 'Prot_coordinates': 6, NOW: longestORF id; ORF_length
    orfInfo = pd.read_table(orf_bed, skiprows=1, header=None) # skip first row
    orfInfo['ORF_length'] = orfInfo.iloc[:,3].str.extract('len:(\d*)')
    orfInfo.rename(columns = {0:'Transcript_id',3: 'Longest_ORF_Info'}, inplace=True)
    orfInfo['Longest_ORF_id'] = orfInfo['Longest_ORF_Info'].str.extract('ID=([^;]*)')
    lorfInfo = orfInfo.iloc[orfInfo.groupby(['Transcript_id']).apply(lambda x: x['ORF_length'].idxmax())]
    longORFinfo  = lorfInfo.loc[:,['Transcript_id', 'Longest_ORF_Info', 'Longest_ORF_id', 'ORF_length']]
    longORFinfo.set_index(['Transcript_id'], inplace=True)
    longORFinfo.drop_duplicates(keep='first', inplace=True)
    return(longORFinfo)

def get_blast_info(query_lengths, blastFile, blastType='blastx', blastDB = 'sp', bestWords=False):
    blastInfo = pd.read_table(blastFile, header=0)
    blastInfo.loc[:,'hit_id'] = blastInfo.loc[:, 'subject_id'].astype(str) + '_' +  blastInfo.loc[:, 'full_name'].astype(str)
    if bestWords:
        bestWordsDF = blastInfo.copy()
        query_lengthsDF = query_lengths.to_frame() #calculate query coverage
        if blastType == 'blastx':
            query_lengthsDF.rename(columns={'Transcript_Length': 'query_length'}, inplace=True)
        else:
            query_lengthsDF.rename(columns={'ORF_Length': 'query_length'}, inplace=True)
        bestWordsDF = bestWordsDF.merge(query_lengthsDF,  how='outer',left_on='query_id',right_index=True).dropna()
        bestWordsDF.loc[:,'query_coverage'] = bestWordsDF.loc[:,'alignment_length'].astype(float)/ bestWordsDF['query_length'].astype(float)
        # get the words from each entry
        bestWordsDF.loc[:,'words'] = bestWordsDF.loc[:,'full_name'].str.strip().str.lower().str.replace(',','').str.replace(':','').str.replace('-like','').str.replace(' isoform','_isoform').str.replace('\[\S*', '').str.replace('\S*\]', '')
        bestWordsDF.loc[:,'words'] = bestWordsDF.loc[:,'words'].str.replace('and','').str.replace('or', '').str.replace('similar', '').str.replace('protein','').str.replace('hypothetical', '').str.replace('quality','').str.replace('unnamed','').str.replace('product','').str.split()
        bestWordsDF.loc[:,'multiplier'] = bestWordsDF.loc[:,'percent_identity'].astype(float) * bestWordsDF.loc[:,'query_coverage'].astype(float)
       # this is no *quite* doing what we want --> want to use all 20 matches to get scores, then assess that way. right now, we're just picking the top total # words really
        bestWordsDF.loc[:,'wordDt']= bestWordsDF.apply(wordsToScoreDt, axis = 1)
        bestWordsDF.loc[:,'score'] = bestWordsDF.apply(get_total_score, axis = 1)
        #get the best word score
        bestWordScore = bestWordsDF.iloc[bestWordsDF.groupby([cNames[0]]).apply(lambda x: x['score'].idxmax())]
        bestWordScore.set_index(cNames[0], inplace=True)
        bestWordScore.drop_duplicates(keep='first', inplace=True)
        #here we need to add a column for best word hit name
        #infoToAdd['best_word'] = bestWordScore.loc[:, 'hit_id']
    blastInfo = blastInfo.iloc[blastInfo.groupby(['query_id']).apply(lambda x: x['evalue'].idxmin())]
    blastInfo.set_index('query_id', inplace=True)
    blastInfo.drop_duplicates(keep='first',inplace=True)
    blastInfo= blastInfo.loc[:, ['hit_id', 'subject_length']]
    query_name, db_name = get_colnames(blastType,blastDB)
    hit_name = db_name + '_' + blastType
    # change column names to reflect the database + the blast type
    blastInfo.rename(columns={'query_id': query_name, 'hit_id': hit_name, 'subject_length': hit_name + '_length'}, inplace=True)
    return(blastInfo)

def wordsToScoreDt(df_row):
    wordlist = list(set(df_row['words']))
    numWords = len(wordlist)
    wordDt = dict(zip(wordlist, numWords*([df_row['multiplier']])))
    return wordDt

def get_total_score(df_row):
    wordlist = list(set(df_row['words']))
    numWords = len(wordlist)
    total_score = numWords* df_row['multiplier']
    return total_score

def get_colnames(blastType,db):
    query_id,db_id = '',''
    if blastType == 'blastx':
        query_id = 'Transcript_id'
    else:
        query_id = 'ORF_id'
    if db == 'sp':
        db_id = 'swissprot'
    elif db == 'ur90':
        db_id == 'uniref90'
    elif db == 'nr':
        db_id = 'nr'
    return query_id, db_id

def split_pfam_list(pfamRow):
    pfamRow[22] = pfamRow[1] + ':'  + ' iE-val__' + pfamRow[12]  + ': ' + ' '.join(pfamRow[22:])
    pfamRow = '\t'.join(pfamRow[:23])
    return(pfamRow)

def get_pfam_info(pfam_file):
    pfamD = pd.read_table(pfam_file, header=None, comment='#')
    #pfamD= pfamD.iloc[:,0].str.split('\s*').apply(split_pfam_list)
    pfamD= pfamD.iloc[:,0].str.split().apply(split_pfam_list)
    pfamDF = pfamD.to_frame().iloc[:,0].apply(lambda x: pd.Series(x.split('\t')))
    ## to have the description included, replace '1' with '22' in both lines below
    groupORF = pfamDF.groupby(3)[1].apply(lambda x: ','.join(x)).reset_index()
    groupORF.rename(columns={3:'ORF_id', 1:'PFAM'}, inplace=True)
    groupORF.set_index('ORF_id', inplace=True)
    return(groupORF)

def get_signalp(signalp_file):
     sigInfo = pd.read_table(signalp_file, header=None, comment='#')
     sigInfo = sigInfo.rename(columns={0:'ORF_id'})
     sigInfo['signalp'] =  sigInfo.iloc[:,2] + ':' + sigInfo.iloc[:,3:5].apply(lambda x: '-'.join(x.astype(str)), axis=1) + '__score:' +  sigInfo.iloc[:,5].astype(str)
     infoToAdd = sigInfo.loc[:,['ORF_id', 'signalp']]
     infoToAdd.set_index('ORF_id', inplace=True)
     return(infoToAdd)

def get_tmhmm(tmhmm_file, AA_threshold = 1):
    tmInfo = pd.read_table(tmhmm_file, header=None, comment='#')
    tmInfo = tmInfo.rename(columns={0:'ORF_id'})
    tmInfo['tmhmm'] =  tmInfo.iloc[:,[2,4,5]].apply(lambda x: '__'.join(x.astype(str)), axis=1)
    tmInfo['expAA'] =  tmInfo.ix[:,2].str.extract('ExpAA=(\d+\.\d{1,2})').astype(float)
    infoToAdd = tmInfo.loc[tmInfo['expAA'] > AA_threshold, ['ORF_id', 'tmhmm']]
    infoToAdd.set_index('ORF_id', inplace=True)
    return(infoToAdd)


def main(args):
    #initial setup
    initDF = pd.read_table(args.geneTransMap, index_col=1, header=None)
    initDF.rename(columns = {0:'Gene_id'}, inplace=True)
    transcript_lengths = get_transcript_length(args.fasta)
    initDF = initDF.join(transcript_lengths)
    blastXinfo = get_blast_info(transcript_lengths, args.spX, blastDB='sp', blastType='blastx')
    initDF = initDF.join(blastXinfo)
    orfInfo = get_orf_info(args.transdecoder_bed)
    orfLengths = orfInfo.loc[:, ['Longest_ORF_id', 'ORF_length']].set_index('Longest_ORF_id')
    initDF = initDF.join(orfInfo)
    if args.ur90X is not None:
        initDF = initDF.join(get_blast_info(transcript_lengths,args.ur90X, blastDB='ur90', blastType='blastx'))
    if args.nrX is not None:
        initDF = initDF.join(get_blast_info(transcript_lengths, args.nrX, blastDB='nr', blastType='blastx')) #bestWords = True (when this is working)
    #orf-level annotations
    blastpDF = get_blast_info(orfLengths, args.spP, blastDB='sp', blastType='blastp')
    initDF = initDF.merge(blastpDF, how='left', left_on='Longest_ORF_id', right_index=True, copy=False)
    pfamDF = get_pfam_info(args.pfam)
    initDF = initDF.merge(pfamDF, how='left', left_on='Longest_ORF_id', right_index=True, copy=False)
    if args.ur90P is not None:
        ur90blastpDF = get_blast_info(orfLengths,args.ur90P, blastDB='ur90', blastType='blastp')
        initDF = initDF.merge(ur90blastpDF, how='left', left_on='Longest_ORF_id', right_index=True, copy=False)
    if args.nrP is not None:
        nrblastpDF = get_blast_info(orfLengths, args.nrP, blastDB='nr', blastType='blastp') #bestWords = True (when this is working)
        initDF = initDF.merge(nrblastpDF, how='left', left_on='Longest_ORF_id', right_index=True, copy=False)
    if args.signalP is not None:
        sigpDF = get_signalp(args.signalP)
        initDF = initDF.merge(sigpDF, how='left', left_on='Longest_ORF_id', right_index=True, copy=False)
    if args.tmhmm is not None:
        tmhmmDF = get_tmhmm(args.tmhmm)
        initDF = initDF.merge(tmhmmDF,how='left', left_on='Longest_ORF_id', right_index=True, copy=False)
    initDF = get_keggInfo(initDF, args.sp2ko, args.ko2path)
    initDF = get_eggNOG(initDF, args.sp2nog,args.nog2function)
    initDF = addGO(initDF, args.sp2goentrez)
#    initDF = initDF.drop('spHitX', 1)
    #change this --> within the functions where we create these:
    #initDF.drop('spHitX',axis=1, inplace=True)
    #initDF.drop('spHitP',axis=1, inplace=True)
#    summary = initDF.describe().transpose()
#    summary.to_csv(args.outfile+'_annotation_summary', sep='\t')
    initDF.to_csv(args.outfile + '_annotation.txt', index_label='Transcript_id', sep='\t', na_rep='.')
    summary = initDF.count()
#    summary['Transcript_id'] =initDF['Transcript_id'].nunique()
#    summary['Gene_id'] =initDF['Gene_id'].nunique()
#    summary['eggNOG_function'] =initDF['eggNOG_function'].nunique()
#    summary['Kegg_Pathway'] = initDF['Kegg_Pathway'].nunique()
#    summary['Kegg_Orthology '] = initDF['Kegg_Orthology'].nunique()
#    summary['swissprot_blastx_unique'] = initDF['swissprot_blastx'].nunique()
#    summary['swissprot_blastp_unique'] = initDF['swissprot_blastp'].nunique()
#    summary['PFAM_unique'] = initDF['PFAM'].nunique()
    summary.to_json(args.outfile + '_annotation_summary.json')



# now map sp id's to other databases
def get_keggInfo(initDF, sp_to_ko, koToPath):
    spKO = pd.read_table(sp_to_ko, header=None, index_col=0, names= ['ko', 'Kegg_Orthology'])
    spKO.drop('ko', axis=1, inplace=True)
    ko2pathD = pd.read_table(koToPath, header=None, names= ['Kegg_Orthology', 'Kegg_Pathway', 'type'])
    ko2pathD.drop('type', axis=1, inplace=True)
    ko2pathD['Kegg_Orthology'] = ko2pathD['Kegg_Orthology'].str.extract('ko:(\S+)')
    ko2pathD['Kegg_Pathway'] = ko2pathD['Kegg_Pathway'].str.extract('path:(\S+)')
    ko2pathD = ko2pathD.groupby(['Kegg_Orthology'])['Kegg_Pathway'].apply(lambda x: ','.join(x)).reset_index()
    ko2pathD.set_index('Kegg_Orthology', inplace=True)
    spKO = spKO.merge(ko2pathD, how='left', left_on='Kegg_Orthology', right_index=True, copy=False)
    initDF['spHitX'] = initDF['swissprot_blastx'].str.extract('sp\|(\S*)\|')
    initDF = initDF.merge(spKO, how='left', left_on='spHitX', right_index=True, copy=False)
    #initDF['spHitP'] = initDF['swissprot_blastp'].str.extract('sp\|(\S*)\|')
    #initDF = initDF.merge(spKO, how=left, left_on='spHitP', right_index=True, copy=False)
    return (initDF)

def addGO(initDF, sp2goentrez):
    sp2goD = pd.read_table(sp2goentrez, header=None, index_col = 0,usecols=[0,6], names=['spHitX','GO'])
    #names=['UniProtKB-ID','GeneID(EntrezGene)','RefSeq','GI','PDB','GO','UniRef100','UniRef90','UniRef50','UniParc','PIR','NCBI-taxon','MIM','UniGene','PubMed','EMBL','EMBL-CDS','Ensembl','Ensembl_TRS', 'Ensembl_PRO','Additional PubMed'])
    initDF = initDF.merge(sp2goD, how='left', left_on='spHitX', right_index=True, copy=False)
    return(initDF)


def get_eggNOG(initDF, sp2nog, nog2function):
    sp2nogD = pd.read_table(sp2nog, header=None, index_col = 0, names=['db', 'eggNOG'])
    sp2nogD.drop('db', axis=1, inplace=True)
    nog2funcD = pd.read_table(nog2function, header=None, names= ['db','eggNOG', 'num', 'num2','eggNOG_function','info'])
    nog2funcD = nog2funcD.loc[:,['eggNOG', 'eggNOG_function']]
    nog2funcD.set_index('eggNOG', inplace=True)
    sp2nogD = sp2nogD.merge(nog2funcD, how='left', left_on='eggNOG', right_index=True, copy=False)
    sp2nogD.fillna('', inplace=True)
    sp2nogD.reset_index(level=0, inplace=True)
    groupNOG = sp2nogD.groupby('index')[['eggNOG', 'eggNOG_function']].agg(lambda x: ''.join(x)).reset_index()
    # this results in some ",S" columns, and likely some " , " columns --> an issue.
    #sp2nogD.groupby('index').apply(lambda x: ','.join(x)).reset_index()  #[['eggNOG', 'eggNOG_function']].head()
#    groupNOG = sp2nogD.groupby(3)[1].apply(lambda x: ','.join(x)).reset_index()
    #groupNOG.rename(columns={3:'spHitX', 1:'eggNOG'}, inplace=True)
    groupNOG.rename(columns={'index':'spHitX'}, inplace=True)
    groupNOG.set_index('spHitX', inplace=True)
    initDF = initDF.merge(groupNOG, how='left', left_on='spHitX', right_index=True, copy=False)
    return(initDF)


if __name__ == '__main__':
    psr = argparse.ArgumentParser(description='combine annotations into a table')
    psr.add_argument('--fasta',help='fasta assembly file')
    psr.add_argument('--geneTransMap',help='tab separated gene to transcript conversion file')
    psr.add_argument('--spX',help='swissprot blastx results (long form, after extension to include full hit length and name)')
    psr.add_argument('--ur90X',help='uniref90 blastx results (long form, after extension to include full hit length and name)')
    psr.add_argument('--transdecoder_bed',help='transdecoder BED format output of ORFs')
    psr.add_argument('--spP',help='swissprot blastp results (long form, after extension to include full hit length and name)')
    psr.add_argument('--ur90P',help='uniref90 blastp results (long form, after extension to include full hit length and name)')
    psr.add_argument('--pfam',help='PFAM domain hits resulting from HMMER mapping of ORFs to pfam-A database')
    psr.add_argument('--signalP',help='signal peptide predictions from signalp (short form tabular output)')
    psr.add_argument('--tmhmm',help='transmembrane domain predictions from tmhmm (tabular output)')
    ###
    psr.add_argument('--nrX',help='NR blastx results (long form, after extension to include full hit length and name)')
    psr.add_argument('--nrP',help='NR blastp results (long form, after extension to include full hit length and name)')
    #psr.add_argument('--closest_blastx',metavar='contig2closest',help='what is this')
    #conversion files
    #IDMAPPING
    psr.add_argument('--sp2ko',help='swissprot to kegg orthology conversion')
    psr.add_argument('--sp2nog',help='swissprot id to eggNOG orthology conversion')
    psr.add_argument('--sp2ortho',help='swissprot to orthoDB')
    psr.add_argument('--sp2bioc',help='swissprot to BioCyc')
    #IDMAPPING SELECTED
    psr.add_argument('--sp2goentrez',help='swissprot to GO Entrez')
    #OTHER
    psr.add_argument('--ko2path',help='Kegg Orthology to Kegg Pathways Mapping File')
    psr.add_argument('--nog2function',help='eggNOG orthology to eggNOG functional classification')
    #NOT IN USE RIGHT NOW
    psr.add_argument('--sp2enzyme',help='Swissprot id to kegg enzyme orthology')
    psr.add_argument('--enzyme2path',help='kegg enzyme to kegg pathway conversion file')
    psr.add_argument('--pfam2enzyme',help='pfam domain to kegg enzyme conversion')
    psr.add_argument('--go2path',help='GO term to kegg pathway conversion')
    psr.add_argument('--go2slim',help='GO term to GO-slim term conversion')

    #outfile name
    psr.add_argument('--outfile',metavar='outfile',help='output filename')

    args = psr.parse_args()
    main(args)
