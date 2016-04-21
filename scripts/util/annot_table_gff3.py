'''
TODO - 
    check to make sure protein -> rna indexing conversions are correct
    add support for signalP and tmhmm
    add support for secondary files as appropriate
    find more efficient way to join sorted DataFrames
'''


import argparse
import pandas as pd
import csv
import os
from itertools import count
import re

gff_colnames = ['seqid', 'source', 'type', 'start', 'end',
                'score', 'strand', 'phase', 'attributes']

blast_colnames = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                  'gapopen', 'qstart', 'qend', 'sstart', 'send',
                  'evalue', 'bitscore']

blast_ex_colnames = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                     'gapopen', 'qstart', 'qend', 'sstart', 'send',
                     'evalue', 'bitscore', 'stitle', 'slen']

ID_GEN = count()


def fasta_to_gff3(fasta):
    with open(fasta, 'r') as f:
        lens = []
        keys = []
        length = 0
        transcriptid = None
        for line in f:
            if line.startswith('>'):
                if length > 0:
                    lens.append(length)
                    keys.append(transcriptid)
                    length = 0
                transcriptid = line.split(' ')[0][1:]
            else:
                length += len(line.strip())
        lens.append(length)
        keys.append(transcriptid)
    retDF = pd.DataFrame(columns=gff_colnames)
    retDF['seqid'] = keys
    retDF['source'] = ['Trinity'] * len(keys)  # should allow rnaspades as well
    retDF['type'] = ['mRNA'] * len(keys)
    retDF['start'] = ['1'] * len(keys)
    retDF['end'] = lens
    retDF['attributes'] = ['ID:{0!s}'.format(next(ID_GEN)) for k in range(len(keys))]
    return retDF


def orf_to_gff3(orf_gff3):
    retDF = pd.read_table(orf_gff3, comment='#', names=gff_colnames)
    def remap_id(row):
        attr = row.attributes
        attr = attr.split(';')
        attr = [x if(not re.match('ID=', x)) else 'ID={0!s}'.format(next(ID_GEN)) for x in attr]
        return ';'.join(attr)
    retDF['attributes'] = retDF.apply(remap_id, axis=1)
    return retDF


def blast_to_gff3(blast_file, blast_type, database=''):
    assert blast_type in ['BLASTX', 'BLASTP']
    ftype = 'protein_match' if(blast_type == 'BLASTP') else 'translated_nucleotide_match'
    blastDF = pd.read_table(blast_file, names=blast_ex_colnames, header=0)
    blastDF = blastDF.drop_duplicates(subset=['qseqid'], keep='first')
    retDF = pd.DataFrame()
    if(blast_type == 'BLASTP'):
        mod_fields = ['qstart', 'qend', 'sstart', 'send']
        blastDF[mod_fields] *= 3
        blastDF[mod_fields] -= 2
        retDF['seqid'] = blastDF.apply(lambda row: row.qseqid.split('|')[0], axis=1)
    else:
        retDF['seqid'] = blastDF['qseqid']
    retDF['source'] = [blast_type] * len(blastDF)
    retDF['type'] = [ftype] * len(blastDF)
    retDF['start'] = blastDF['qstart']
    retDF['end'] = blastDF['qend']
    retDF['score'] = blastDF['evalue']
    retDF['strand'] = ['.'] * len(blastDF)
    retDF['phase'] = ['.'] * len(blastDF)

    def blast_build_attr(row):
        data = []
        data.append('ID={0!s}:{1!s}'.format(next(ID_GEN), row.qseqid))
        data.append('Target={0!s} {1!s} {2!s}'.format(row.sseqid, row.sstart, row.send))
        data.append('Target_Title={0!s}'.format(row.stitle))
        if(database != ''):
            data.append('Database={0}'.format(database))
        return ';'.join(data)

    retDF['attributes'] = blastDF.apply(blast_build_attr, axis=1)
    return retDF


def hmmscan_row_to_list(row):
    row = row[0]
    fields = row.split()
    num_fields = 22
    ret = fields[:num_fields]
    ret.append(' '.join(fields[num_fields:]))
    return pd.Series(ret)


def hmmscan_to_gff3(hmmscan, database='PFAM'):
    hmmscanDF = pd.read_table(hmmscan, header=None, comment='#')
    hmmscanDF = hmmscanDF.apply(hmmscan_row_to_list, axis=1)
    num_cols = hmmscanDF.shape[1]
    mod_fields = [num_cols-3, num_cols-4]
    hmmscanDF = hmmscanDF.convert_objects(convert_numeric=True)
    hmmscanDF[mod_fields] *= 3
    hmmscanDF[mod_fields] -= 2
    retDF = pd.DataFrame()
    retDF['seqid'] = hmmscanDF.apply(lambda row: row[3].split('|')[0], axis=1)
    retDF['source'] = ['HMMER'] * len(hmmscanDF)
    retDF['type'] = ['protein_hmm_match'] * len(hmmscanDF)
    retDF['start'] = hmmscanDF[num_cols-4]
    retDF['end'] = hmmscanDF[num_cols-3]
    retDF['score'] = hmmscanDF[11]
    retDF['strand'] = ['.'] * len(hmmscanDF)
    retDF['phase'] = ['.'] * len(hmmscanDF)

    def build_attr(row):
        data = []
        data.append('ID={0}:{1}'.format(next(ID_GEN), row[3]))
        data.append('Name={0}'.format(row[0]))
        data.append('Target={0} {1} {2} +'.format(row[0], row[15], row[16]))
        data.append('Note={0}'.format(row[num_cols-1]))
        data.append('accuracy={0}'.format(row[num_cols-2]))
        if database:
            data.append('Dbxref="{0}:{1}"'.format(database, row[1]))
        return ';'.join(data)

    retDF['attributes'] = hmmscanDF.apply(build_attr, axis=1)

    return retDF


def signalp_to_gff3(infiles):
    pass


def tmhmm_to_gff3(infiles):
    pass


def file_check(f):
    return f is not None and os.path.exists(f)


def main(args):
    gff3_entries = []
    gff3_entries.append(pd.DataFrame(columns=gff_colnames))
    if(file_check(args.fasta)):
        gff3_entries.append(fasta_to_gff3(args.fasta))
    if(file_check(args.transdecoder_gff3)):
        gff3_entries.append(orf_to_gff3(args.transdecoder_gff3))
    if(file_check(args.spX)):
        gff3_entries.append(blast_to_gff3(args.spX, 'BLASTX', 'uniprot_sprot'))
    if(file_check(args.spP)):
        gff3_entries.append(blast_to_gff3(args.spP, 'BLASTP', 'uniprot_sprot'))
    if(file_check(args.ur90X)):
        gff3_entries.append(blast_to_gff3(args.ur90X, 'BLASTX', 'uniref90'))
    if(file_check(args.ur90P)):
        gff3_entries.append(blast_to_gff3(args.ur90P, 'BLASTP', 'uniprot90'))
    if(file_check(args.nrX)):
        gff3_entries.append(blast_to_gff3(args.nrX, 'BLASTX', 'nr'))
    if(file_check(args.nrP)):
        gff3_entries.append(blast_to_gff3(args.nrP, 'BLASTP', 'nr'))
    if(file_check(args.pfam)):
        gff3_entries.append(hmmscan_to_gff3(args.pfam))
    df = pd.concat(gff3_entries)
    df = df.sort('seqid')
    df.to_csv(args.outfile, sep='\t', na_rep='.', index=False,
              columns=gff_colnames, header=False,
              quoting=csv.QUOTE_NONE)
    return df


if __name__ == '__main__':
    psr = argparse.ArgumentParser(description=(
        'combine annotations into a table'))
    psr.add_argument('--fasta', help='fasta assembly file')
    psr.add_argument('--geneTransMap', help=(
        'tab separated gene to transcript conversion file'))
    psr.add_argument('--spX', help=(
        'swissprot blastx results (long form, after extension to include '
        'full hit length and name)'))
    psr.add_argument('--ur90X', help=(
        'uniref90 blastx results (long form, after extension to include '
        'full hit length and name)'))
    psr.add_argument('--nrX', help=(
        'NR blastx results (long form,  after extension to include full hit '
        'length and name)'))
    psr.add_argument('--transdecoder_gff3', help=(
        'transdecoder gff3 format output of ORFs'))
    psr.add_argument('--spP', help=(
        'swissprot blastp results (long form, after extension to include full '
        'hit length and name)'))
    psr.add_argument('--ur90P', help=(
        'uniref90 blastp results (long form, after extension to include full '
        'hit length and name)'))
    psr.add_argument('--nrP', help=(
        'NR blastp results (long form,  after extension to include full hit '
        'length and name)'))
    psr.add_argument('--pfam', help=(
        'PFAM domain hits resulting from HMMER mapping of ORFs to pfam-A '
        'database'))
    psr.add_argument('--signalP', help=(
        'signal peptide predictions from signalp (short form tabular output)'))
    psr.add_argument('--tmhmm', help=(
        'transmembrane domain predictions from tmhmm (tabular output)'))
    # outfile name
    psr.add_argument('--outfile', metavar='outfile', help='output filename')

    args = psr.parse_args()
    main(args)
