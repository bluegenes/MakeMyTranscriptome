import sys
import re
import os
import argparse
import json
from itertools import chain,combinations

edger_keys = ['contig_id','logFC','logCPM','PValue','FDR']
deseq_keys = ['contig_id','baseMean','logFC','lfcSE','stat','pvalue','FDR']
key_maps = [edger_keys,deseq_keys]


def de_file_parser(f,key_map):
	rf = open(f)
	rf.readline()
	num_fields = len(key_map)
	for line in rf:
		line = line.strip().split()
		yield_val = {key_map[i]:line[i] for i in range(num_fields)}
		yield yield_val
	rf.close()


def handle_de_output(counts_table,de_files,threshold,key_map,outfile):
	ct_file = open(counts_table)
	ct_header = ct_file.readline().strip().split('\t')[1:]
	countsD = {entries[0]:entries[2:] for entries in (line.strip().split("\t") for line in ct_file)}	
	ct_file.close()
	de_dicts = {f:{entry['contig_id']:entry for entry in de_file_parser(f,key_map)} for f in de_files}
	de_counts = {key:0 for key in chain(*[combinations(de_files,x) for x in range(len(de_files)+1)])}
	contig_ids = set(key for key in chain(*[de_dicts[f] for f in de_files]))
	outfile = open(outfile,'w')
	header = ['contig']
	for f in de_files:
		header.extend([os.path.basename(f)+'_'+key for key in key_map if key!='contig_id'])
	header.extend(ct_header)
	outfile.write('\t'.join(header)+'\n')
	for contig in contig_ids:
		sig_set = []
		for f in de_files:
			if(contig not in de_dicts[f] or de_dicts[f][contig]['FDR']=='NA' ):
				continue
			if(float(de_dicts[f][contig]['FDR'])<threshold):
				sig_set.append(f)
		for key in chain(*[combinations(sig_set,x) for x in range(len(sig_set)+1)]):
			de_counts[key]+=1
		if(len(sig_set)>0):
			temp_strs = [contig]
			for f in de_dicts:
				if(contig not in de_dicts[f]):
					temp_strs.extend(['.' for key in key_map if key!='contig_id'])
				else:
					temp_strs.extend([de_dicts[f][contig][key] for key in key_map if key!='contig_id'])
			temp_strs.extend(countsD[contig])
			outfile.write('\t'.join(temp_strs)+'\n')
	print(json.dumps([{'count':de_counts[k],'set':k} for k in sorted(de_counts)],indent=4))


if(__name__=='__main__'):
	parser = argparse.ArgumentParser(description="handles de_output files")
	parser.add_argument('-i','--input',help='a comma seperated list of de_files')
	parser.add_argument('-k','--key',help='use to specify defile type. "0" for edger and "1" for deseq.',type=int)
	parser.add_argument('-c','--counts_table',help='The counts_to_table file for the de_files')
	parser.add_argument('-t','--threshold',help='FDR threshold value',type=float)
	parser.add_argument('-o','--outfile',help='The file to write significant hits to')
	args = parser.parse_args()
	args.input=args.input.split(',')
	args.key = key_maps[args.key]
	handle_de_output(args.counts_table,args.input,args.threshold,args.key,args.outfile)



