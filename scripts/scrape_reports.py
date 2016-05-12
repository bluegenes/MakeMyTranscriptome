import argparse
import json
import os
import sys
import pandas as pd
import assembly_report as ar

assembly_dir_identifiers = ['log_files', 'assembly_files', 'annotation_files', 'expression_files']
relative_path_report = 'log_files/report'


def _scrape_assembly(root, depth=1, link_flag=False, force_report=False):
    if(depth < 0):
        return {}
    dirs = set(os.listdir(root))
    # ask is assembly directory
    if(all(d in dirs for d in assembly_dir_identifiers)):
        fj = os.path.join(root, relative_path_report + '.json')
        fc = os.path.join(root, relative_path_report + '.csv')
        # if report doesn't exist, create it. If 'force', recreate it.
        if (not os.path.isfile(fj) or force_report):
	    ar.create_report(root, json_target=fj, csv_target=fc)
	    print(fc)
        with open(fj) as rfj:
	    with open(fc) as rfc:
                ret = {}
		df = {}
                ret[os.path.basename(root)] = json.load(rfj)
		try:
		    df[os.path.basename(root)] = pd.Series(pd.read_csv(rfc, index_col=0, header=None,squeeze=True), name=os.path.basename(root))
		    #import pdb;pdb.set_trace()
		except:
                    df = {}
		    #import pdb;pdb.set_trace()
        return ret, df
    else:
        dirs = [os.path.join(root, d) for d in dirs]
        dirs = [d for d in dirs if(not os.path.islink(d) or link_flag)]
        dirs = [d for d in dirs if(os.path.isdir(d))]
        ret = {}
	df = {}
        for d in dirs:
	    info, dfInfo = _scrape_assembly(d, depth-1, link_flag, force_report)
            ret.update(info)
            df.update(dfInfo)
        return ret, df


def scrape_assembly_reports(roots, depth=1, link_flag=False, force_report=False):
    roots = [os.path.abspath(r) for r in roots]
    ret = {}
    df =  {}
    for r in roots:
        json_info, csv_info =_scrape_assembly(r, depth, link_flag, force_report)
	ret.update(json_info) #update json
	df.update(csv_info)
    return ret,df


if(__name__ == '__main__'):
    parser = argparse.ArgumentParser(description=(
        "Descritpion: A script that will search for and collect information "
        "from MMT assemblies."))
    parser.add_argument('-r', '--roots', default='.', help=(
        'a list of comma seperated paths to directories that should be '
        'searched. Default is current directory'))
    parser.add_argument('-d', '--depth', default=1, type=float, help=(
        'the dpeth to which the search will be carried out. Default is 0. Use '
        '-1 or "inf" to specify a complete search.'))
    parser.add_argument('-l', '--links', action='store_true', help=(
        'use this flag to signify that links should be followed during '
        'recursion.'))
    parser.add_argument('-f', '--force', action='store_true', help=(
        'force creating of a new assembly report'))
    parser.add_argument('-o', '--output', help=(
        'The location to write all json data to. Default is stdout.'))
    args = parser.parse_args()
    if(args.depth < 0):
        args.depth = float('inf')
    args.roots = args.roots.split(',')
    json_data, df_data = scrape_assembly_reports(args.roots, args.depth, args.links, args.force)
    out = sys.stdout if(args.output is None) else open(args.output, 'w')
    json.dump(json_data, out, sort_keys=True, indent=4)
    out_csv = args.output.rsplit('.')[0] + '.csv'
    #import pdb;pdb.set_trace()
    #infoDF = pd.DataFrame.from_dict(df_data)
    infoDF = pd.DataFrame(df_data).transpose()
    #import pdb;pdb.set_trace()
    infoDF.to_csv(out_csv, sep=',') #, sort_keys=True, indent=4)
    #print('')
    out.close()
