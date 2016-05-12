import os
import glob
import argparse
import parse_qual_metrics as parsers
import json
import time
import pandas as pd

relative_paths = {'transrate': 'quality_files/transrate/',
                  'busco': 'quality_files',
                  'fastqc': 'assembly_files',
                  'cegma': 'quality_files',
                  'history': 'log_files',
                  'summary': 'annotation_files'}

buscoDirBase = "run_busco_*"
buscoShortSum = 'short_summary*'
fastqc_pre_trim = "fastqc_pre*"
fastqc_post_trim = "fastqc_post*"
fastqc_final = "fastqc_final*"

cegmaFile = "*.completeness_report"
buscoFile = "short_summary*"
transrateFile = "*assemblies.csv"
transrateReadCountF = "*read_count.txt"
detonateFile = "*.score"
fastqcFile = "fastqc_data.txt"
summaryFile = 'annotation_summary.json'


def get_busco_info(assembly_dir):
    pattern = os.path.join(assembly_dir, relative_paths['busco'], buscoDirBase, buscoShortSum)
    busco_summaries = glob.glob(pattern)
    ret = {}
    for p in busco_summaries:
        ret.update(parsers.get_busco_info(p))
    return ret


def get_fastqc_data(assembly_dir):
    ret = {}
    fastqc_dir_pattern = os.path.join(assembly_dir, relative_paths['fastqc'], fastqc_pre_trim)
    ret['fastqc_pre_trimming'] = parsers.fastqc_parser(fastqc_dir_pattern, 'fastqc_data.txt')
    fastqc_dir_pattern = os.path.join(assembly_dir, relative_paths['fastqc'], fastqc_post_trim)
    ret['fastqc_post_trimming'] = parsers.fastqc_parser(fastqc_dir_pattern, 'fastqc_data.txt')
    fastqc_dir_pattern = os.path.join(assembly_dir, relative_paths['fastqc'], fastqc_final)
    ret['fastqc_final'] = parsers.fastqc_parser(fastqc_dir_pattern, 'fastqc_data.txt')
    return ret


def get_transrate_info(assembly_dir):
    trans_dir = os.path.join(assembly_dir, relative_paths['transrate'])
    ret = parsers.transrate_parser(trans_dir, transrateFile)
    return ret


def get_cegma_info(assembly_dir):
    return parsers.cegma_parser(os.path.join(assembly_dir, relative_paths['cegma']), cegmaFile)


def get_history_info(assembly_dir):
    files = glob.glob(os.path.join(assembly_dir, relative_paths['history'], 'history.json'))
    if(len(files) > 0):
        hist_data = parsers.get_history(files[0])
        for t in hist_data:
            hist_data[t]['date'] = time.strftime("%a %b %d %H:%M:%S %Y", hist_data[t]['date'])
        return hist_data
    else:
        return {}


def get_summary_info(assembly_dir):
    files = glob.glob(os.path.join(assembly_dir, relative_paths['summary'], summaryFile))
    if(len(files) > 0):
        f = open(files[0])
        data = json.load(f)
        f.close()
        return data


def get_data(assembly_dir):
    data = {}
    data['busco'] = get_busco_info(assembly_dir)
    data['transrate'] = get_transrate_info(assembly_dir)
    data.update(get_fastqc_data(assembly_dir))
    data['cegma'] = get_cegma_info(assembly_dir)
    data['task_info'] = get_history_info(assembly_dir)
    data['annot_summary'] = get_summary_info(assembly_dir)
    return data


def create_truncated_csv(dataD):
    s1  = pd.Series({'fastqc_pre_trim': dataD['fastqc_pre_trimming'].items()})
    s2  = pd.Series({'fastqc_post_trim': dataD['fastqc_post_trimming'].items()})
    s3  = pd.Series({'fastqc_final': dataD['fastqc_final'].items()})
    s4 = pd.Series(dataD['busco'],name= 'busco')
    s5 = pd.Series(dataD['cegma'],name= 'cegma')
    s6 = pd.Series(dataD['transrate'],name= 'transrate')
    s7 = pd.Series(dataD['annot_summary'],name= 'annot_summary')
    #s8 = pd.Series(dataD['task_info'],name= 'task_info')
    #dataDF = pd.concat([s1,s2,s3,s4,s5,s6,s7], axis=0)
    dataDF = pd.concat([s1,s2,s3,s4,s5,s6], axis=0)
    return dataDF
#    s1  = pd.Series({'fastqc_pre_trim']: [x]+y for (x,y) in dataD['fastqc_pre_trimming'].items()}, name='fastqc_pre_trimming')
#    s2  = pd.Series({'fastqc_post_trim_'+ x: y for (x,y) in dataD['fastqc_post_trimming'].items()}, name='fastqc_post_trimming')
#    s3  = pd.Series({'fastqc_final_'+ x: y for (x,y) in dataD['fastqc_final'].items()}, name='fastqc_final')


def create_report(assembly_dir, json_target=None, csv_target=None, human_target=None):
    data = get_data(assembly_dir)
    if(json_target is None):
        json_target = os.path.join(assembly_dir, 'log_files', 'report.json')
    if(csv_target is None):
        csv_target = os.path.join(assembly_dir, 'log_files', 'report.csv')
    f = open(json_target, 'w')
    json.dump(data, f, sort_keys=True, indent=4)
    f.close()
    f = open(csv_target, 'w') 
    dataF = create_truncated_csv(data)
    dataF.to_csv(csv_target, sep=',')
    f.close()


if(__name__ == '__main__'):
    parser = argparse.ArgumentParser(description=(
        "Description: A script that will parse an assembly directory "
        "and output a report compiling some of the data."))
    parser.add_argument('-a', '--assembly_dir', help='The assembly directory.')
    parser.add_argument('-o', '--out', help='The output name for the report')
    parser.add_argument('-j', '--json', help=(
        'The location that the json encoded version of the report will be printed.'))
    parser.add_argument('-c', '--csv', help=(
        'The location that the csv encoded version of the report will be printed.'))
    parser.add_argument('-q', '--qualityDir', default='quality_files', help=(
        'optional: alternative name for quality directory'))
    parser.add_argument('-t', '--transrateDir', default='transrate', help=(
        'optional: alternative name for transrate directory'))
    parser.add_argument('--transratePostAssemblyDir', default='transrate_post_assembly', help=(
        'optional: alternative name for transrate directory'))
    parser.add_argument('--assembly_filesDir', default='assembly_files', help=(
        'optional: alternative name for transrate directory'))
    args = parser.parse_args()
    create_report(args.assembly_dir, args.json, args.csv, args.out)
