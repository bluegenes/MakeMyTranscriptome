import os
import glob
import argparse
import parse_qual_metrics as parsers
import json


relative_paths = {'transrate': 'quality_files/transrate/',
                  'busco': 'assembly_files',
                  'fastqc': 'assembly_files',
                  'cegma': 'quality_files',
                  'history': 'log_files'}

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
        return parsers.get_history(files[0])
    else:
        return {}


def get_data(assembly_dir):
    data = {}
    data['busco'] = get_busco_info(assembly_dir)
    data['transrate'] = get_transrate_info(assembly_dir)
    data.update(get_fastqc_data(assembly_dir))
    data['cegma'] = get_cegma_info(assembly_dir)
    data['task_info'] = get_history_info(assembly_dir)
    return data


def create_report(assembly_dir, json_target, human_target):
    data = get_data(assembly_dir)
    f = open(json_target, 'w')
    json.dump(data, f, sort_keys=True, indent=4)
    f.close()


if(__name__ == '__main__'):
    parser = argparse.ArgumentParser(description=(
        "Description: A script that will parse an assembly directory "
        "and output a report compiling some of the data."))
    parser.add_argument('-a', '--assembly_dir', help='The assembly directory.')
    parser.add_argument('-o', '--out', help='The output name for the report')
    parser.add_argument('-j', '--json', help=(
        'The location that the json encoded version of the report will be printed.'))
    parser.add_argument('-q', '--qualityDir', default='quality_files', help=(
        'optional: alternative name for quality directory'))
    parser.add_argument('-t', '--transrateDir', default='transrate', help=(
        'optional: alternative name for transrate directory',))
    parser.add_argument('--transratePostAssemblyDir', default='transrate_post_assembly', help=(
        'optional: alternative name for transrate directory'))
    parser.add_argument('--assembly_filesDir', default='assembly_files', help=(
        'optional: alternative name for transrate directory'))
    args = parser.parse_args()
    main(args.assembly_dir, args.json, args.out)
