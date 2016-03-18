import argparse
import json
import os
import sys

assembly_dir_identifiers = ['log_files', 'assembly_files', 'annotation_files', 'expression_files']
relative_path_report = 'log_files/report.json'


def _scrape_assembly(root, depth=1, link_flag=False):
    if(depth < 0):
        return {}
    dirs = set(os.listdir(root))
    # ask is assembly directory
    if(all(d in dirs for d in assembly_dir_identifiers)):
        f = os.path.join(root, relative_path_report)
        # ask if report already exists.
        if(os.path.isfile(f)):
            f = open(f)
            ret = {}
            ret[os.path.basename(root)] = json.load(f)
            f.close()
            return ret
        else:
            return {}
    else:
        dirs = [os.path.join(root, d) for d in dirs]
        dirs = [d for d in dirs if(not os.path.islink(d) or link_flag)]
        dirs = [d for d in dirs if(os.path.isdir(d))]
        ret = {}
        for d in dirs:
            ret.update(_scrape_assembly(d, depth-1, link_flag))
        return ret


def scrape_assembly_reports(roots, depth=1, link_flag=False):
    roots = [os.path.abspath(r) for r in roots]
    ret = {}
    for r in roots:
        ret.update(_scrape_assembly(r, depth, link_flag))
    return ret


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
    parser.add_argument('-o', '--output', help=(
        'The location to write all json data to. Default is stdout.'))
    args = parser.parse_args()
    if(args.depth < 0):
        args.depth = float('inf')
    args.roots = args.roots.split(',')
    data = scrape_assembly_reports(args.roots, args.depth, args.links)
    out = sys.stdout if(args.output is None) else open(args.output, 'w')
    json.dump(data, out, sort_keys=True, indent=4)
    print('')
    out.close()
