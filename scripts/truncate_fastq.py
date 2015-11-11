import argparse
import random
import sys


def fq_parser(fastq):
    f = open(fastq)
    ret = []
    count = 0
    for line in f:
        ret.append(line.rstrip('\n'))
        count += 1
        if(count % 4 == 0):
            yield ret
            ret = []
    f.close()


def fq_writer(f, entry):
    for line in entry:
        f.write(line)
        f.write('\n')


def main(fastq, target=None, length=50, cut_type='left'):
    outfile = sys.stdout if(target is None) else open(target, 'w')
    for entry in fq_parser(fastq):
        seq = entry[1]
        scores = entry[3]
        if(cut_type == 'left'):
            seq = seq[:length]
            scores = scores[:length]
        elif(cut_type == 'right'):
            seq = seq[-1*length:]
            scores = scores[-1*length:]
        elif(cut_type == 'random' and len(seq) > length):
            offset = random.randint(0, len(seq)-length)
            seq = seq[offset:offset+length]
            scores = scores[offset:offset+length]
        entry[1] = seq
        entry[3] = scores
        fq_writer(outfile, entry)
    outfile.close()


if(__name__ == '__main__'):
    parser = argparse.ArgumentParser(description=(
        'Generates a virtual set of reads based on existing fastq files by '
        'truncating reads. Default behavoir is to keep the first --length bases '
        'of each read. Flags can be used to capture other sections instead.'))
    parser.add_argument(
        'fastq', help='The path to a fastq file containing reads that need to be truncated.')
    parser.add_argument('-l', '--length', default=50, type=int,
                        help='The length to be truncated to. Default is 50.')
    parser.add_argument('-t', '--target', help='The output file. Default is stdout.')
    parser.add_argument('-random', action='store_true', help=(
        'Use this flag to signify that truncation should occur via randomly '
        'chosen windows rather than left or right truncations'))
    parser.add_argument('-right', action='store_true', help=(
        'Use this flag to signify that the last --length bases of each read '
        'should be kept instead.'))
    args = parser.parse_args()
    cut_type = 'left'
    if(args.random):
        cut_type = 'random'
    elif(args.right):
        cut_type = 'right'
    main(args.fastq, args.target, args.length, cut_type)


