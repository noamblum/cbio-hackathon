#!/usr/bin/env python
import argparse
import pyBigWig

parser = argparse.ArgumentParser("Convert BigWig file to BED file")
parser.add_argument("input_file", help="BigWig file to convert", type=str)
parser.add_argument("output_file", help="Output file path", type=str)
parser.add_argument("-t", "--threshold", help="Threshold for conversion", type=int, default=10, required=False)
args = parser.parse_args()

threshold = args.threshold
bw_path = args.input_file
of_path = args.output_file

with pyBigWig.open(bw_path) as bw, open(of_path, "w") as of:
    for chrom, len in bw.chroms().items():
        intervals = bw.intervals(chrom)
        for interval in intervals:
            if abs(interval[2]) > threshold:
                of.write("{}\t{}\t{}\n".format(chrom, interval[0], interval[1]))
