#!/usr/bin/env python
import argparse
import pyBigWig
from numpy import average

parser = argparse.ArgumentParser("Convert BigWig file to BED file")
parser.add_argument("input_file", help="BigWig file to convert", type=str)
parser.add_argument("output_file", help="Output file path", type=str)
parser.add_argument("-t", "--threshold", help="Threshold for conversion", type=float, default=10, required=False)
args = parser.parse_args()

threshold = args.threshold
bw_path = args.input_file
of_path = args.output_file

id = 1

with pyBigWig.open(bw_path) as bw, open(of_path, "w") as of:
    for chrom, len in bw.chroms().items():
        intervals = bw.intervals(chrom)
        start = None
        end = None
        vals = []
        weights = []
        for interval in intervals:
            if abs(interval[2]) > threshold:
                if start is None:
                    start = interval[0]
                    end = interval[1]
                    vals.append(abs(interval[2]))
                    weights.append(interval[1] - interval[0])
                elif interval[0] == end:
                    end = interval[1]
                    vals.append(abs(interval[2]))
                    weights.append(interval[1] - interval[0])
                else:
                    sign = "+" if interval[2] > 0 else "-"
                    of.write(f"{chrom}\t{start}\t{end}\t{id}\t{average(vals, weights=weights)}\t{sign}\n")
                    start = interval[0]
                    end = interval[1]
                    vals = [abs(interval[2])]
                    weights = [interval[1] - interval[0]]
                    id += 1
