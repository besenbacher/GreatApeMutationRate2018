#!/usr/bin/env python
import datapaths
import sys
from bx.seq.twobit import *
import string
import argparse
parser = argparse.ArgumentParser(description='Check_reference')
parser.add_argument('two_bit_file', type=str,
                    help='reference genome in 2bit format')
args = parser.parse_args()

t = TwoBitFile(open(args.two_bit_file, 'rb'))

for line in sys.stdin:
    L = line.split()
    chrom = L[0]
    if chrom == "CHROM":
        print '\t'.join(L) + '\t' + "repeat"
        continue
    start = int(L[1]) - 1
    ref = t[chrom][start:start+1]
    assert ref.upper()==L[2]
    if ref == ref.upper():
        repeat = "F"
    else:
        repeat = "T"
    print '\t'.join(L) + '\t' + repeat
