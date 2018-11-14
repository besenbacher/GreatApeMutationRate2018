#!/usr/bin/env python
import datapaths
import sys
from bx.seq.twobit import *
import string

#six different mutations:
#AT -> CG,CG,TA
#CG -> GC,TA,AT

radius = 1
genome = sys.argv[1]
t = TwoBitFile(open(datapaths.TwoBit(genome)))

complement = string.maketrans('ATCGN', 'TAGCN')
def reverse_complement(s):
    return s.translate(complement)[::-1]

for line in sys.stdin:
    L = line.split()
    chrom = L[0]
    start = int(L[1]) - 1
    context = t[chrom][start-radius:start+radius+1]

    if len(context) < 2*radius+1:
        continue

    if context[radius]==context[radius].upper():
        repeat = 'F'
    else:
        repeat = 'T'

    context = context.upper()

    if context[radius] in ['T','G']:
       context = reverse_complement(context)

    if context[radius] == 'A':
        context_type = 'Weak'
    elif context[radius] == 'C':
        if context[radius+1] == 'G':
            context_type = 'CpG'
        else:
            context_type = 'nonCpG-Strong'
    else:
        continue
    print '\t'.join(L) + '\t' + context_type + '\t' + repeat
