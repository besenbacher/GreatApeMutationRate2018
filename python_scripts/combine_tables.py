#!/bin/env python

import sys

D = {}
for line in sys.stdin:
    L = line.split()
    count = int(L[-1])
    key = tuple(L[:-1])
    if key not in D:
        D[key] = 0
    D[key] += count

for key in D:
    print ' '.join(key), D[key]

    
