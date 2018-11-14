#!/usr/bin/python

D = {}
import sys

for line in sys.stdin:
    line = line.strip()
    if line not in D:
        D[line] = 0
    D[line]+=1

for line in D:
    print line, D[line]
