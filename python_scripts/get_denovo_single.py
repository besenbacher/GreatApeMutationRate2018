#!/bin/env python
import sys
import vcf
import random
import argparse
from bx.seq.twobit import *
import string
from dr_utils_ver2 import *
from collections import Counter

parser = argparse.ArgumentParser(description='''
Find unfiltered denovo mutations and variants to test callability
''')

parser.add_argument('family_description', type=str,
                    help='family description file')
parser.add_argument('outdir', type=str,
                    help='output directory')
parser.add_argument('--two_bit', type=argparse.FileType('rb'),
                    help='reference genome in 2bit format')
parser.add_argument('child', type=str,
                    help='child')
parser.add_argument('--postfix', type=str, default = '',
                    help='add this post fix to created files')
parser.add_argument('--verbose', action='store_true', help="print progress output")
parser.add_argument('--var_type', type=str, default="SNV", help="should it look for SNV or indel variants")
parser.add_argument('--no_test_data', action='store_true',
                    help="program should not create het_test and homeref_test files")
parser.add_argument('--PoO_data', default=[], type=str, nargs='+',
                    help='parent of origin data files')
args = parser.parse_args()

families = read_families(args.family_description)
families = [x for x in families if x.child == args.child]
vcf_reader = vcf.Reader(sys.stdin)
outdir =  args.outdir
make_test_data = not args.no_test_data
make_indel = args.var_type == "indel"
radius = 1
tb = TwoBitFile(args.two_bit)

poo = {}
for poo_file in args.PoO_data:
    f = open(poo_file)
    for line in f:
        L = line.split()
        poo[tuple(L[0].split(';'))] = L[1:]
    f.close()

n_poo_vars = 0
for x in poo:
    n_poo_vars = len(poo[x])
    break

if n_poo_vars == 0:
    POO_header = []
else:
    POO_header = ['POO', 'POO_method', 'POO_LR', 'P_prob', 'M_prob', 'mtype', 'cqual', 'POO_nreads', 'POO_nj']

assert n_poo_vars == len(POO_header)

q_measures = ['BaseQRankSum', 'FS','MQ', 'MQRankSum', 'QD', 'ReadPosRankSum', 'SOR', 'AN', 'ClippingRankSum', 'ExcessHet']
list_q_measures = []
values = ['GT', 'AD1', 'AD2', 'DP', 'GQ', 'PL00', 'PL01', 'PL11', 'SAC1', 'SAC2', 'SAC3', 'SAC4']
min_other_GQ = 30
min_test_parent_GQ = 30
min_test_parent_DP = 10
min_test_lines = 5000
homref_DP_count = Counter()
het_DP_count = Counter()

ignore_chromosomes = ['chrX', 'chrY', 'chrM']

denovo_out = open(outdir + '/denovo_raw_' + args.postfix + '.dat', 'w')
if make_test_data:
    het_test_out = open(outdir + '/het_test_' + args.postfix + '.dat', "w")
    homoref_test_out = open(outdir + '/homoref_test_' + args.postfix + '.dat', 'w')

complement = string.maketrans('ATCGN', 'TAGCN')
def reverse_complement(s):
    return s.translate(complement)[::-1]

def get_context(chrom,pos):
    context = tb[chrom][pos-radius-1:pos+radius].upper()
    if context[radius] in ['T','G']:
        context = reverse_complement(context)
    return context

def max_altAD_other_families(x, family):
    maxAD = 0
    otherAC = 0
    sumAD = 0
    maxOtherGQ = 0
    fampns = set(x for x in family)
    for call in x.samples:
        if call.sample not in fampns:
            try:
                if call.data.PL[0] >= min_other_GQ:
                    otherAC += 1
                maxOtherGQ = max(maxOtherGQ, call.data.PL[0])
                if call.data.AD[1] > maxAD:
                    maxAD = call.data.AD[1]
                sumAD += call.data.AD[1]
            except:
                continue
    return otherAC, maxAD, sumAD, maxOtherGQ

def has_good_GQ(x, childname):
    for call in x.samples:
        if call.gt_type != 0 and call.sample != childname and call.data.GQ > min_test_parent_GQ:
            return True
    return False
            #return call.data.GQ

def print_line(outfile, x, f, denovo=True):
    good = True
    try:
        context = get_context(x.CHROM, x.POS)
    except:
        context = "NNN"
        print "context_error"

    line = [x.CHROM, x.POS, x.REF, x.ALT[0], f.child, int(x.is_transition), x.nucl_diversity, x.QUAL, replace_empty(x.FILTER), x.INFO["MLEAC"][0], context]
    for val in q_measures:
        if val not in x.INFO:
            line.append("NA")
        else:
            line.append(x.INFO[val])

        for val, n in list_q_measures:
            if val not in x.INFO:
                line.extend(["NA"]*n)
            else:
                line.extend(x.INFO[val])

    for individual in f.trio() :
        if x.genotype(individual).gt_type == None:
            good = False
            break
        try:
            line.append(x.genotype(individual).gt_type)
            if (x.genotype(individual).data.AD == None or
                len(x.genotype(individual).data.AD)!=2 or
                x.genotype(individual).data.DP == None or
                x.genotype(individual).data.GQ == None or
                x.genotype(individual).data.PL == None or
                x.genotype(individual).data.SAC == None or
                len(x.genotype(individual).data.PL)!= 3):
                good = False
                break
            line.append(x.genotype(individual).data.AD[0])
            line.append(x.genotype(individual).data.AD[1])
            line.append(x.genotype(individual).data.DP)
            line.append(x.genotype(individual).data.GQ)
            line.append(x.genotype(individual).data.PL[0])
            line.append(x.genotype(individual).data.PL[1])
            line.append(x.genotype(individual).data.PL[2])
            line.append(x.genotype(individual).data.SAC[0])
            line.append(x.genotype(individual).data.SAC[1])
            line.append(x.genotype(individual).data.SAC[2])
            line.append(x.genotype(individual).data.SAC[3])
        except AttributeError:
            good = False
            break
    gc_n = 0
    gc_ac = 0
    if denovo:
        for individual in f.grandchildren:
            if x.genotype(individual).gt_type == None:
                continue
            try:
                gc_ac += x.genotype(individual).gt_type
                gc_n +=1
            except AttributeError:
                continue
    if good:
        if denovo:
            line.extend(max_altAD_other_families(x, f))
            line.append(gc_n)
            line.append(gc_ac)
            if (f.child, x.CHROM, str(x.POS), x.REF, str(x.ALT[0])) in poo:
                line.extend(poo[(f.child, x.CHROM, str(x.POS), x.REF, str(x.ALT[0]))])
            else:
                #print (f.child, x.CHROM, str(x.POS), x.REF, str(x.ALT[0]))
                line.extend(['NA']*n_poo_vars)
        print >> outfile, '\t'.join(str(x) for x in line)


header = ['CHROM','POS','REF','ALT','CHILD','ts','div','QUAL','FILTER','AC','CONTEXT'] + \
    q_measures + \
    ['\t'.join([x[0] + "." +str(i+1) for i in range(x[1])]) for x in list_q_measures] + \
    [individual+ '.' + value for individual in ['FATHER','MOTHER','CHILD'] for value in values]

print >> denovo_out, '\t'.join(header + ['OTHER.AC','MAX.OTHER.AD','SUM.OTHER.AD','MAX.OTHER.PL00','GC.N','GC.AC'] + POO_header)
if make_test_data:
    print >> het_test_out, '\t'.join(header)
    print >> homoref_test_out, '\t'.join(header)

chrom = "NA"
pos = "NA"
for x in vcf_reader:
    #I only look at biallelic variants!
    if len(x.INFO["MLEAC"])!=1:
        continue
    if make_indel:
        if x.is_snp:
            continue
    else:
        if not x.is_snp:
            continue
    if x.CHROM in ignore_chromosomes:
        continue
    if chrom != x.CHROM:
        chrom = x.CHROM
        print chrom
    if pos != int(x.POS/100000):
        pos = int(x.POS/100000)
        print x.POS
    for f in families:
        good = True
        try:
            for individual in f.trio():
                if (x.genotype(individual).data.GQ == None):
                    good = False
                    break
        except AttributeError:
            continue
        if not good:
            continue
        if (#x.INFO["MLEAC"][0] <=5 and
            x.genotype(f.father).data.GQ >= 10 and
            x.genotype(f.mother).data.GQ >= 10 and
            x.genotype(f.child).data.GQ >= 10 and
            x.genotype(f.child).gt_type == 1 and
            x.genotype(f.father).gt_type == 0 and
            x.genotype(f.mother).gt_type == 0):
            print_line(denovo_out, x, f)
        if make_test_data:
            #if x.ID is None: # or (not record.FILTER is None and len(record.FILTER)>0):
            #    continue
            if (#x.INFO["MLEAC"][0] >= 5 and
                #x.INFO["MLEAC"][0] <= 2 and
                x.genotype(f.father).data.GQ >= min_test_parent_GQ and
                x.genotype(f.mother).data.GQ >= min_test_parent_GQ and
                x.genotype(f.father).data.GQ >= min_test_parent_DP and
                x.genotype(f.mother).data.GQ >= min_test_parent_DP and
                has_good_GQ(x, f.child) and
                #get_singleton_GQ(x) >= min_test_parent_GQ and
                x.genotype(f.father).gt_type == 0 and
                x.genotype(f.mother).gt_type == 0):
                if (homref_DP_count[x.genotype(f.child).data.DP] < min_test_lines):
                    homref_DP_count[x.genotype(f.child).data.DP] += 1
                    print_line(homoref_test_out, x, f, denovo=False)
            if (x.genotype(f.father).data.GQ >= min_test_parent_GQ and
                x.genotype(f.mother).data.GQ >= min_test_parent_GQ and
                x.genotype(f.father).data.GQ >= min_test_parent_DP and
                x.genotype(f.mother).data.GQ >= min_test_parent_DP and
                ((x.genotype(f.father).gt_type == 0 and
                 x.genotype(f.mother).gt_type == 2) or
                (x.genotype(f.father).gt_type == 2 and
                 x.genotype(f.mother).gt_type == 0))):
                if (het_DP_count[x.genotype(f.child).data.DP] < min_test_lines):
                    het_DP_count[x.genotype(f.child).data.DP] += 1
                    print_line(het_test_out, x, f, denovo=False)


denovo_out.close()
if make_test_data:
    het_test_out.close()
    homoref_test_out.close()
