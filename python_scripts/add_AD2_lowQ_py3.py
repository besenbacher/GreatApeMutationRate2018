import sys
import argparse
import pysam
from math import log
from math import exp

parser = argparse.ArgumentParser(description='''
Assign parent of origin to heterozygous variants in a child using read-backed phasing
''')

#parser.add_argument('variants', type=str, help='text file with variant positions.')
parser.add_argument('child_bamfile', type=str, help='fathers bamfile')
parser.add_argument('father_bamfile', type=str, help='fathers bamfile')
parser.add_argument('mother_bamfile', type=str, help='mothers bamfile')
parser.add_argument('--verbose', action='store_true', help="print progress output")

args = parser.parse_args()

childbam = pysam.Samfile(args.child_bamfile, 'rb')
fatherbam = pysam.Samfile(args.father_bamfile, 'rb')
motherbam = pysam.Samfile(args.mother_bamfile, 'rb')


def get_AD2(bamfile, chrom, pos, alt):
    '''Extract a specific column of the alignment and get the number of reads
    containing the alternative read'''
    pileup = bamfile.pileup(chrom, pos-1, pos)
    alt_MQs = []
    for col in pileup:
        if col.pos == pos-1:
            #for read in col.pileups:
            #    print(read.query_position)
            #    print(dir(read))
            alt_MQs = [read.alignment.mapping_quality for read in col.pileups if not read.query_position is None and read.alignment.seq[read.query_position] == alt]
    MQ0 = sum(1 for x in alt_MQs if x==0)
    MQlt20 = sum(1 for x in alt_MQs if x>0 and x<20)
    MQgt20lt40 = sum(1 for x in alt_MQs if x>=20 and x<40)
    MQgt40 = sum(1 for x in alt_MQs if x>=40)
    return (MQ0, MQlt20, MQgt20lt40, MQgt40)

header = sys.stdin.readline()

print(header.strip() +
      " CHILD.AD2.MQ0 CHILD.AD2.MQ1t20 CHILD.AD2.MQ20t39 CHILD.AD2.MQ40p" +
      " FATHER.AD2.MQ0 FATHER.AD2.MQ1t20 FATHER.AD2.MQ20t39 FATHER.AD2.MQ40p" +
      " MOTHER.AD2.MQ0 MOTHER.AD2.MQ1t20 MOTHER.AD2.MQ20t39 MOTHER.AD2.MQ40p")

for line in sys.stdin:
    L = line.split()
    chrom, pos, ref, alt, child = L[:5]
    print( line.strip() + '\t' +\
        '\t'.join(str(x) for x in get_AD2(childbam, chrom, int(pos), alt)) + '\t' + \
        '\t'.join(str(x) for x in get_AD2(fatherbam, chrom, int(pos), alt)) + '\t' + \
        '\t'.join(str(x) for x in get_AD2(motherbam, chrom, int(pos), alt)))
