import sys
ddir='/home/besen/MutationRates/NewDeNovoCalling'
f = open(ddir + '/avg_coverage.txt')
cov = {}
for line in f:
    child,c = line.split()
    cov[child] = c
f.close()

def read_family_description(fname):
    f = open(fname)
    res = []
    for line in f:
        L = line.split()
        res.append((L[0],L[1],L[2]))
    f.close()
    return res

CHIMP_FAMILIES = read_family_description(ddir +'/family_description/chimpanzees.txt')
GORILLA_FAMILIES = read_family_description(ddir +'/family_description/gorillas.txt')
ORANGUTAN_FAMILIES = read_family_description(ddir +'/family_description/orangutans.txt')

c2p = {}
for child, father, mother in CHIMP_FAMILIES + GORILLA_FAMILIES + ORANGUTAN_FAMILIES:
    if child in cov and father in cov and mother in cov:
        c2p[child] = (cov[child], cov[father], cov[mother])

header = sys.stdin.readline()
print(header.strip(), 'CHILD.COV', 'FATHER.COV', 'MOTHER.COV', sep='\t')

for line in sys.stdin:
    L = line.split()
    chrom, pos, ref, alt, child = L[:5]
    print(line.strip(), '\t'.join(c2p[child]), sep='\t')
