from gwf import *
from gwf_templates import *
import itertools

### Workflow for variant calling. ###

ORANGUTAN_REF = 'ponAbe2'
CHIMP_REF     = 'panTro5' 
GORILLA_REF   = 'gorGor4'

ORANGUTANS = ['Ogan',  'Aris', 'Schubbi']
CHIMPS     = ['Frits', 'Carolina', 'Carl', 'Simliki', 'ERR466113','ERR466114','ERR466115','ERR466116','ERR466117','ERR466118','ERR466119','ERR466120','ERR466121']
GORILLAS   = ['Banjo', 'Mimi', 'Mutasi', 'Mawenzi', 'Efata', 'Undi']

ALL_APES = ORANGUTANS + CHIMPS + GORILLAS

def read_family_description(fname):
    f = open(fname)
    res = []
    for line in f:
        L = line.split()
        res.append((L[0],L[1],L[2]))
    f.close()
    return res

CHIMP_FAMILIES = read_family_description('family_description/chimpanzees.txt')
GORILLA_FAMILIES = read_family_description('family_description/gorillas.txt')
ORANGUTAN_FAMILIES = read_family_description('family_description/orangutans.txt')
#Order: (Child, Father, Mother)

FAMILIES = []
FAMILIES.extend(['gorillas'] + list(family) for family in GORILLA_FAMILIES)
FAMILIES.extend(['chimpanzees'] + list(family) for family in CHIMP_FAMILIES)
FAMILIES.extend(['orangutans'] + list(family) for family in ORANGUTAN_FAMILIES)


#First we need to index the reference genomes with Picard
#refGenomes should be in ./ref-genomes

target('picard_index_orang')   << picard_index_reference(refGenome=ORANGUTAN_REF)
target('picard_index_chimp')   << picard_index_reference(refGenome=CHIMP_REF)
target('picard_index_gorilla') << picard_index_reference(refGenome=GORILLA_REF)

## Use GATK to process bam files before calling.

#To work with the bam files in GATK we first need them indexed...
# bamfiles should be in ./merged-bams

for individual in ALL_APES:
    target('samtools_index_bam_' + individual) << samtools_index_bam(individual=individual)
    target('samtools_index_filtered_bam_' + individual) << samtools_index_filtered_bam(individual=individual)


#Next, identify regions that needs to be locally re-aligned.
# intermediary output files from gatk will be put in ./gatk_files/

for individual in CHIMPS:
    target('collect_realign_regions_' + individual) << collect_realign_regions(refGenome=CHIMP_REF, individual=individual)

for individual in ORANGUTANS:
    target('collect_realign_regions_' + individual) << collect_realign_regions(refGenome=ORANGUTAN_REF, individual=individual)

for individual in GORILLAS:
    target('collect_realign_regions_' + individual) << collect_realign_regions(refGenome=GORILLA_REF, individual=individual)


#Realign the bam files:

for individual in CHIMPS:
    target('realign_' + individual) << realign(refGenome=CHIMP_REF, individual=individual)

for individual in ORANGUTANS:
    target('realign_' + individual) << realign(refGenome=ORANGUTAN_REF, individual=individual)

for individual in GORILLAS:
    target('realign_' + individual) << realign(refGenome=GORILLA_REF, individual=individual)


#Recalibrate base quality scores:

orang_known = ['../KnownVariants/abelii.2012-08-17.snps.dropsamples.autos.vqsr99.vcf',
               '../KnownVariants/abelii.2012-08-17.snps.dropsamples.chrX.vcf',
               '../KnownVariants/pygmaeus.2012-08-17.snps.dropsamples.autos.vqsr99.vcf',
               '../KnownVariants/pygmaeus.2012-08-17.snps.dropsamples.chrX.vcf']

chimp_known = ['../faststorage/KnownVariants/pantro.snps.2012-11-25.autos.99vqsr.beagle_' + CHIMP_REF +'_lift.vcf',
               '../faststorage/KnownVariants/gatk.allchimpanzee.2012-05-30.snps.chrX.vqsr99_' + CHIMP_REF + '_lift.vcf']
gorilla_known = ['../faststorage/KnownVariants/emitall-var.dropsamples.autos.vqsr99.2012-04-26.abfilter_' + GORILLA_REF  + '_lift_sorted.vcf']

for individual in CHIMPS:
    target('calc_recalibrate_info_' + individual) << \
        calc_recalibrate_info(individual=individual,
                              refGenome=CHIMP_REF,
                              knownlist=chimp_known)
    target('recalibrate_' + individual) << recalibrate(refGenome=CHIMP_REF, individual=individual)
    target('filter_bam_' + individual) << filter_bam_files(individual=individual,
                                                           refGenome=CHIMP_REF)


for individual in ORANGUTANS:
    target('calc_recalibrate_info_' + individual) << \
        calc_recalibrate_info(individual=individual,
                              refGenome=ORANGUTAN_REF,
                              knownlist=orang_known)
    target('recalibrate_' + individual) << recalibrate(refGenome=ORANGUTAN_REF, individual=individual)
    target('filter_bam_' + individual) << filter_bam_files(individual=individual,
                                                           refGenome=ORANGUTAN_REF)

for individual in GORILLAS:
    target('calc_recalibrate_info_' + individual) << \
        calc_recalibrate_info(individual=individual,
                              refGenome=GORILLA_REF,
                              knownlist=gorilla_known)
    target('recalibrate_' + individual) << recalibrate(refGenome=GORILLA_REF, individual=individual)
    target('filter_bam_' + individual) << filter_bam_files(individual=individual,
                                                          refGenome=GORILLA_REF)


#Call variants using haplotype caller

def get_chromosomes(refname):
    f = open('ref-genomes/'+ refname + '.fa.fai')
    chroms = []
    for line in f:
        chrom = line.split()[0]
        #if 'random' not in chrom and
        if 'Un' not in chrom:
            chroms.append(chrom)
    f.close()
    return chroms

def is_autosome(chrom):
    return (not chrom.startswith('chrX') and
            not chrom.startswith('chrY') and
            not chrom.startswith('chrUn') and
            not chrom.startswith('chrM'))

def get_regions(refname):
    f = open('split_files/' + refname + '_split_1000_1000000.txt')
    L = []
    for line in f:
        chrom, start, end = line.split()
        L.append((chrom, int(start), int(end)))
    f.close()
    return L

CHIMP_CHROMOSOMES = get_chromosomes(CHIMP_REF)
ORANGUTAN_CHROMOSOMES = get_chromosomes(ORANGUTAN_REF)
GORILLA_CHROMOSOMES = get_chromosomes(GORILLA_REF)
ORANGUTAN_REGIONS = get_regions(ORANGUTAN_REF)
CHIMP_REGIONS = get_regions(CHIMP_REF)
GORILLA_REGIONS = get_regions(GORILLA_REF)

ORANGUTAN_AUTOSOMES = [x for x in ORANGUTAN_CHROMOSOMES if is_autosome(x)]
GORILLA_AUTOSOMES = [x for x in GORILLA_CHROMOSOMES if is_autosome(x)]
CHIMP_AUTOSOMES = [x for x in CHIMP_CHROMOSOMES if is_autosome(x)]

ORANGUTAN_AUTOSOME_REGIONS = [x for x in ORANGUTAN_REGIONS if is_autosome(x[0])]
GORILLA_AUTOSOME_REGIONS = [x for x in GORILLA_REGIONS if is_autosome(x[0])]
CHIMP_AUTOSOME_REGIONS = [x for x in CHIMP_REGIONS if is_autosome(x[0])]

# CHIMP_REF_2BIT should be path to chimp reference in 2bit format.
# 2bit files can be downloaded from UCSC:
# http://hgdownload.soe.ucsc.edu/goldenPath/panTro5/bigZips/panTro5.2bit
CHIMP_REF_2BIT = '~/Data/2bit/' + CHIMP_REF + '.2bit'
GORILLA_REF_2BIT = '~/Data/2bit/' +  GORILLA_REF + '.2bit'
ORANGUTAN_REF_2BIT = '~/Data/2bit/' + ORANGUTAN_REF + '.2bit'



for chrom, start, end in ORANGUTAN_REGIONS:
    target('hapcaller_bo_orangutan_' + chrom + '_' + str(start)) << \
      haplotype_caller_region_bo(refGenome=ORANGUTAN_REF, specie="orangutans",
                                 individuals=ORANGUTANS, chrom=chrom, start=start, end=end)


for chrom, start, end in GORILLA_REGIONS:
    target('hapcaller_bo_gorilla_' + chrom + '_' + str(start)) << \
      haplotype_caller_region_bo(refGenome=GORILLA_REF, specie="gorillas",
                                 individuals=GORILLAS, chrom=chrom, start=start, end=end)


for chrom, start, end in CHIMP_REGIONS:
    target('hapcaller_bo_chimps_' + chrom + '_' + str(start)) << \
      haplotype_caller_region_bo(refGenome=CHIMP_REF, specie="chimpanzees",
                               individuals=CHIMPS, chrom=chrom, start=start, end=end)


target('combine_regions_chimps') << combine_regions(refname=CHIMP_REF, specie="chimpanzees", regions=CHIMP_REGIONS)
target('combine_regions_gorillas') << combine_regions(refname=GORILLA_REF, specie="gorillas", regions=GORILLA_REGIONS)
target('combine_regions_orangutans') << combine_regions(refname=ORANGUTAN_REF, specie="orangutans", regions=ORANGUTAN_REGIONS)

for individual in ALL_APES:
    target('get_rg_' + individual) << get_readgroup(individual=individual)

for individual in set(itertools.chain(*CHIMP_FAMILIES)):
    for chrom, start, end in CHIMP_REGIONS:
        target('split_bo_' +individual + '_' +chrom +'_' + str(start)) << \
          split_bamout_region(individual=individual, chrom=chrom, start=start, end=end, specie="chimpanzees")

for individual in set(itertools.chain(*GORILLA_FAMILIES)):
  for chrom, start, end in GORILLA_REGIONS:
      target('split_bo_' +individual + '_' +chrom +'_' + str(start)) << \
          split_bamout_region(individual=individual, chrom=chrom, start=start, end=end, specie="gorillas")

for individual in set(itertools.chain(*ORANGUTAN_FAMILIES)):
    for chrom, start, end in ORANGUTAN_REGIONS:
        target('split_bo_' +individual + '_' +chrom +'_' + str(start)) << \
            split_bamout_region(individual=individual, chrom=chrom, start=start, end=end, specie="orangutans")


minQ =20
for child, father, mother in CHIMP_FAMILIES:
    for chrom in CHIMP_AUTOSOMES:
        target('family_depth_chimp_' + chrom +'_' + child) << \
            get_family_coverage_context_wr(father=father, mother=mother, child=child,
                                           minQ=minQ, minq=10, refname=CHIMP_REF, chrom=chrom)
    target('combine_coverage_' + child + '_ autosomes') << \
        combine_coverage_context_wr(child=child, chromosomes=CHIMP_AUTOSOMES,
                                    outname='autosomes', minQ=minQ, minq=10)

for child, father, mother in GORILLA_FAMILIES:
    for chrom in GORILLA_AUTOSOMES:
        target('family_depth_gorilla_' + chrom +'_' + child) << \
            get_family_coverage_context_wr(father=father, mother=mother, child=child,
                                           minQ=minQ, minq=10, refname=GORILLA_REF, chrom=chrom)
    target('combine_coverage_' + child + '_ autosomes') << \
        combine_coverage_context_wr(child=child, chromosomes=GORILLA_AUTOSOMES,
                                    outname='autosomes', minQ=minQ, minq=10)

for child, father, mother in ORANGUTAN_FAMILIES:
    for chrom in ORANGUTAN_AUTOSOMES:
        target('family_depth_orangutan_' + chrom +'_' + child) << \
            get_family_coverage_context_wr(father=father, mother=mother, child=child,
                                           minQ=minQ, minq=10, refname=ORANGUTAN_REF, chrom=chrom)
    target('combine_coverage_' + child + '_ autosomes') << \
        combine_coverage_context_wr(child=child, chromosomes=ORANGUTAN_AUTOSOMES,
                                    outname='autosomes', minQ=minQ, minq=10)


### Parent of Origin Assignment ###
# This part assigns Parent of Origin for all het variants in a child
# If you don't have bam-out bamfiles you should just use the filtered bam-files

for child, father, mother in ORANGUTAN_FAMILIES:
    for chrom, start, end in ORANGUTAN_REGIONS:
        target('samtools_index_bamout_' + child + '_' + chrom +'_' + str(start) ) << \
            samtools_index_bamout(individual=child, chrom=chrom, start=start, end=end)
        target('run_poo_' + child + '_' + chrom +'_' + str(start)) << \
            run_poo_region(child=child, father=father, mother=mother, chrom=chrom, start=start, end=end, specie='orangutans')
    target('combine_poo_regions_'+ child) << \
        combine_poo_region(child=child, regions=ORANGUTAN_AUTOSOME_REGIONS, outname='autosomes')


for child, father, mother in GORILLA_FAMILIES:
    for chrom, start, end in GORILLA_REGIONS:
        target('samtools_index_bamout_' + child + '_' + chrom +'_' + str(start) ) << \
            samtools_index_bamout(individual=child, chrom=chrom, start=start, end=end)
        target('run_poo_' + child + '_' + chrom +'_' + str(start)) << \
            run_poo_region(child=child, father=father, mother=mother, chrom=chrom, start=start, end=end, specie='gorillas')
    target('combine_poo_regions_'+ child) << \
        combine_poo_region(child=child, regions=GORILLA_AUTOSOME_REGIONS, outname='autosomes')

for child, father, mother in CHIMP_FAMILIES:
    for chrom, start, end in CHIMP_REGIONS:
        target('samtools_index_bamout_' + child + '_' + chrom +'_' + str(start) ) << \
            samtools_index_bamout(individual=child, chrom=chrom, start=start, end=end)
        target('run_poo_' + child + '_' + chrom +'_' + str(start)) << \
            run_poo_region(child=child, father=father, mother=mother, chrom=chrom, start=start, end=end, specie='chimpanzees')
    target('combine_poo_regions_'+ child) << \
        combine_poo_region(child=child, regions=CHIMP_AUTOSOME_REGIONS, outname='autosomes')


### Make input files for rate estimation script ###

# This part creates three files for each specie:
# denovo_raw_SNV.dat: contains all variants where the child is heterozygous and the parents are homozygous for the reference.
# het_test_SNV.dat: file for testing the callability of heterozygous variants
# homoref_test_SNV.dat: file for testing the callability of homozygous variants

for child,father,mother in CHIMP_FAMILIES:
    target('vcf2denovo_dat_' + child) << \
        vcf2denovo_dat_child(specie='chimpanzees', vtype='SNV', refname=CHIMP_REF_2BIT, child=child)
for child,father,mother in GORILLA_FAMILIES:
    target('vcf2denovo_dat_' + child) << \
        vcf2denovo_dat_child(specie='gorillas', vtype='SNV', refname=GORILLA_REF_2BIT, child=child)
for child,father,mother in ORANGUTAN_FAMILIES:
    target('vcf2denovo_dat_' + child) << \
        vcf2denovo_dat_child(specie='orangutans', vtype='SNV', refname=ORANGUTAN_REF_2BIT, child=child)


#This part adds extra columns to the files created above

knownD = {
    'gorillas':gorilla_known,
    'orangutans':orang_known,
    'chimpanzees':chimp_known,
}

refD = {
    'gorillas':GORILLA_REF_2BIT,
    'orangutans':ORANGUTAN_REF_2BIT,
    'chimpanzees':CHIMP_REF_2BIT,
}

for specie in knownD:
    target('make_known_' + specie) << \
        make_known(specie=specie, caller='gatk', known=knownD[specie])

for specie, child, father, mother in FAMILIES:
    target('add_known_' + child) << \
        add_known(specie=specie, caller='gatk', known=knownD[specie], child=child)


minQ=20
for specie, child, father, mother in FAMILIES:
    for var_type in ['SNV']:
        target('get_sam_depth_variants_' + child + '_' + var_type + '_' + str(minQ)) << \
          get_sam_depth_variants(specie=specie,
                                 child=child,
                                 father=father,
                                 mother=mother,
                                 vtype=var_type,
                                 caller="gatk",
                                 minQ=minQ,
                                 minq=10)


vtype='SNV'
for specie, child, father, mother in FAMILIES:
    for ftype in ['denovo_raw', 'het_test', 'homoref_test']:
        ftype += '_' + vtype
        target('merge_and_join_'+ child + '_' + ftype) << \
            merge_and_join(child=child, specie=specie, caller="gatk",
                           ftype=ftype, minQ=20, minq=10, extra='')
    for ftype,extra in [('denovo_raw','_w_DPS20_w_known'), ('het_test','_w_DPS20'),('homoref_test','_w_DPS20')]:
        ftype += '_' + vtype
        target('addAD2_lowMQ_'+ child + '_' + ftype + '_' + str(minQ)) << \
            addAD2_lowMQ(child=child, father=father, mother=mother,
                         specie=specie, caller="gatk", ftype=ftype, extra=extra)
    for ftype,extra in [('denovo_raw','_w_DPS20_w_known_w_lowMQ'), ('het_test','_w_DPS20_w_lowMQ'),('homoref_test','_w_DPS20_w_lowMQ')]:
        ftype += '_' + vtype
        target('addPC_'+ child + '_' + ftype + '_' + str(minQ)) << \
            add_parents_cov(child=child, father=father, mother=mother,
                            specie=specie, caller="gatk", ftype=ftype, extra=extra)

    for ftype,extra in [('denovo_raw','_w_DPS20_w_known_w_lowMQ_w_pc'), ('het_test','_w_DPS20_w_lowMQ_w_pc'),('homoref_test','_w_DPS20_w_lowMQ_w_pc')]:
        ftype += '_' + vtype
        target('add_repeat_'+ child + '_' + ftype + '_' + str(minQ)) << \
            add_repeat(child=child, father=father, mother=mother, refgenome=refD[specie],
                       specie=specie, caller="gatk", ftype=ftype, extra=extra)
    
    for ftype,extra in [('denovo_raw','_w_DPS20_w_known_w_lowMQ_w_pc_w_repeat'), ('het_test','_w_DPS20_w_lowMQ_w_pc_w_repeat'),('homoref_test','_w_DPS20_w_lowMQ_w_pc_w_repeat')]:
        ftype += '_' + vtype
        target('add_segdup_'+ child + '_' + ftype + '_' + str(minQ)) << \
            add_segdup(child=child, father=father, mother=mother, refgenome=refD[specie],
                       specie=specie, caller="gatk", ftype=ftype, extra=extra)

### This part finds de novo mutations and calculates rates given a set of cutoffs. ###

calc_rate_family = \
    template(input=['dat_files/{caller}/{specie}/denovo_raw_SNV_{child}_w_DPS20_w_known_w_lowMQ_w_pc_w_repeat_w_sd.dat',
                    'dat_files/{caller}/{specie}/het_test_SNV_w_DPS20_w_lowMQ_w_pc_w_repeat_w_sd.dat',
                    'dat_files/{caller}/{specie}/homoref_test_SNV_w_DPS20_w_lowMQ_w_pc_w_repeat_w_sd.dat'],
             output=["new_results_split_lowMQAD2_w_pc_w_repeat_noSD/{specie}/{child}/HOM.GQ_{HOMGQ}_HET.GQ_{HETGQ}_P.AD2_{AD2}_MAX.AR_{AR}_minDP_{minDP}_maxDP_{maxDP}_lowMQAD2_{lowMQAD2}/avg_rate.txt"],
             account="MutationRates",walltime="2:00:00",memory="6g") << \
'''
mkdir -p new_results_split_lowMQAD2_w_pc_w_repeat_noSD/{specie}/{child}/HOM.GQ_{HOMGQ}_HET.GQ_{HETGQ}_P.AD2_{AD2}_MAX.AR_{AR}_minDP_{minDP}_maxDP_{maxDP}_lowMQAD2_{lowMQAD2}/
./R_scripts/calc_rate_family.R {specie} {child} new_results_split_lowMQAD2_w_pc_w_repeat_noSD/{specie}/{child}/HOM.GQ_{HOMGQ}_HET.GQ_{HETGQ}_P.AD2_{AD2}_MAX.AR_{AR}_minDP_{minDP}_maxDP_{maxDP}_lowMQAD2_{lowMQAD2}/ {HOMGQ} {HETGQ} {AD2} {AR} {minDP} {maxDP} {lowMQAD2} {minRC} {maxRC} {minRPRS} {maxRPRS}
'''

RPRS = {'chimpanzees':(-2.5,2.5), 'gorillas': (-2.5,2.5), 'orangutans':(-2.5,2.5)}

AD2 = 0
lowMQAD2 = 1
minRC = 0.0
maxDP = 1000
AR = 1.0

# Test a grid over different quality parameters
for GQ in range(20,100,5):
    for minDP in [5,10]:
        for maxRC in [1.5, 1.7, 1.9]:
            for specie, child, father, mother in FAMILIES:
                minRPRS,maxRPRS = RPRS[specie]
                HOMGQ = GQ
                HETGQ = GQ
                target('rate_%s_%d_%d_%d_%.2f_%d_%d_%d'% (child, HOMGQ, HETGQ, AD2, AR, minDP, maxDP, lowMQAD2)) << \
                    calc_rate_family(caller='gatk', specie=specie, child=child,
                                     HOMGQ=HOMGQ, HETGQ=HETGQ, AD2=AD2, AR=AR, 
                                     minDP=minDP, maxDP=maxDP,lowMQAD2=lowMQAD2, 
                                     minRPRS=minRPRS, maxRPRS=maxRPRS, minRC=minRC, maxRC=maxRC)


