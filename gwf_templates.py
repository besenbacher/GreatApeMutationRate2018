from gwf import *

### Helper function to add extra 'input' and 'output' files to a options dict
def add_options(options, extra):
    for inout in ['input', 'output']:
        if inout in extra:
            if inout in options:
                if type(options[inout]) != list:
                    options[inout] = [options[inout]]
                if type(extra[inout]) != list:
                    extra[inout] = [extra[inout]]
                options[inout] = list(set(options[inout] + extra[inout]))
            else:
                options[inout] = extra[inout]
    return options


picard_index_reference = \
    template(input='ref-genomes/{refGenome}.fa',
             output=['ref-genomes/{refGenome}.dict','ref-genomes/{refGenome}.fa.fai'],
             account="DanishPanGenome",walltime="1:00:00") << '''
samtools faidx ref-genomes/{refGenome}.fa
picard CreateSequenceDictionary R=ref-genomes/{refGenome}.fa O=ref-genomes/{refGenome}.dict
'''

samtools_index_bam = \
    template(input='merged_bam_files/{individual}.bam',
             output='merged_bam_files/{individual}.bam.bai',
            walltime='59:00', account='MutationRates') << '''
samtools index merged_bam_files/{individual}.bam
'''

samtools_index_filtered_bam = \
    template(input='filtered_bam_files/{individual}.bam',
             output='filtered_bam_files/{individual}.bam.bai',
             walltime='11:00:00', account='MutationRates') << '''
samtools index filtered_bam_files/{individual}.bam
'''

samtools_index_bamout = \
    template(input='bamout_files/{individual}/{chrom}.{start}.{end}.bam',
            output='bamout_files/{individual}/{chrom}.{start}.{end}.bam.bai',
            cores=1, account='MutationRates', walltime='0:59:00',memory='1g') << '''
samtools index bamout_files/{individual}/{chrom}.{start}.{end}.bam
'''

collect_realign_regions = \
  template(input=['ref-genomes/{refGenome}.fa', 'ref-genomes/{refGenome}.dict',
                  'merged_bam_files/{individual}.bam', 'merged_bam_files/{individual}.bam.bai'],
           output=['new_gatk_files/{individual}.realign.intervals'],
           cores=16, account="DanishPanGenome") << '''
gatk -nt 16 -T RealignerTargetCreator -R ref-genomes/{refGenome}.fa \
     -I merged_bam_files/{individual}.bam \
     -o new_gatk_files/{individual}.realign.intervals
'''

realign = template(input=['ref-genomes/{refGenome}.fa', 'ref-genomes/{refGenome}.dict',
                          'merged_bam_files/{individual}.bam', 'merged_bam_files/{individual}.bam.bai',
                          'new_gatk_files/{individual}.realign.intervals'],
                   output=['new_gatk_files/{individual}.realigned.bam','new_gatk_files/{individual}.realigned.bai'],
                   cores=16, account="DanishPanGenome") << '''
gatk -T IndelRealigner -R ref-genomes/{refGenome}.fa \
     --filter_bases_not_stored \
     -I merged_bam_files/{individual}.bam \
     -targetIntervals new_gatk_files/{individual}.realign.intervals \
     -o new_gatk_files/{individual}.realigned.bam
'''


_calc_recalibrate_info = template(input=['ref-genomes/{refGenome}.fa', 'ref-genomes/{refGenome}.dict',
                                         'new_gatk_files/{individual}.realigned.bam', 'new_gatk_files/{individual}.realigned.bai'],
                                  output=['new_gatk_files/{individual}.recalibration_report.grp'],
                                  cores=16, account="DanishPanGenome") << '''
gatk -T BaseRecalibrator \
     -nct 16 \
     -R ref-genomes/{refGenome}.fa \
     {knownlist}  \
     -I new_gatk_files/{individual}.realigned.bam \
     -o new_gatk_files/{individual}.recalibration_report.grp
'''

def calc_recalibrate_info(**arguments):
    arguments['knownlist'] = ' '.join('-knownSites ' + x for x in arguments['knownlist'])
    return _calc_recalibrate_info(**arguments)

recalibrate = template(input=['ref-genomes/{refGenome}.fa', 'ref-genomes/{refGenome}.dict',
                              'new_gatk_files/{individual}.realigned.bam', 'new_gatk_files/{individual}.realigned.bai',
                              'new_gatk_files/{individual}.recalibration_report.grp'],
                       output=['new_gatk_files/{individual}.recalibrated.bam', 'new_gatk_files/{individual}.recalibrated.bai'],
                       cores=16, memory="50g",account="DanishPanGenome") << '''
java -Djava.io.tmpdir=/scratch/$PBS_JOBID/tmp \
     -Djava.awt.headless=true \
     -Xmx50g \
     -jar /com/extra/GATK/3.8/jar-bin/GenomeAnalysisTK.jar \
     -T PrintReads \
     -nct 16 \
     -R ref-genomes/{refGenome}.fa \
     -I new_gatk_files/{individual}.realigned.bam \
     -BQSR new_gatk_files/{individual}.recalibration_report.grp \
     -o new_gatk_files/{individual}.recalibrated.bam
'''


filter_bam_files = \
    template(\
        input=['new_gatk_files/{individual}.recalibrated.bam',
               'new_gatk_files/{individual}.recalibrated.bai',
               'ref-genomes/{refGenome}.fa'],
        output=['filtered_bam_files/{individual}.bam'],
        account="DanishPanGenome",
        walltime="12:00:00", cores=16) << \
'''
gatk \
 -R  ref-genomes/{refGenome}.fa \
 -T PrintReads \
 -I new_gatk_files/{individual}.recalibrated.bam \
 -o filtered_bam_files/{individual}.bam \
 -nct 16 \
 --read_filter BadCigar \
 --read_filter DuplicateRead \
 --read_filter FailsVendorQualityCheck \
 --read_filter HCMappingQuality \
 --read_filter MappingQualityUnavailable \
 --read_filter NotPrimaryAlignment \
 --read_filter UnmappedRead \
 --filter_bases_not_stored \
 --filter_mismatching_base_and_quals
'''



_haplotype_caller_region_bo = template(input=['ref-genomes/{refGenome}.fa'],
                                   output=['new_vcf_files/{specie}/haplocaller.raw.{chrom}.{start}.{end}.vcf', 'new_vcf_files/{specie}/haplocaller.raw.{chrom}.{start}.{end}.vcf.idx', 'new_gatk_files/bamout/{specie}/{chrom}.{start}.{end}.bam'],
                                   cores=1,
                                   memory="16g", account="DanishPanGenome",
                                   walltime="160:00:00") <<  '''
java -Djava.io.tmpdir=/scratch/$PBS_JOBID/tmp \
     -Djava.awt.headless=true \
     -Xmx16g \
     -jar /com/extra/GATK/3.8/jar-bin/GenomeAnalysisTK.jar \
     -T HaplotypeCaller \
     -nct 1 \
     -R 'ref-genomes/{refGenome}.fa' \
     {includelist} \
     -A DepthPerSampleHC \
     -A Coverage \
     -A HaplotypeScore \
     -A StrandAlleleCountsBySample \
     -L {chrom}:{start}-{end} \
     -bamout new_gatk_files/bamout/{specie}/{chrom}.{start}.{end}.bam \
     -o new_vcf_files/{specie}/haplocaller.raw.{chrom}.{start}.{end}.vcf
'''


def haplotype_caller_region_bo(**arguments):
    assert 'individuals' in arguments
    bamfiles = ['new_gatk_files/'+individual+'.recalibrated.bam' for individual in arguments['individuals']]
    extra_options = {'input':bamfiles}
    arguments['includelist'] = ' '.join('-I ' + x for x in bamfiles)
    (options, spec) = _haplotype_caller_region_bo(**arguments)
    return (add_options(options,extra_options), spec)


_combine_regions = template(input=[],
                                output=['new_vcf_files/{specie}.haplocaller.raw.auto.vcf'],
                                walltime="10:00:00", account="DanishPanGenome")  << '''#
vcf-concat $(cat split_files/{refname}_split_1000_1000000.txt | fgrep -v chrX | fgrep -v chrUn | fgrep -v chrM | fgrep -v chrY | fgrep -v hap | ~/Scripts/gorsort.sh | awk -vORS=" " '{{print "new_vcf_files/{specie}/haplocaller.raw."$1"."$2"."$3".vcf"}}') > new_vcf_files/{specie}.haplocaller.raw.auto.vcf
'''

def combine_regions(**arguments):
    vcf_files = ['new_vcf_files/' + arguments["specie"] + '/haplocaller.raw.' + chrom + '.' + str(start) + '.' + str(end) + '.vcf' for (chrom,start,end) in arguments['regions']]
    extra_options = {'input':vcf_files}
    (options, spec) = _combine_regions(**arguments)
    return (add_options(options, extra_options), spec)


get_readgroup = \
    template(input=['merged_bam_files/{individual}.rg.txt'],
             output=['read_groups/{individual}.rg.txt'],
             account="DanishPanGenome",
             walltime="0:59:00",memory="1g") << '''
cat merged_bam_files/{individual}.rg.txt | cut -f2 | cut -d: -f2 > read_groups/{individual}.rg.txt
'''

split_bamout_region = \
    template(input=['read_groups/{individual}.rg.txt','new_gatk_files/bamout/{specie}/{chrom}.{start}.{end}.bam'],
             output=['bamout_files/{individual}/{chrom}.{start}.{end}.bam','bamout_files/{individual}/{chrom}.{start}.{end}.bam.bai'],
             walltime="10:00:00",
             memory="4g") << '''
mkdir -p bamout_files/{individual}
samtools view -R read_groups/{individual}.rg.txt -b -o bamout_files/{individual}/{chrom}.{start}.{end}.bam new_gatk_files/bamout/{specie}/{chrom}.{start}.{end}.bam
samtools index bamout_files/{individual}/{chrom}.{start}.{end}.bam
'''


#get count of coverage combinations in a trio stratified by 3bp-context and repeat status
get_family_coverage_context_wr = \
  template(input=['filtered_bam_files/{father}.bam','filtered_bam_files/{father}.bam.bai',
                  'filtered_bam_files/{mother}.bam','filtered_bam_files/{mother}.bam.bai',
                  'filtered_bam_files/{child}.bam','filtered_bam_files/{child}.bam.bai'],
           output=['new_family_coverage_wr/{child}/minQ{minQ}_minq{minq}_type_{chrom}.txt'],
           account="DanishPanGenome", walltime='48:00:00') << '''
mkdir -p new_family_coverage_wr/{child}/
samtools depth -Q{minQ} -q{minq} filtered_bam_files/{father}.bam filtered_bam_files/{mother}.bam filtered_bam_files/{child}.bam -r{chrom}| awk '$3>=5 && $3<=200 && $4>=5 && $4<=200 && $5>=5 && $5<=200' | ./python_scripts/add_context_type_w_repeat.py {refname} | cut -f3- | ./python_scripts/table.py > new_family_coverage_wr/{child}/minQ{minQ}_minq{minq}_type_{chrom}.txt
'''

_combine_coverage_context_wr = \
  template(output=['new_family_coverage_wr/{child}_minQ{minQ}_minq{minq}_type_{outname}.txt'],
          memory='8g',
          account='DanishPanGenome',
          walltime='10:00:00') << \
'''
cat {input_list} | ./python_scripts/combine_tables.py > new_family_coverage_wr/{child}_minQ{minQ}_minq{minq}_type_{outname}.txt
'''

def combine_coverage_context_wr(**arguments):
    chromosomes = arguments['chromosomes']
    input_files = ['new_family_coverage_wr/{child}/minQ{minQ}_minq{minq}_type'.format(**arguments) + '_' + chrom + '.txt' for chrom in chromosomes]
    arguments['input_list'] = ' '.join(input_files)
    extra_options = {'input':input_files}
    (options, spec) = _combine_coverage_context_wr(**arguments)
    return (add_options(options, extra_options), spec)


run_poo_region = \
    template(input=['new_vcf_files/{specie}/haplocaller.raw.{chrom}.{start}.{end}.vcf',
                    'bamout_files/{child}/{chrom}.{start}.{end}.bam',
                    'bamout_files/{child}/{chrom}.{start}.{end}.bam.bai'],
             output=['poo_data/{child}/{chrom}.{start}.{end}.txt'],
             walltime='10:00:00', account='DanishPanGenome',memory='1g') << '''
mkdir -p poo_data/{child}/
python python_scripts/vcf_parent_of_origin_py3.py new_vcf_files/{specie}/haplocaller.raw.{chrom}.{start}.{end}.vcf {father} {mother} {child} bamout_files/{child}/{chrom}.{start}.{end}.bam > poo_data/{child}/{chrom}.{start}.{end}.txt
'''

_combine_poo = \
  template(output=['poo_data/{child}_{outname}.txt'],
          memory='8g',
          account='DanishPanGenome',
          walltime='10:00:00') << \
'''
cat {input_list} | sort -k1,1 > poo_data/{child}_{outname}.txt
'''

def combine_poo_region(**arguments):
    regions = arguments['regions']
    input_files = ['poo_data/{child}/'.format(**arguments) + '{}.{}.{}.txt'.format(chrom,start,end) for chrom,start,end in regions]
    arguments['input_list'] = ' '.join(input_files)
    extra_options = {'input':input_files}
    (options, spec) = _combine_poo(**arguments)
    return (add_options(options, extra_options), spec)



vcf2denovo_dat_child = \
    template(input=["new_vcf_files/{specie}.haplocaller.raw.auto.vcf",
                    "family_description/{specie}.txt",
                    "poo_data/{child}_autosomes.txt"],
             output=["dat_files/gatk/{specie}/denovo_raw_{vtype}_{child}.dat",
                     "dat_files/gatk/{specie}/het_test_{vtype}_{child}.dat",
                     "dat_files/gatk/{specie}/homoref_test_{vtype}_{child}.dat"],
                 walltime='11:59:00', memory='2g',account='MutationRates') << \
'''
mkdir -p dat_files/gatk/{specie}/
cat new_vcf_files/{specie}.haplocaller.raw.auto.vcf | python python_scripts/get_denovo_single.py family_description/{specie}.txt dat_files/gatk/{specie}/ {refname} {child} --var_type {vtype} --postfix {vtype}_{child} --PoO_data poo_data/{child}_autosomes.txt
'''


_make_known = \
    template(input=[],
            output=['tmp/{specie}_known.txt'],
            walltime='10:00:00', memory='4g', account='DanishPanGenome') << '''

cat {known_files} | awk '$1!~/#/{{print $1"_"$2"_"$5,"TRUE"}}' | sort -u -k1,1 > tmp/{specie}_known.txt
'''

add_known = template(input=['dat_files/{caller}/{specie}/denovo_raw_SNV_{child}_w_DPS20.dat',
                             'tmp/{specie}_known.txt'],
                      output=['dat_files/{caller}/{specie}/denovo_raw_SNV_{child}_w_DPS20_w_known.dat'],
                      walltime='10:00:00', memory='4g', account='DanishPanGenome') << '''
n=$(head -1 dat_files/{caller}/{specie}/denovo_raw_SNV_{child}_w_DPS20.dat | awk '{{print NF}}')
cat dat_files/{caller}/{specie}/denovo_raw_SNV_{child}_w_DPS20.dat | awk '{{print $1"_"$2"_"$4,$0}}' | sort -k1,1 | join -a1 - tmp/{specie}_known.txt | cut -d" " -f2- | awk '{{if (NF=='$n') {{if ($1=="CHROM") print $0,"known_variant"; else print $0,"FALSE"}} else print $0}}' | sort -gk2 > dat_files/{caller}/{specie}/denovo_raw_SNV_{child}_w_DPS20_w_known.dat
'''

def make_known(**arguments):
    extra_options = {'input':arguments["known"]}
    arguments["known_files"] = ' '.join(arguments["known"])
    (options, spec) = _make_known(**arguments)
    return (add_options(options, extra_options), spec)


get_sam_depth_variants = \
  template(input = \
           ['dat_files/{caller}/{specie}/denovo_raw_{vtype}_{child}.dat',
            'dat_files/{caller}/{specie}/het_test_{vtype}_{child}.dat',
            'dat_files/{caller}/{specie}/homoref_test_{vtype}_{child}.dat',
            'filtered_bam_files/{father}.bam',
            'filtered_bam_files/{mother}.bam',
            'filtered_bam_files/{child}.bam'],
           output = \
            ["sam_depth/{specie}/depth_Q{minQ}_q{minq}_denovo_raw_{vtype}.{child}.txt",
             "sam_depth/{specie}/depth_Q{minQ}_q{minq}_het_test_{vtype}.{child}.txt",
             "sam_depth/{specie}/depth_Q{minQ}_q{minq}_homoref_test_{vtype}.{child}.txt"],
           account="MutationRates", walltime="110:00:00",memory="12g") << \
'''
mkdir -p sam_depth/{specie}
samtools depth -b <(cat dat_files/{caller}/{specie}/denovo_raw_{vtype}_{child}.dat | awk '$5=="{child}" {{print $1,$2-1,$2}}' | gsort -k1,1 -k2,2n) -Q{minQ} -q{minq} filtered_bam_files/{child}.bam filtered_bam_files/{father}.bam filtered_bam_files/{mother}.bam | awk '{{print $1"_"$2"_{child}",$3,$4,$5}}' | sort -k1,1 > sam_depth/{specie}/depth_Q{minQ}_q{minq}_denovo_raw_{vtype}.{child}.txt
samtools depth -b <(cat dat_files/{caller}/{specie}/homoref_test_{vtype}_{child}.dat | awk '$5=="{child}" {{print $1,$2-1,$2}}' | gsort -k1,1 -k2,2n) -Q{minQ} -q{minq} filtered_bam_files/{child}.bam filtered_bam_files/{father}.bam filtered_bam_files/{mother}.bam | awk '{{print $1"_"$2"_{child}",$3,$4,$5}}' | sort -k1,1 > sam_depth/{specie}/depth_Q{minQ}_q{minq}_homoref_test_{vtype}.{child}.txt
samtools depth -b <(cat dat_files/{caller}/{specie}/het_test_{vtype}_{child}.dat | awk '$5=="{child}" {{print $1,$2-1,$2}}' | gsort -k1,1 -k2,2n) -Q{minQ} -q{minq} filtered_bam_files/{child}.bam filtered_bam_files/{father}.bam filtered_bam_files/{mother}.bam | awk '{{print $1"_"$2"_{child}",$3,$4,$5}}' | sort -k1,1 > sam_depth/{specie}/depth_Q{minQ}_q{minq}_het_test_{vtype}.{child}.txt
'''

merge_and_join = \
    template(
        input=['dat_files/{caller}/{specie}/{ftype}_{child}{extra}.dat',
               'sam_depth/{specie}/depth_Q{minQ}_q{minq}_{ftype}.{child}.txt'],
        output=['dat_files/{caller}/{specie}/{ftype}_{child}{extra}_w_DPS{minQ}.dat'],
        memory='32g', account='MutationRates', walltime='11:00:00') << \
'''
head -1 dat_files/{caller}/{specie}/{ftype}_{child}{extra}.dat | awk -v OFS="\t" '{{print $0,"CHILD.DPS{minQ}","FATHER.DPS{minQ}","MOTHER.DPS{minQ}"}}' > dat_files/{caller}/{specie}/{ftype}_{child}{extra}_w_DPS{minQ}.dat
cat dat_files/{caller}/{specie}/{ftype}_{child}{extra}.dat | awk '$1!="CHROM"{{print $1"_"$2"_"$5, $0}}' | sort -k1,1 | join - <(sort -k1,1 sam_depth/{specie}/depth_Q{minQ}_q{minq}_{ftype}.{child}.txt) | sed 'y/ /\t/' | cut -f2- >>  dat_files/{caller}/{specie}/{ftype}_{child}{extra}_w_DPS{minQ}.dat
'''

addAD2_lowMQ = \
    template(
        input=['dat_files/{caller}/{specie}/{ftype}_{child}{extra}.dat'],
        output=['dat_files/{caller}/{specie}/{ftype}_{child}{extra}_w_lowMQ.dat'],
        memory='4g', account='MutationRates', walltime='24:00:00') << \
'''
cat dat_files/{caller}/{specie}/{ftype}_{child}{extra}.dat | python python_scripts/add_AD2_lowQ_py3.py new_gatk_files/{child}.recalibrated.bam new_gatk_files/{father}.recalibrated.bam new_gatk_files/{mother}.recalibrated.bam > dat_files/{caller}/{specie}/{ftype}_{child}{extra}_w_lowMQ.dat
'''

add_parents_cov = \
    template(
        input=['dat_files/{caller}/{specie}/{ftype}_{child}{extra}.dat'],
        output=['dat_files/{caller}/{specie}/{ftype}_{child}{extra}_w_pc.dat'],
        memory='4g', account='MutationRates', walltime='24:00:00') << \
'''
cat dat_files/{caller}/{specie}/{ftype}_{child}{extra}.dat | python python_scripts/add_individual_coverage.py > dat_files/{caller}/{specie}/{ftype}_{child}{extra}_w_pc.dat
'''

add_repeat = \
    template(
        input=['dat_files/{caller}/{specie}/{ftype}_{child}{extra}.dat'],
        output=['dat_files/{caller}/{specie}/{ftype}_{child}{extra}_w_repeat.dat'],
        memory='4g', account='MutationRates', walltime='24:00:00') << \
'''
cat dat_files/{caller}/{specie}/{ftype}_{child}{extra}.dat | python python_scripts/add_repeat_status.py ~/Data/2bit/{refgenome}.2bit > dat_files/{caller}/{specie}/{ftype}_{child}{extra}_w_repeat.dat
'''

add_segdup = \
    template(
	input=['dat_files/{caller}/{specie}/{ftype}_{child}{extra}.dat'],
    output=['dat_files/{caller}/{specie}/{ftype}_{child}{extra}_w_sd.dat'],
    memory='4g', account='MutationRates', walltime='24:00:00') << \
'''
cat dat_files/{caller}/{specie}/{ftype}_{child}{extra}.dat | head -1 | awk -vOFS="\t" '{{print $0,"SD"}}' > dat_files/{caller}/{specie}/\
{ftype}_{child}{extra}_w_sd.dat
cat dat_files/{caller}/{specie}/{ftype}_{child}{extra}.dat | awk '$1!="CHROM"' | gsort -k1,1 -k2,2n | ~/Scripts/bed_filter_pos.py segdup/{refgenome}_combined.bed >> dat_files/{caller}/{specie}/{ftype}_{child}{extra}_w_sd.dat
'''

