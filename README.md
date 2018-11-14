# Great Ape Mutation Rate workflow

This workflow was implemented using the old pre-1.0 version of gwf that has now been retired and is no longer available for download. I will keep this repository as it is to document the analysis that was used in the article. But I am working on simplifying the mutation rate calling workflow and updating it the current version of gwf, so if you are interested in running this workflow on a new data set you should send me an email asking about the progress of that updated version.

### Input files: ###
* Pedigree of samples
* a bam-file for each sample
* reference genome in fasta format and 2bit format

### Necessary software: ###
* samtools
* GATK
* Python, including the following packages:
    * bx-python
    * [PyVCF](http://pyvcf.readthedocs.io/en/latest/)
    * 
* R, including the following packages:
    * dplyr
    