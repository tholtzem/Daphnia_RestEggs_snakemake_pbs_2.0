
#rule Index_realignedBAM:
#  input:
#    realigned = 'realigned/{sample}.realigned.bam'
#  output:
#    idx = 'realigned/{sample}.realigned.bai'
#  log: 'log/{sample}.realigned.bam.bailog'
#  threads: 2
#  message: """--- Indexing realigned BAM files with samtools ---"""
#  shell:
#    """
#    samtools index {input.realigned} {output.idx} 2> {log}
#    """


rule get_bcf:
  input:
    ref = config["ref_HiC"],
    bamlist = 'list/ancestry/parents_hybrids_LC_realignedBAM.list',
    sites = 'angsd/LC/{prefix}_globalSNP.list',
    chroms = 'angsd/LC/{prefix}_globalSNP.chr'
  output:
    touch('ancestry/LC/get_bcf_{prefix}.done')
  log: 'log/LC/get_bcf_{prefix}.log'
  threads: 24
  message:
    """ Call genotypes from realigned bam files for hybrids and parental species using estimated allele frequencies from genotype likelihoods as a prior and output in bcf format """
  shell:
    """
    module load angsd/0.938
    /apps/uibk/bin/sysconfcpus -n 240 angsd -b {input.bamlist} -ref {input.ref} -out ancestry/LC/{wildcards.prefix} -GL 2 -doPost 1 -doMajorMinor 3 -doMaf 1 -doBcf 1 --ignore-RG 0 -doGeno 1 -doCounts 1 -geno_minDepth 10 -sites {input.sites} -rf {input.chroms} 2> {log}
    """


rule bcf2vcf:
  input:
    'ancestry/LC/get_bcf_{prefix}.done'
  output:
    'ancestry/LC/{prefix}.vcf.gz'
  log:
    'log/LC/{prefix}_BCF2VCF.log'
  message: """ --- Convert bcf 2 uncompressed vcf for downstream analysis --- """
  threads: 12
  shell:
    """
    bcftools index -f --csi {wildcards.prefix}.bcf --threads {threads} &&
    bcftools convert --threads {threads} -Oz -o {output} {wildcards.prefix}.bcf 2> {log}
    """


#rule list_parents_hybrids:
#  output:
#    'list/parents_hybrids.list'
#  log: 'log/parents_hybrids_list.log'
#  threads: 12
#  message: """--- Create a list of parental species and hybrids, one sample per line ---"""
#  shell:
#    """
#    path=(list/pop_list)
#   hybrids=$(cat $path/hybrids_LC.txt)
#    LONG=$(cat $path/longispina.txt | cut -f1 -d'.' | cut -f2 -d'/')
#    GAL=$(cat $path/galeata.txt | cut -f1 -d'.' | cut -f2 -d'/')
#    CUC=$(cat $path/cucullata.txt | cut -f1 -d'.' | cut -f2 -d'/')
#    echo $hybrids $LONG $GAL $CUC | sed -z 's/ /\n/g;s/,$/\n/' > {output} 2> {log}
#    """


#rule subset_vcf:
#  input:
#    samples = 'list/parents_hybrids.list',
#    vcf = 'ancestry/LC/{prefix}.vcf.gz'
#  output:
#    vcf = 'ancestry/parents_hybrids.vcf.gz'
#  log: 'log/parents_hybrids.log'
#  threads: 12
#  message: """--- Subset vcf using a file with a list of samples---"""
#  shell:
#    """
#    bcftools view -S {input.samples} {input.vcf} -Oz -o {output.vcf} 2> {log}
#    """

rule stats_DP_vcf:
  input:
    'ancestry/LC/{prefix}.vcf.gz'
  output:
    INFO_DP = 'ancestry/LC/{prefix}_INFODP.txt',
    n_sites = 'ancestry/LC/{prefix}_nbrSites.txt'
  log: 'log/{prefix}_statsVCF.log'
  threads: 12
  message: """ Count sites of raw vcf and for each site print DP values """
  shell:
    """
    bcftools view -H {input} | wc -l > {output.n_sites} &&
    bcftools query {input} -f'%DP\n' > {output.INFO_DP} 2> {log} 
    """


rule hardfilter_vcf:
  input:
    'ancestry/LC/{prefix}.vcf.gz'
  output:
    vcf = 'ancestry/LC/{prefix}_hf.vcf.gz',
    n_sites = 'ancestry/LC/{prefix}_hf_nbrSites.txt'
  log: 'log/{prefix}_hf.log'
  threads: 12
  message: """ Filter average site-level (maxDepth) """
  shell:
    """
    bcftools filter -i 'INFO/DP<7530' -Oz -o {output} {input} --threads {threads} &&
    bcftools view -H {output.vcf} | wc -l > {output.n_sites} 2> {log}
    """


rule filter_genotype_DP:
  input:
    'ancestry/LC/{prefix}.vcf.gz'
  output:
    vcf = 'ancestry/LC/{prefix}_hf_DP10.vcf.gz',
    n_sites = 'ancestry/LC/{prefix}_hf_DP10_nbrSites.txt'
  log: 'log/{prefix}_hf_DP10.log'
  threads: 12
  message: """ --- Filter VCF for individual genotypes for read depth --- """
  shell:
    """
    bcftools filter -S . -e 'FMT/DP<10' -Oz -o {output} {input} &&
    bcftools view -H {output.vcf} | wc -l > {output.n_sites} 2> {log}
    """


rule imiss:
  input:
    'ancestry/LC/{prefix}_hf_DP10.vcf.gz',
  output:
    'ancestry/LC/{prefix}_hf_DP10.imiss'
  log: 'log/{prefix}_hf_DP10.imiss.log'
  threads: 12
  message: """ --- Identify individuals with a high amount of missing data --- """
  shell:
    """
    bcftools stats -s - {input} | grep -E ^PSC | cut -f3,14 > {output} 2> {log}
    """


rule plot_imiss:
  input:
    imiss = 'ancestry/LC/{prefix}_hf_DP10.imiss',
    n_sites = 'ancestry/LC/{prefix}_hf_DP10_nbrSites.txt'
  output:
    'ancestry/LC/{prefix}_hf_DP10.imiss.pdf'
  log: 'log/plot_{prefix}_hf_DP10.imiss.log'
  threads: 4
  message: """--- Identify individuals with a high amount of missing data ---"""
  shell:
    """
    N_SITES=$(cat {input.n_sites})
    Rscript scripts/plot_imiss.R {input.imiss} {input.n_sites} {output} 2> {log}
    """

#rule remove_imiss:
#  input:
#    'ancestry/LC/{prefix}_hf_DP10.vcf.gz',
#  output:
#    'ancestry/LC/{prefix}_hf_DP10.imissRM.vcf.gz'
#  log: 'log{prefix}_hf_DP10.imissRM.imissRM.log'
#  message: """ Remove individual with high missingness from data set """
#  shell:
#    """
#    bcftools view -s ^FP2_1 -O z -o {output} {input}
#    """


rule lmiss:
  input:
    'ancestry/LC/{prefix}_hf_DP10.vcf.gz'
    #'ancestry/LC/{prefix}_hf_DP10.imissRM.vcf.gz'
  output:
    vcf = 'ancestry/LC/{prefix}_hf_DP10_lmiss20.vcf.gz',
    n_sites = 'ancestry/LC/{prefix}_hf_DP10_lmiss20_nbrSites.txt'
    #'ancestry/LC/{prefix}_hf_DP10.imissRM_lmiss20.vcf.gz'
  log: 'log/{prefix}_hf_DP10.imissRM_lmiss20.log'
  message: """ Remove sites with more than 20 % missingness """
  shell:
    """
    bcftools filter -e 'F_MISSING > 0.2' -Oz -o {output.vcf} {input} &&
    bcftools view -H {input} | wc -l > {output.n_sites} 2> {log}
    """




    
