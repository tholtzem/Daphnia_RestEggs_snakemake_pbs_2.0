
rule Index_realignedBAM:
  input:
    realigned = 'realigned/{sample}.realigned.bam'
  output:
    idx = 'realigned/{sample}.realigned.bai'
  log: 'log/{sample}.realigned.bai.log'
  threads: 2
  message: """--- Indexing realigned BAM files with samtools ---"""
  shell:
    """
    samtools index {input.realigned} {output.idx} 2> {log}
    """


rule get_bcf:
  input:
    ref = config["ref_HiC"],
    bamlist = 'list/ancestry/parents_hybrids_LC_realignedBAM.list',
    #bamlist = 'depth/stats/LC_realignedBAM_df1.list',
    sites = 'angsd/LC/{prefix}_globalSNP.list',
    chroms = 'angsd/LC/{prefix}_globalSNP.chr'
  output:
    touch('ancestry/LC/data/get_bcf_{prefix}.done')
  log: 'log/LC/get_bcf_{prefix}.log'
  threads: 36
  message:
    """ Call genotypes from realigned bam files for hybrids and parental species using estimated allele frequencies from genotype likelihoods as a prior and output in bcf format """
  shell:
    """
    module load angsd/0.938
    /apps/uibk/bin/sysconfcpus -n 360 angsd -b {input.bamlist} -ref {input.ref} -out ancestry/LC/data/{wildcards.prefix} -GL 2 -doPost 1 -doMajorMinor 3 -doMaf 1 -doBcf 1 --ignore-RG 0 -doGeno 1 -doCounts 1 -geno_minDepth 10 -sites {input.sites} -rf {input.chroms} 2> {log}
    """


rule bcf2vcf:
  input:
    'ancestry/LC/data/get_bcf_{prefix}.done'
  output:
    'ancestry/LC/data/{prefix}.vcf.gz'
  log:
    'log/LC/{prefix}_BCF2VCF.log'
  message: """ --- Convert bcf 2 uncompressed vcf for downstream analysis --- """
  threads: 12
  shell:
    """
    bcftools index -f --csi ancestry/LC/data/{wildcards.prefix}.bcf --threads {threads} &&
    bcftools convert -Oz -o {output} ancestry/LC/{wildcards.prefix}.bcf --threads {threads} 2> {log}
    """


rule stats_DP_vcf:
  input:
    vcf = 'ancestry/LC/data/{prefix}.vcf.gz'
  output:
    INFO_DP = 'ancestry/LC/data/{prefix}_INFODP.txt',
    n_sites = 'ancestry/LC/data/{prefix}_nbrSites.txt'
  log: 'log/{prefix}_statsVCF.log'
  threads: 12
  message: """ Count sites of raw vcf and for each site print DP values """
  shell:
    """
    bcftools view -H {input.vcf} | wc -l > {output.n_sites} &&
    bcftools query {input.vcf} -f'%DP\n' > {output.INFO_DP} 2> {log} 
    """


rule plot_INFODP:
  input:
    args1 = 'ancestry/LC/data/{prefix}_INFODP.txt',
  output:
    args2 = 'ancestry/LC/data/{prefix}_INFODP.pdf',
    args3 = 'ancestry/LC/data/{prefix}_INFODP.stats.txt'
  log: 'log/{prefix}_plot_INFODP.log'
  threads: 12
  message: """ Plot INFO/DP values """
  shell:
    """
    Rscript scripts/plot_INFODP.R {input.args1} {output.args2} {output.args3} 2> {log}
    """


rule hardfilter_vcf:
  input:
    stats = 'ancestry/LC/data/{prefix}_INFODP.stats.txt',
    vcf = 'ancestry/LC/data/{prefix}.vcf.gz'
  output:
    vcf = 'ancestry/LC/data/{prefix}_hf.vcf.gz',
    n_sites = 'ancestry/LC/data/{prefix}_hf.nbrSites.txt'
  log: 'log/{prefix}_hf.log'
  threads: 12
  message: """ Filter average site-level (maxDepth) """
  shell:
    """
    bcftools filter -i 'INFO/DP<2088' -Oz -o {output.vcf} {input.vcf} &&
    bcftools view -H {output.vcf} | wc -l > {output.n_sites} 2> {log}
    """


rule filter_genotype_DP:
  input:
    'ancestry/LC/data/{prefix}_hf.vcf.gz'
  output:
    vcf = 'ancestry/LC/data/{prefix}_hf.DP{genoDP}.vcf.gz',
    n_sites = 'ancestry/LC/data/{prefix}_hf.DP{genoDP}.nbrSites.txt'
  log: 'log/{prefix}_hf.DP{genoDP}.log'
  threads: 12
  message: """ --- Filter VCF for individual genotypes for read depth --- """
  shell:
    """
    bcftools filter -S . -e 'FMT/DP<10' -Oz -o {output.vcf} {input} &&
    bcftools view -H {output.vcf} | wc -l > {output.n_sites} 2> {log}
    """


rule imiss:
  input:
    'ancestry/LC/data/{prefix}_hf.DP{genoDP}.vcf.gz'
  output:
    'ancestry/LC/data/{prefix}_hf.DP{genoDP}.imiss'
  log: 'log/{prefix}_hf.DP{genoDP}.imiss.log'
  threads: 12
  message: """ --- Identify individuals with a high amount of missing data --- """
  shell:
    """
    bcftools stats -s - {input} | grep -E ^PSC | cut -f3,14 > {output} 2> {log}
    """


rule plot_imiss:
  input:
    imiss = 'ancestry/LC/data/{prefix}_hf.DP{genoDP}.imiss',
    n_sites = 'ancestry/LC/data/{prefix}_hf.DP{genoDP}.nbrSites.txt'
  output:
    pdf = 'ancestry/LC/data/{prefix}_hf.DP{genoDP}.imiss.pdf',
    tsv = 'ancestry/LC/data/{prefix}_hf.DP{genoDP}.imiss.tsv'
  log: 'log/plot_{prefix}_hf.DP{genoDP}.imiss.log'
  threads: 4
  message: """--- Identify individuals with a high amount of missing data ---"""
  shell:
    """
    N_SITES=$(cat {input.n_sites})
    Rscript scripts/plot_imiss.R {input.imiss} $N_SITES {output.pdf} {output.tsv} 2> {log}
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
    'ancestry/LC/data/{prefix}_hf.DP{genoDP}.vcf.gz'
  output:
    vcf = 'ancestry/LC/data/{prefix}_hf.DP{genoDP}_vmiss20.vcf',
    n_sites = 'ancestry/LC/data/{prefix}_hf.DP{genoDP}_nbrSites.vmiss20.txt'
  log: 'log/{prefix}_hf.DP{genoDP}_vmiss20.log'
  message: """ Remove sites with more than 20 % missingness, output uncompressed vcf for downstream analysis """
  shell:
    """
    bcftools filter -e 'F_MISSING > 0.2' -Ov -o {output.vcf} {input} &&
    bcftools view -H {input} | wc -l > {output.n_sites} 2> {log}
    """
