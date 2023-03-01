
rule Index_realignedBAM:
  input:
    realigned = 'realigned/{sample}.realigned.bam'
  output:
    idx = 'realigned/{sample}.realigned.bai'
  log: 'log/{sample}.realigned.bam.bailog'
  threads: 2
  message: """--- Indexing realigned BAM files with samtools ---"""
  shell:
    """
    samtools index {input.realigned} {output.idx} 2> {log}
    """


rule get_bcf:
  input:
    ref = config["ref_rapid"],
    bamlist = 'depth/stats/LC_REF_realignedBAM_df1.list',
    sites = 'angsd/LC_REF/angsd_GL2_minInd221_maf0.05_minDepth221_maxDepth10006_globalSNP.list',
    chroms = 'angsd/LC_REF/angsd_GL2_minInd221_maf0.05_minDepth221_maxDepth10006_globalSNP.chr'
  output:
    touch('ancestry/get_bcf.done')
  log: 'log/get_bcf.log'
  threads: 24
  message:
    """ Call genotypes using estimated allele frequencies from genotype likelihoods as a prior, output in bcf format """
  shell:
    """
    module load angsd/0.938
    /apps/uibk/bin/sysconfcpus -n 240 angsd -b {input.bamlist} -ref {input.ref} -out ancestry/angsd_GL2_minInd221_maf0.05_minDepth221_maxDepth10006 -GL 2 -doPost 1 -doMajorMinor 3 -doMaf 1 -doBcf 1 --ignore-RG 0 -doGeno 1 -doCounts 1 -geno_minDepth 10 -sites {input.sites} -rf {input.chroms} 2> {log}
    """


rule bcf2vcf:
  input:
    touched = 'ancestry/get_bcf.done'
  output:
    'ancestry/angsd_GL2_minInd221_maf0.05_minDepth221_maxDepth10006.vcf.gz'
  log:
    'log/angsd_GL2_minInd221_maf0.05_minDepth221_maxDepth10006_BCF2VCF.log'
  message: """ --- Convert bcf 2 uncompressed vcf for downstream analysis --- """
  threads: 12
  shell:
    """
    BCF=(ancestry/angsd_GL2_minInd221_maf0.05_minDepth221_maxDepth10006.bcf)
    bcftools index -f --csi $BCF --threads {threads} &&
    bcftools convert --threads {threads} -Oz -o {output} $BCF 2> {log}
    """


rule list_parents_hybrids:
  output:
    'list/parents_hybrids.list'
  log: 'log/parents_hybrids_list.log'
  threads: 12
  message: """--- Create a list of parental species and hybrids, one sample per line ---"""
  shell:
    """
    path=(list/pop_list)
    hybrids=$(cat $path/hybrids_LC.txt)
    LONG=$(cat $path/longispina.txt | cut -f1 -d'.' | cut -f2 -d'/')
    GAL=$(cat $path/galeata.txt | cut -f1 -d'.' | cut -f2 -d'/')
    CUC=$(cat $path/cucullata.txt | cut -f1 -d'.' | cut -f2 -d'/')
    echo $hybrids $LONG $GAL $CUC | sed -z 's/ /\n/g;s/,$/\n/' > {output} 2> {log}
    """


rule subset_vcf:
  input:
    samples = 'list/parents_hybrids.list',
    vcf = 'ancestry/angsd_GL2_minInd221_maf0.05_minDepth221_maxDepth10006.vcf.gz'
  output:
    vcf = 'ancestry/parents_hybrids.vcf.gz'
  log: 'log/parents_hybrids.log'
  threads: 12
  message: """--- Subset vcf using a file with a list of samples---"""
  shell:
    """
    bcftools view -S {input.samples} {input.vcf} -Oz -o {output.vcf} 2> {log}
    """


#rule filter_QUAL_DP:
#  input:
#    vcf = 'ancestry/parents_hybrids.vcf.gz'
#  output:
#    'ancestry/parents_hybrids_Q30_DP10.vcf.gz'
#  log: 'log/parents_hybrids_Q30_DP10.log'
#  threads: 4
#  message: """--- Filter VCF for site quality and individual genotypes for read depth and genotype quality ---"""
#  shell:
#    """
#    bcftools view -e 'QUAL<30' {input} | bcftools filter -S . -e 'FMT/DP<10' -O z -o {output} 2> {log}
#    """
#
#
#rule imiss:
#  input:
#    'ancestry/parents_hybrids_Q30_DP10.vcf.gz'
#  output:
#    'ancestry/parents_hybrids_Q30_DP10_alleles.imiss'
# log: 'log/parents_hybrids_Q30_DP10_alleles.imiss.log'
# threads: 4
#  message: """--- Identify individuals with a high amount of missing data ---"""
#  shell:
#    """
#    bcftools stats -s - {input} | grep -E ^PSC | cut -f3,14 > {output} 2> {log}
#    """
##, output unzipped vcf for downstream analysis output unzipped vcf for downstream analysiss 
