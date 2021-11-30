rule angsd_index_sites:
  input:
    'list/LDpruned_snps.list'
    #'list/LDpruned_snps_dgal1.list'
  output:
    touch('list/index_SNP.done')
    #touch('list/index_snp_dgal1.done')
  log:
    'log/index_SNP.done'
    #'log/index_snp_dgal1.log'
  threads: 12
  message:
    """ Index SNP list using angsd """
  shell:
    """
    angsd sites index {input} 2> {log}
    """

rule getChrom_from_sites:
  input:
    'list/LDpruned_snps.list'
    #'list/LDpruned_snps_dgal1.list'
  output: 
    touch('list/getChrom_from_sites.done')
    #touch('list/getChrom_dgal1_from_sites.done')
  log:
    'log/getChrom_from_sites.done'
    #'log/getdgal1_from_sites.done'
  message:
    """ Get the chromosomes/scaffolds for which we have sites we want to analyse """
  shell:
    """
     cut -f1 {input} | sort | uniq > list/LDpruned_Chrom.list
    """

rule index_bam:
  input:
    'realigned/{sample}.bam'
  output:
    'realigned/{sample}.bam.bai'
  log: 'log/{sample}.index_bam.log'
  message: """ --- Indexing bam files --- """
  shell:
    """
    samtools index -b {input} 2> {log}
    """

rule angsd_GL2_LDpruned:
  input:
    ref = config["ref"],
    #touched = 'list/index_snp_dgal1.done',
    #touched2 = 'list/getChrom_dgal1_from_sites.done'
    touched = 'list/index_SNP.done',
    touched2 = 'list/getChrom_from_sites.done'
  output:
    touch('angsd/angsd_LC_GL2_maf0018_LDpruned.done')
    #touch('angsd/angsd_LC_GL2_maf0018_LDpruned_dgal1.done')
  log:
    'log/angsd_LC_GL2_maf0018_LDpruned.log'
    #'log/angsd_LC_GL2_maf0018_LDpruned_dgal1.log'
  threads: 24
  message:
    """ Calculate genotype likelihoods on predefined LD pruned sites only using angsd (-doMajorMinor 3 -sites)  """
  shell:
    """
    angsd -b 'list/realignedBAM_depth10.list' -anc {input.ref} -out angsd/angsd_LC_GL2_maf0018_LDpruned -GL 2 -doGlf 2 -doMajorMinor 3 -doMaf 1 -doCounts 1 -sites 'list/LDpruned_snps.list' -rf list/LDpruned_Chrom.list 2> {log}
    """

#module load angsd-0.921-gcc-4.8-b7bamor


rule PCAngsd_LDpruned:
  input:
    touched = 'angsd/angsd_LC_GL2_maf0018_LDpruned.done'
    #GL2 = 'angsd/angsd_LC_GL2_LDpruned.beagle.gz'
  output:
    touch('pcangsd/PCAngsd_GL2_maf0018_LDpruned.done')
  log: 'log/PCAngsd_GL2_maf0018_LDpruned_covmat.log'
  threads: 12
  message:
    """ Estimate covariance matrix from GL using PCAngsd """
  shell:
    """
    module load pcangsd/1.01
    pcangsd -beagle angsd/angsd_LC_GL2_maf0018_LDpruned.beagle.gz -o pcangsd/PCAngsd_GL2_maf0018_LDpruned_covmat 2> {log}
    """

rule plot_covMat_LDpruned:
  input:
    touched = 'pcangsd/PCAngsd_GL2_maf0018_LDpruned.done'
    #covMat = 'pcangsd/PCAngsd_GL2_LDpruned_covmat.cov'
  output:
    pdf = 'pcangsd/PCAngsd_GL2_maf0018_LDpruned.pdf'
  log: 'log/PCAngsd_GL2_LDpruned_plotCovMat.log'
  threads: 12
  message:
    """ Estimate covariance matrix from GL using PCAngsd """
  shell:
    """
    Rscript scripts/plot_covMat_LDpruned.R pcangsd/PCAngsd_GL2_maf0018_LDpruned_covmat.cov {output.pdf} 2> {log}
    """
