rule angsd_index_sites:
  input:
    'ngsLD/{sets}/LDpruned_snps_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.list'
  output:
    touch('ngsLD/{sets}/index_SNPs_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done')
  log:
    'log/{sets}/index_SNPs_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.log'
  threads: 12
  message:
    """ Index SNP list using angsd """
  shell:
    """
    module load angsd/0.938
    /apps/uibk/bin/sysconfcpus -n {threads} angsd sites index {input} 2> {log}
    """



rule getChrom_from_sites:
  input:
    'ngsLD/{sets}/LDpruned_snps_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.list',
    'ngsLD/{sets}/index_SNPs_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done'
  output:
    'ngsLD/{sets}/LDpruned_snps_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.chr'
  log:
    'log/{sets}/getChrom_from_LDpruned_snps_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.log'
  message:
    """ Get a list of chromosomes/scaffolds """
  shell:
    """
     cut -f1 {input} | sort | uniq > {output}
    """


#rule index_bam:
#  input:
#    'realigned/{sample}.bam'
#  output:
#    'realigned/{sample}.bam.bai'
#  log: 'log/{sample}.index_bam.log'
#  message: """ --- Indexing bam files --- """
#  shell:
#    """
#    samtools index -b {input} 2> {log}
#    """


rule getbeagle_LDpruned:
  input:
    ref = config["ref_rapid"],
    bamlist = 'depth/stats/{sets}_realignedBAM_df1.list',
    touched = 'ngsLD/{sets}/index_SNPs_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done',
    sites = 'ngsLD/{sets}/LDpruned_snps_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.list',
    chroms = 'ngsLD/{sets}/LDpruned_snps_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.chr'
  output:
    touch('angsd/{sets}/LDpruned_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done')
  log:
    'log/{sets}/LDpruned_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.log'
  threads: 24
  message:
    """ Calculate genotype likelihoods on predefined LD pruned sites only using angsd (-doMajorMinor 3 -sites)  """
  shell:
    """
    module load angsd/0.938
    /apps/uibk/bin/sysconfcpus -n 24 angsd -b {input.bamlist} -out angsd/{wildcards.sets}/LDpruned_angsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth} -anc {input.ref} -GL {wildcards.GL} -doGlf 2 -doMajorMinor 3 -doMaf 1 -sites {input.sites} -rf {input.chroms} 2> {log}
    """
