rule PCAngsd_LDpruned:
  input:
    touched = 'angsd/{sets}/LDpruned_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done'
  output:
    touch('pcangsd/{sets}/PCAngsd_LDpruned_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done')
  log: 'log/{sets}/PCAngsd_LDpruned_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.log'
  threads: 24
  message:
    """ Estimate covariance matrix from GL using PCAngsd """
  shell:
    """
    module load pcangsd/1.01

    for i in {{3..10}}; do
       pcangsd -beagle angsd/{wildcards.sets}/LDpruned_angsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth}.beagle.gz -o pcangsd/{wildcards.sets}/PCAngsd_LDpruned_angsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth} -minMaf {wildcards.minMaf} -admix -admix_K $i -e 10 -threads {threads} 2> {log}
    done
    """



#rule plot_covMat_LDpruned:
#  input:
#    touched = 'pcangsd/{sets}/PCAngsd_LDpruned_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done',
#    metadata = 'list/Da_Resteggs_metadata_cleaned_sorted_{sets}.tsv'
#  output:
#    pdf = 'pcangsd/{sets}/PCA_LDpruned_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.pdf'
#  log: 'log/{sets}/PCA_LDpruned_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_plotCovMat.log'
#  threads: 12
#  message:
#    """ Plot the estimated covariance matrix from PCAngsd in R """
#  shell:
#    """
#    Rscript scripts/plot_covMat_LDpruned.R pcangsd/{wildcards.sets}/PCAngsd_LDpruned_angsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth}.cov {input.metadata} {output.pdf} 2> {log}
#    """
#
#
#
#rule PCAngsd_admix_LD:
#  input:
#    touched = 'angsd/{sets}/LDpruned_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done'
#  output:
#    touch('pcangsd/{sets}/PCAngsd_LDpruned_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_admix_K{K}.done')
#  log: 'log/{sets}/PCAngsd_LDpruned_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_admix_K{K}.log'
#  threads: 24 
#  message:
#    """ Infer admixture proportions using PCAngsd """
#  shell:
#    """
#    module load pcangsd/1.01
#    pcangsd -beagle angsd/{wildcards.sets}/LDpruned_angsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth}.beagle.gz -o pcangsd/{wildcards.sets}/PCAngsd_LDpruned_angsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth}_admix_K{wildcards.K} -admix -admix_K {wildcards.K} -minMaf 0.05 -e 10 -threads {threads} 2> {log}
#    """
#
#
#
#rule PCAngsd_kinship:
#  input:
#    touched = 'angsd/{sets}/LDpruned_angsd_GL2_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done'
#  output:
#    touch('pcangsd/PCAngsd_LDpruned_angsd_GL2_minInd184_maf0.05_minDepth184_maxDepth8656_kinship.done')
#  log: 'log/PCAngsd_LDpruned_angsd_GL2_minInd184_maf0.05_minDepth184_maxDepth8656_kinship.log'
#  threads: 12 
#  message:
#    """ Infer admixture proportions using PCAngsd """
#  shell:
#    """
#    module load pcangsd/1.01
#    BEAGLE=(angsd/LDpruned_angsd_GL2_minInd184_maf0.05_minDepth184_maxDepth8656.beagle.gz)
#    pcangsd -beagle $BEAGLE -kinship -minMaf 0.05 -o pcangsd/PCAngsd_LDpruned_angsd_GL2_minInd184_maf0.05_minDepth184_maxDepth8656_kinship 2> {log}
#    """
#
#
#
#rule PCAngsd_inbreedingSites:
#  input:
#    touched = 'angsd/LDpruned_angsd_GL2_minInd184_maf0.05_minDepth184_maxDepth8656.done'
#  output:
#    touch('pcangsd/PCAngsd_LDpruned_angsd_GL2_minInd184_maf0.05_minDepth184_maxDepth8656_inbreeding.done')
#  log: 'log/PCAngsd_LDpruned_angsd_GL2_minInd184_maf0.05_minDepth184_maxDepth8656_inbreeding.log'
#  threads: 12 
#  message:
#    """ Infer admixture proportions using PCAngsd """
#  shell:
#    """
#    module load pcangsd/1.01
#    BEAGLE=(angsd/LDpruned_angsd_GL2_minInd184_maf0.05_minDepth184_maxDepth8656.beagle.gz)
#    pcangsd -beagle $BEAGLE -inbreedSites -minMaf 0.05 -o pcangsd/PCAngsd_LDpruned_angsd_GL2_minInd184_maf0.05_minDepth184_maxDepth8656_inbreeding -sites_save 2> {log}
#    """
#
#
#rule PCAngsd_selection:
#  input:
#    touched = 'angsd/LDpruned_angsd_GL2_minInd184_maf0.05_minDepth184_maxDepth8656.done'
#  output:
#    touch('pcangsd/PCAngsd_LDpruned_angsd_GL2_minInd184_maf0.05_minDepth184_maxDepth8656_kinship.done')
#  log: 'log/PCAngsd_LDpruned_angsd_GL2_minInd184_maf0.05_minDepth184_maxDepth8656_kinship.log'
#  threads: 12 
#  message:
#    """ Infer admixture proportions using PCAngsd """
#  shell:
#    """
#    module load pcangsd/1.01
#    BEAGLE=(angsd/LDpruned_angsd_GL2_minInd184_maf0.05_minDepth184_maxDepth8656.beagle.gz)
#    pcangsd -beagle $BEAGLE -selection -minMaf 0.05 -o pcangsd/PCAngsd_LDpruned_angsd_GL2_minInd184_maf0.05_minDepth184_maxDepth8656_selection -sites_save  2> {log}
#    """
#
#
