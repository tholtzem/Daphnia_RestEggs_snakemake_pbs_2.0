
rule PCAngsd_covmat:
  input: 
    touched = 'angsd/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done' 
  output:
    touch('pcangsd/{sets}/PCAngsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done')
  log: 'log/{sets}/PCAngsd_covmat_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.log'
  threads: 24
  message:
    """ Estimate covariance matrix from GL using PCAngsd """
  shell:
    """
    module load pcangsd/1.01
    pcangsd -beagle angsd/{wildcards.sets}/angsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth}.beagle.gz -o pcangsd/{wildcards.sets}/PCAngsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth} -minMaf {wildcards.minMaf} -e 10 -threads {threads} 2> {log}
    """


#rule plot_covMat:
#  input:
#    #touched = 'pcangsd/PCAngsd_GL2_{IND}.done',
#    #covMat ='pcangsd/PCAngsd_GL2_{IND}_covmat.cov'
#    touched = 'pcangsd/{sets}/PCAngsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done',
#    #metadata = 'list/samples184_metadata.tsv'
#   metadata = 'list/Da_Resteggs_metadata_cleaned_sorted_{sets}.tsv'
# output:
#    pdf = 'pcangsd/{sets}/PCAngsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_covmat.pdf'
#  log: 'log/{sets}/PCAngsd_plotCovaMat_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.log'
#  threads: 12
#  message:
#    """ Estimate covariance matrix from GL using PCAngsd """
#  shell:
#    """
#    Rscript scripts/plot_covMat.R pcangsd/{wildcards.sets}/PCAngsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth}.cov {input.metadata} {output.pdf} 2> {log}
#    """

rule PCAngsd_admix_unpruned:
  input:
    touched = 'angsd/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done',
    BEAGLE = 'angsd/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.beagle.gz'
  output:
    touch('pcangsd/{sets}/PCAngsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_unpruned_admix_K{K}.done')
  log: 'log/{sets}/PCAngsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_unpruned_admix_K{K}.log'
  threads: 24 
  message:
    """ Infer admixture proportions using PCAngsd (unpruned data) """
  shell:
    """
    module load pcangsd/1.01 
    pcangsd -beagle {input.BEAGLE} -o pcangsd/{wildcards.sets}/PCAngsd_LDpruned_angsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth}_admix_K{wildcards.K} -admix -admix_K {wildcards.K} -minMaf 0.05 -e 10 -threads {threads} 2> {log}
    """

