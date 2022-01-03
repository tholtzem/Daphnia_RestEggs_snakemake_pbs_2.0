
rule PCAngsd2:
  input: 
    touched = 'angsd/angsd_GL2_cutoffs_maf0018_nInd55.done'
  output:
    touch('pcangsd/PCAngsd_GL2_cutoffs_maf0018_nInd55.done')
  log: 'log/PCAngsd_GL2_cutoffs_maf0018_nInd55_covmat.log'
  threads: 12
  message:
    """ Estimate covariance matrix from GL using PCAngsd """
  shell:
    """
    module load pcangsd/1.01
    pcangsd -beagle angsd/angsd_LC_GL2_cutoff_maf0018_nInd55.beagle.gz -o pcangsd/PCAngsd_GL2_cutoffs_maf0018_nInd55_covmat 2> {log}
    """

rule plot_covMat:
  input:
    #touched = 'pcangsd/PCAngsd_GL2_{IND}.done',
    #covMat ='pcangsd/PCAngsd_GL2_{IND}_covmat.cov'
    touched = 'pcangsd/PCAngsd_GL2_cutoffs_maf0018_nInd55.done' 
  output:
    pdf = 'pcangsd/PCAngsd_GL2_cutoffs_maf0018_nInd55_covmat.pdf'
  log: 'log/PCAngsd_GL2_cutoffs_maf0018_plotCovMat_nInd55.log'
  threads: 12
  message:
    """ Estimate covariance matrix from GL using PCAngsd """
  shell:
    """
    Rscript scripts/plot_covMat.R pcangsd/PCAngsd_GL2_cutoffs_maf0018_nInd55_covmat.cov {output.pdf} 2> {log}
    """

