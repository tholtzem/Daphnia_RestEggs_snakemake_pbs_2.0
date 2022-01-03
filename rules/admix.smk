rule PCAngsd_admix:
  input:
    'angsd/angsd_LC_GL2_maf0018_LDpruned.done'
  output:
    touch('pcangsd/PCAngsd_GL2_maf0018_LDpruned_admix_K3.done')
  log: 'log/PCAngsd_GL2_maf0018_LDpruned_admix_K3.log'
  threads: 12 
  message:
    """ Infer admixture proportions using PCAngsd """
  shell:
    """
    module load pcangsd/1.01
    pcangsd -beagle angsd/angsd_LC_GL2_maf0018_LDpruned.beagle.gz -admix -admix_K 3 -o pcangsd/PCAngsd_GL2_maf0018_LDpruned_admix_K3 2> {log}
    """

