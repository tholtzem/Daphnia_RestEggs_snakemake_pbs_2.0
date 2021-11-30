rule angsd_GL2:
  input:
    ref = config["ref"]
  output:
    #touch('angsd/GL2_cutoffs_{IND}.done')
    touch('angsd/angsd_GL2_cutoffs_maf0018_nInd55.done')
  log: 'log/angsd_GL2_cutoffs_maf0018_nInd55.log'
  threads: 24
  message:
    """ Calculate genotype likelihoods for bam files having a read depth gretaer than 10 and the cutoff -nInd 55 and -minMaf 0.018 (1.8 %) while recording how many sites are predicted to be variable for each scenario, """
  shell:
    """
    angsd -b 'list/realignedBAM_depth10.list' -ref {input.ref} -out angsd/angsd_LC_GL2_cutoff_maf0018_nInd55 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -minMapQ 20 -minQ 20 -minInd 55 -setMinDepth 400 -setMaxDepth 2000 -GL 2 -doGlf 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minMaf 0.018 -doCounts 1 2> {log}
    """

