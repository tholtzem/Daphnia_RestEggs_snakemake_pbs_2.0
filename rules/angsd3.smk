rule angsd_GL3:
  input:
    ref = config["ref"]
  output:
    touch('angsd/angsd_GL3_cutoffs_55_maf0018.done')
  log: 'log/angsdGL3_cutoffs_55_maf0018.log'
  threads: 24
  message:
    """ Calculate genotype likelihoods for bam files having a read depth gretaer than 10 and try varying the cutoff -nInd and record how many sites are predicted to be variable for each scenario """
  shell:
    """
    angsd -b 'list/realignedBAM_depth10.list' -ref {input.ref} -out angsd/angsd_LC_GL3_cutoff.nInd55_maf0018 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -minMapQ 20 -minQ 20 -minInd 55 -setMinDepth 400 -setMaxDepth 2000 -GL 2 -doGlf 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -doCounts 1 2> {log}
    """

