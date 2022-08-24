rule get_beagle_allSNPs:
  input:
    ref = config["ref_rapid"],
    bamlist = 'depth/stats/{sets}_realignedBAM_df1.list'
  output:
    touch('angsd/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done')
  log: 'log/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.log'
  threads: 24
  message:
    """ Calculate genotype likelihoods and call SNPs for all bam files using angsd setting various cutoffs """
  shell:
    """
    module load angsd/0.935
    /apps/uibk/bin/sysconfcpus -n 24 angsd -b {input.bamlist} -ref {input.ref} -out angsd/{wildcards.sets}/angsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth} -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minMapQ 30 -minQ 20 -skipTriallelic 1 -minInd {wildcards.IND} -setMinDepth {wildcards.MinDepth} -setMaxDepth {wildcards.MaxDepth} -GL {wildcards.GL} -doGlf 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minMaf {wildcards.minMaf} -doCounts 1 2> {log}
    """
