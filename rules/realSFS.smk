# Note: The estimation of the SFS is a 2-step procedure: first, we generate a ".saf.idx" file (site allele frequency likelihood) using angsd (-anc -doSaf 1), followed by an optimization of the .saf.idx file which will estimate the Site frequency spectrum (SFS) ".sfs" using realSFS.  Folding should now be done in realSFS and not in the saf file generation (-fold 1). 

rule angsd_saf_samples:
  input:
    ref = config["ref_rapid"],
    bam = 'realigned/{sample}.realigned.bam'
  output:
    touch('saf/saf_samples/{sample}.GL2.saf.idx.done')
  log:
    'log/saf_samples/{sample}.GL2.saf.idx.log'
  threads: 12
  message:
    """ Compute site allele frequency likelihood (.saf.idx) for single samples using angsd (required for the global heterozygosity estimate for single samples)"""
  shell:
    """
    module load angsd/0.938
    /apps/uibk/bin/sysconfcpus -n 12 angsd -i {input.bam} -out saf/saf_samples/{wildcards.sample}.GL2 -anc {input.ref} -ref {input.ref} -doSaf 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minQ 20 -minMapQ 30 -GL 2 2> {log}
    """


rule real_SFS_samples:
  input:
    'saf/saf_samples/{sample}.GL2.saf.idx.done'
  output:
    'saf/saf_samples/{sample}.GL2.est.ml'
  log:
    'log/saf_samples/{sample}.GL2.est.ml.log'
  threads: 12
  message:
    """ Optimize .saf.idx and estimate the site frequency spectrum (SFS) using realSFS and fold (required for the global heterozygosity estimate for single samples) """
  shell:
    """
    module load angsd/0.938
    /apps/uibk/bin/sysconfcpus -n 12 realSFS saf/saf_samples/{wildcards.sample}.GL2.saf.idx -fold 1 > {output} 2> {log}
    """

rule angsd_SAF_PRE_LONG:
  input:
    ref = config["ref_rapid"],
    bamlist = 'list/pop_list/df_PRE_long.txt'
  output:
    touch('saf/POPS/df_PRE_long.GL2.saf.idx.done')
  log:
   'log/saf_POPS/df_PRE_long.GL2.saf.idx.log'
  threads: 12
  message:
    """ Compute site allele frequency likelihood (.saf.idx) for each population using angsd """
  shell:
    """
    module load angsd/0.938
    /apps/uibk/bin/sysconfcpus -n 12 angsd -b {input.bamlist} -ref {input.ref} -anc {input.ref} -out saf/POPS/df_PRE_long.GL2 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minQ 20 -minMapQ 30 -minInd 8 -setMinDepth 8 -setMaxDepth 277 -GL 2 -doCounts 1 -doSaf 1 2> {log}
    """

rule angsd_SAF_df_POST_long_ALL:
  input:
    ref = config["ref_rapid"],
    bamlist = 'list/pop_list/df_POST_long_ALL.txt'
  output:
    touch('saf/POPS/df_POST_long_ALL.GL2.saf.idx.done')
  log:
    'log/saf_POPS/df_POST_long_ALL.GL2.saf.idx.log'
  threads: 12
  message:
    """ Compute site allele frequency likelihood (.saf.idx) for each population using angsd """
  shell:
    """
    module load angsd/0.938
    /apps/uibk/bin/sysconfcpus -n 12 angsd -b {input.bamlist} -ref {input.ref} -anc {input.ref} -out saf/POPS/df_POST_long_ALL.GL2 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minQ 20 -minMapQ 30 -minInd 36 -setMinDepth 36 -setMaxDepth 1250 -GL 2 -doCounts 1 -doSaf 1 2> {log}
   """

rule angsd_SAF_df_POST_long_LC:
  input:
    ref = config["ref_rapid"],
    bamlist = 'list/pop_list/df_POST_long_LC.txt'
  output:
   touch('saf/POPS/df_POST_long_LC.GL2.saf.idx.done')
  log:
    'log/saf_POPS/df_POST_long_LC.GL2.saf.idx.log'
  threads: 12
  message:
    """ Compute site allele frequency likelihood (.saf.idx) for each population using angsd """
  shell:
    """
    module load angsd/0.938
    /apps/uibk/bin/sysconfcpus -n 12 angsd -b {input.bamlist} -ref {input.ref} -anc {input.ref} -out saf/POPS/df_POST_long_LC.GL2 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minQ 20 -minMapQ 30 -minInd 21 -setMinDepth 21 -setMaxDepth 729 -GL 2 -doCounts 1 -doSaf 1 2> {log}
    """


rule onedsfs_folded:
  input:
    'saf/POPS/{POPS}.GL2.saf.idx.done',
    safidx = 'saf/POPS/{POPS}.GL2.saf.idx'
  output:
    sfs = 'saf/POPS/{POPS}.GL2.sfs'
  log:
    'log/saf_POPS/{POPS}.sfs.log'
  threads: 12
  message:
    """ Optimize .saf.idx and calculate all pairwise 2dsfs's (joint site frequency spectra) using realSFS and fold for theta estimates"""
  shell:
    """
    module load angsd/0.938
    /apps/uibk/bin/sysconfcpus -n 12 realSFS {input.safidx} -fold 1 > {output.sfs} 2> {log}
    """

rule sfs2theta:
  input:
    safidx = 'saf/POPS/{POPS}.GL2.saf.idx',
    sfs = 'saf/POPS/{POPS}.GL2.sfs'
  output:
    'saf/POPS/{POPS}.sfs2theta.done'
  log:
    'log/saf_POPS/{POPS}.sfs2theta.log'
  threads: 12
  message:
    """ calculate theta estimates for each site """
  shell:
    """
    module load angsd/0.938
    /apps/uibk/bin/sysconfcpus -n 12 realSFS saf2theta {input.safidx} -sfs {input.sfs} -outname saf/POPS/{wildcards.POPS} > {output} 2> {log}
    """


rule theta_stat_chrom:
  input:
    'saf/POPS/{POPS}.sfs2theta.done'
  output:
    touch('saf/POPS/{POPS}.theta_stats_chrom.done')
  log:
    'log/saf_POPS/{POPS}.theta_stats_chrom.log'
  threads: 12
  message:
    """ Estimate Tajimas D and other statistics """
  shell:
    """
    module load angsd/0.938
    /apps/uibk/bin/sysconfcpus -n 12 thetaStat do_stat {wildcards.POPS}.thetas.idx 2> {log}
    """

rule theta_stat_SW:
  input:
    'saf/POPS/{POPS}.sfs2theta.done'
  output:
    touch('saf/POPS/{POPS}.theta_stats_SW.done')
  log:
    'log/saf_POPS/{POPS}.theta_stats_SW.log'
  threads: 12
  message:
    """ Estimate Tajimas D and other statistics using a sliding window analysis """
  shell:
    """
    module load angsd/0.938
    /apps/uibk/bin/sysconfcpus -n 12 thetaStat do_stat {wildcards.POPS}.thetas.idx -win 10000 -step 1000 -outnames {wildcards.POPS}.theta.thetasWindow.gz 2> {log}
    """

rule twodsfs_folded:
  input:
    'saf/POPS/df_PRE_long.GL2.saf.idx.done',
    'saf/POPS/df_POST_long_ALL.GL2.saf.idx.done',
    'saf/POPS/df_POST_long_LC.GL2.saf.idx.done',
    pop1 = 'saf/POPS/df_PRE_long.GL2.saf.idx',
    pop2 = 'saf/POPS/df_POST_long_ALL.GL2.saf.idx',
    pop3 = 'saf/POPS/df_POST_long_LC.GL2.saf.idx'
  output:
    pop1_pop2 = 'saf/POPS/PRElong_ALLlong.sfs',
    pop1_pop3 = 'saf/POPS/PRElong_LClong.sfs'
  log:
    pop1_pop2 = 'log/saf_POPS/PRElong_ALLlong.sfs.log',
    pop1_pop3 = 'log/saf_POPS/PRElong_LClong.sfs.log'
  threads: 12
  message:
    """ Optimize .saf.idx and calculate all pairwise 2dsfs's (joint site frequency spectra) using realSFS and fold for FST """
  shell:
    """
    module load angsd/0.938
    /apps/uibk/bin/sysconfcpus -n 12 realSFS {input.pop1} {input.pop2} -fold 1 > {output.pop1_pop2} 2> {log.pop1_pop2} &&
    /apps/uibk/bin/sysconfcpus -n 12 realSFS {input.pop1} {input.pop3} -fold 1 > {output.pop1_pop3} 2> {log.pop1_pop3} 
    """


rule Fst_index:
  input:
    pop1_pop2 = 'saf/POPS/PRElong_ALLlong.sfs',
    pop1_pop3 = 'saf/POPS/PRElong_LClong.sfs'
  output:
    pop1_pop2 = touch('saf/POPS/PRElong_ALLlong.FstIndex.done'),
    pop1_pop3 = touch('saf/POPS/PRElong_LClong.FstIndex.done')
  log:
    pop1_pop2 = 'log/saf_POPS/PRElong_ALLlong.FstIndex.log',
    pop1_pop3 = 'log/saf_POPS/PRElong_LClong.FstIndex.log'
  threads: 12
  message:
    """ Index fst for easy window analysis """
  shell:
    """
    module load angsd/0.938
    /apps/uibk/bin/sysconfcpus -n 12 realSFS fst index saf/POPS/df_PRE_long.GL2.saf.idx saf/POPS/df_POST_long_ALL.GL2.saf.idx -sfs {input.pop1_pop2} -fstout saf/POPS/PRElong_ALLlong 2> {log.pop1_pop2} &&
    /apps/uibk/bin/sysconfcpus -n 12 realSFS fst index saf/POPS/df_PRE_long.GL2.saf.idx saf/POPS/df_POST_long_LC.GL2.saf.idx -sfs {input.pop1_pop3} -fstout saf/POPS/PRElong_LClong 2> {log.pop1_pop3}
    """

rule Fst_global:
  input:
    pop1_pop2 = 'saf/POPS/PRElong_ALLlong.FstIndex.done',
    pop1_pop3 = 'saf/POPS/PRElong_LClong.FstIndex.done'
  output:
    pop1_pop2 = touch('saf/POPS/PRElong_ALLlong.Fst_Global.done'),
    pop1_pop3 = touch('saf/POPS/PRElong_LClong.Fst_Global.done')
  log:
    pop1_pop2 = 'log/saf_POPS/PRElong_ALLlong.Fst_Global.log',
    pop1_pop3 = 'log/saf_POPS/PRElong_LClong.Fst_Global.log'
  threads: 12
  message:
    """ Get the global fst estimate """
  shell:
    """
    module load angsd/0.938
    /apps/uibk/bin/sysconfcpus -n 12 realSFS fst stats saf/POPS/PRElong_ALLlong.fst.idx 2> {log.pop1_pop2} &&
    /apps/uibk/bin/sysconfcpus -n 12 realSFS fst stats saf/POPS/PRElong_LClong.fst.idx 2> {log.pop1_pop3}
    """



#rule plotHeterozygosity_samples:
#  input:
#    'saf/{sample}.est.ml'
#  output:
#    'saf/{sample}.heterozygosity.txt'
#  log:
#    'log/{sample}.heterozygosity.log'
#  message:
#    """ Plot heterozygosity of single samples """
#  shell:
#    """
#   Rscript scripts/plotHeterozygosity.R {input} {output} 2> {log}
#    """



#
#
#rule plotHeterozygosity_pops:
#  input:
#    'saf/{pops}.saf.est.ml'
#  output:
#    'saf/{pops}.heterozygosity.txt'
#  log:
#    'log/{pops}.heterozygosity.log'
#  message:
#    """ Plot heterozygosity """
#  shell:
#    """
#   Rscript scripts/plotHeterozygosity.R {input} {output} 2> {log}
#    """



