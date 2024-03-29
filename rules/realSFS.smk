# Note: The estimation of the SFS is a 2-step procedure: first, we generate a ".saf.idx" file (site allele frequency likelihood) using angsd (-anc -doSaf 1), followed by an optimization of the .saf.idx file which will estimate the Site frequency spectrum (SFS) ".sfs" using realSFS.  Folding should now be done in realSFS and not in the saf file generation (-fold 1).

rule get_SAF_POP_list:
  input:
    'depth/stats/{POPS}_realignedBAM_df1.list'
  output:
    'list/saf/{POPS}.list'
  log:
   'log/saf/get_SAF_{POPS}_list.log'
  threads: 12
  message:
    """ copy population files for saf """
  shell:
    """
    cp {input} {output} 2> {log}
    """


rule angsd_SAF_POP:
  input:
    ref = config["ref_HiC"],
    bamlist = 'list/saf/{POPS}.list'
  output:
    touch('saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.saf.idx.done')
  log:
   'log/saf_POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.saf.idx.log'
  threads: 12
  message:
    """ Compute site allele frequency likelihood (.saf.idx) for each population using angsd """
  shell:
    """
    module load angsd/0.938
    /apps/uibk/bin/sysconfcpus -n 12 angsd -b {input.bamlist} -ref {input.ref} -anc {input.ref} -out saf/POPS/{wildcards.POPS}_{wildcards.IND}_{wildcards.MinDepth}_{wildcards.MaxDepth}.GL2 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minQ 20 -minMapQ 30 -minInd {wildcards.IND} -setMinDepth {wildcards.MinDepth} -setMaxDepth {wildcards.MaxDepth} -GL 2 -doCounts 1 -doSaf 1 2> {log}
    """


rule onedsfs_folded:
  input:
    'saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.saf.idx.done'
  output:
    sfs = 'saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.sfs'
  log:
    'log/saf_POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.sfs.log'
  threads: 12
  message:
    """ Optimize .saf.idx and calculate all pairwise 2dsfs's (joint site frequency spectra) using realSFS and fold for theta estimates"""
  shell:
    """
    module load singularity/2.x
    module load angsd/0.938
    singularity exec /apps/uibk/angsd/0.938/angsd.sandbox /opt/angsd/misc/realSFS saf/POPS/{wildcards.POPS}_{wildcards.IND}_{wildcards.MinDepth}_{wildcards.MaxDepth}.GL2.saf.idx -fold 1 > {output.sfs} 2> {log}
    """


rule sfs2theta:
  input:
    'saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.saf.idx.done',
    sfs = 'saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.sfs'
  output:
    'saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.sfs2theta.done'
  log:
    'log/saf_POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.sfs2theta.log'
  threads: 12
  message:
    """ calculate theta estimates for each site """
  shell:
    """
    module load singularity/2.x
    module load angsd/0.938
    singularity exec /apps/uibk/angsd/0.938/angsd.sandbox /opt/angsd/misc/realSFS saf2theta saf/POPS/{wildcards.POPS}_{wildcards.IND}_{wildcards.MinDepth}_{wildcards.MaxDepth}.GL2.saf.idx -sfs {input.sfs} -outname saf/POPS/{wildcards.POPS}_{wildcards.IND}_{wildcards.MinDepth}_{wildcards.MaxDepth}.GL2 > {output} 2> {log}
    """


rule theta_stat_chrom:
  input:
    'saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.sfs2theta.done'
  output:
    touch('saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.theta_stats_chrom.done')
  log:
    'log/saf_POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.theta_stats_chrom.log'
  threads: 12
  message:
    """ Estimate Tajimas D and other statistics """
  shell:
    """
    module load singularity/2.x
    module load angsd/0.938
    singularity exec /apps/uibk/angsd/0.938/angsd.sandbox /opt/angsd/misc/thetaStat do_stat saf/POPS/{wildcards.POPS}_{wildcards.IND}_{wildcards.MinDepth}_{wildcards.MaxDepth}.GL2.thetas.idx 2> {log}
    """


rule theta_stat_SW:
  input:
    'saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.sfs2theta.done'
  output:
    'saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.theta_stats_Window10kb.gz.pestPG.pestPG'
  log:
    'log/saf_POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.theta_stats_Window10kb.log'
  threads: 12
  message:
    """ Estimate Tajimas D and other statistics using a sliding window analysis """
  shell:
    """
    module load singularity/2.x
    module load angsd/0.938
    singularity exec /apps/uibk/angsd/0.938/angsd.sandbox /opt/angsd/misc/thetaStat do_stat saf/POPS/{wildcards.POPS}_{wildcards.IND}_{wildcards.MinDepth}_{wildcards.MaxDepth}.GL2.thetas.idx -win 10000 -step 1000 -outnames {output} 2> {log}
    """


rule twodsfs_folded:
  input:
    #pop1 = 'saf/POPS/df_{POP1}.GL2.saf.idx.done',
    #pop2 = 'saf/POPS/df_{POP2}.GL2.saf.idx.done'
    #pop1 = 'saf/POPS/{POP1}_{IND}_{MinDepth}_{MaxDepth}.GL2.saf.idx.done',
    #pop2 = 'saf/POPS/{POP2}_{IND}_{MinDepth}_{MaxDepth}.GL2.saf.idx.done'
  output:
    pop_pair = 'saf/POPS/{POP1}_vs_{POP2}_2D.sfs'
  log:
    pop_pair = 'log/saf_POPS/{POP1}_vs_{POP2}_2D.sfs.log'
  threads: 12
  message:
    """ Optimize .saf.idx and calculate all pairwise 2dsfs's (joint site frequency spectra) using realSFS and fold for FST """
  shell:
    """
    module load singularity/2.x
    module load angsd/0.938
    singularity exec /apps/uibk/angsd/0.938/angsd.sandbox /opt/angsd/misc/realSFS saf/POPS/{wildcards.POP1}*.GL2.saf.idx saf/POPS/{wildcards.POP2}*.GL2.saf.idx -fold 1 > {output.pop_pair} 2> {log.pop_pair}
    """


rule Fst_index:
  input:
    #pop1 = 'saf/POPS/{POP1}.GL2.saf.idx.done',
    #pop2 = 'saf/POPS/{POP2}.GL2.saf.idx.done',
    pop_pair = 'saf/POPS/{POP1}_vs_{POP2}_2D.sfs'
  output:
    pop_pair = touch('saf/POPS/{POP1}_vs_{POP2}.FstIndex.done')
  log:
    pop_pair = 'log/saf_POPS/{POP1}_vs_{POP2}.FstIndex.log'
  threads: 12
  message:
    """ Index fst for easy window analysis """
  shell:
    """
    module load singularity/2.x
    module load angsd/0.938
    singularity exec /apps/uibk/angsd/0.938/angsd.sandbox /opt/angsd/misc/realSFS fst index saf/POPS/{wildcards.POP1}*.GL2.saf.idx saf/POPS/{wildcards.POP1}*.GL2.saf.idx -sfs {input.pop_pair} -fstout saf/POPS/{wildcards.POP1}_vs_{wildcards.POP2} 2> {log.pop_pair}
    """


rule Fst_global:
  input:
    pop_pair = 'saf/POPS/{POP1}_vs_{POP2}.FstIndex.done'
  output:
    pop_pair = touch('saf/POPS/{POP1}_vs_{POP2}.Fst_Global.done')
  log:
    pop_pair = 'saf/POPS/{POP1}_vs_{POP2}.Fst_Global.log'
  threads: 12
  message:
    """ Get the global Fst estimate """
  shell:
    """
    module load singularity/2.x
    module load angsd/0.938
    singularity exec /apps/uibk/angsd/0.938/angsd.sandbox /opt/angsd/misc/realSFS fst stats saf/POPS/{wildcards.POP1}_vs_{wildcards.POP2}.fst.idx > saf/POPS/{wildcards.POP1}_vs_{wildcards.POP2}.fst.global 2> {log.pop_pair}
    """


rule Fst_windows:
  input:
    pop_pair = 'saf/POPS/{POP1}_vs_{POP2}.Fst_Global.done'
  output:
    pop_pair = touch('saf/POPS/{POP1}_vs_{POP2}.Fst_Windows.done')
  log:
    pop_pair = 'log/saf_POPS/{POP1}_vs_{POP2}.Fst_Windows.log'
  threads: 12
  message:
    """ Get the Fst estimates per window """
  shell:
    """
    module load singularity/2.x
    module load angsd/0.938
    singularity exec /apps/uibk/angsd/0.938/angsd.sandbox /opt/angsd/misc/realSFS fst stats2 saf/POPS/{wildcards.POP1}_vs_{wildcards.POP2}.fst.idx -win 10000 -step 1000 > saf/POPS/{wildcards.POP1}_vs_{wildcards.POP2}.fst.Window10kb.txt 2> {log.pop_pair}
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
#
#
#rule angsd_saf_samples:
#  input:
#    ref = config["ref_rapid"],
#    bam = 'realigned/{sample}.realigned.bam'
#  output:
#    touch('saf/saf_samples/{sample}.GL2.saf.idx.done')
#  log:
#    'log/saf_samples/{sample}.GL2.saf.idx.log'
#  threads: 12
#  message:
#    """ Compute site allele frequency likelihood (.saf.idx) for single samples using angsd (required for the global heterozygosity estimate for single samples)"""
#  shell:
#    """
#    module load angsd/0.938
#    /apps/uibk/bin/sysconfcpus -n 12 angsd -i {input.bam} -out saf/saf_samples/{wildcards.sample}.GL2 -anc {input.ref} -ref {input.ref} -doSaf 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minQ 20 -minMapQ 30 -GL 2 2> {log}
#    """
#
#
#rule real_SFS_samples:
#  input:
#    'saf/saf_samples/{sample}.GL2.saf.idx.done'
#  output:
#    'saf/saf_samples/{sample}.GL2.est.ml'
#  log:
#    'log/saf_samples/{sample}.GL2.est.ml.log'
#  threads: 12
#  message:
#    """ Optimize .saf.idx and estimate the site frequency spectrum (SFS) using realSFS and fold (required for the global heterozygosity estimate for single samples) """
#  shell:
#    """
#    module load singularity/2.x
#    module load angsd/0.938
#    singularity exec /apps/uibk/angsd/0.938/angsd.sandbox /opt/angsd/misc/realSFS saf/saf_samples/{wildcards.sample}.GL2.saf.idx -fold 1 > {output} 2> {log}
#    """
