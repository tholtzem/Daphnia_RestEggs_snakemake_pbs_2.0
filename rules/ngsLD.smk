
rule prepare_genoFile:
  input:
    'angsd/angsd_GL2_cutoffs_maf0018_nInd55.done'
  output:
    'ngsLD/angsd_LC_GL2_cutoff_maf0018_nInd55_sub_geno.beagle.gz'
  log: 'log/prepare_genoFile.log'
  threads: 12
  message:
    """ Prepare beagle formatted genotype likelihood file generated from ANGSD (-doGlf 2) by removing the header row and the first three columns (i.e. positions, major allele, minor allele) """
  shell:
    """
    zcat angsd/angsd_LC_GL2_cutoff_maf0018_nInd55.beagle.gz | cut -f 4- | awk 'NR != 1' | awk 'NR % 50 == 0' | gzip  > {output} 2> {log}
    """

rule prepare_posFile:
  input:
    'angsd/angsd_GL2_cutoffs_maf0018_nInd55.done'
  output:
    'ngsLD/angsd_LC_GL2_cutoff_maf0018_nInd55_sub_pos.gz'
  log: 'log/prepare_posFile.log'
  threads: 12
  message:
    """ Prepare beagle formatted genotype likelihood file generated from ANGSD (-doGlf 2) by removing the header row and the first three columns (i.e. positions, major allele, minor allele) """
  shell:
    """
    zcat angsd/angsd_LC_GL2_cutoff_maf0018_nInd55.mafs.gz | cut -f 1,2 |  awk 'NR != 1' | awk 'NR % 50 == 0' | sed 's/:/_/g'| gzip > {output}
    """

rule run_ngsLD:
  input:
   position = 'ngsLD/angsd_LC_GL2_cutoff_maf0018_nInd55_sub_pos.gz',
   geno = 'ngsLD/angsd_LC_GL2_cutoff_maf0018_nInd55_sub_geno.beagle.gz'
  output:
    'ngsLD/angsd_LC_GL2_cutoff_maf0018_nInd55_sub.ld.gz'
  log: 'log/run_ngsLD.log'
  threads: 24
  message:
    """ Estimate LD using ngsLD, which employs a maximum likelihood approach to account for genotype uncertainty in low-coverage whole genome sequencing data. """
  shell:
    """
    module load ngsld/1.1.1
    ngsld --geno {input.geno} --pos {input.position} --probs --n_ind 55 --n_sites 128612 --max_kb_dist 0 --max_snp_dist 0 --rnd_sample 0.008 --n_threads {threads} | gzip --best > {output} 2> {log}
    """


rule run_LDpruning:
  input:
   'ngsLD/angsd_LC_GL2_cutoff_maf0018_nInd55_sub.ld.gz'
  output:
    'ngsLD/angsd_LC_GL2_cutoff_maf0018_nInd55_sub.unlinked.id'
  log: 'log/LD_pruned.log'
  message:
    """ Prune your data, remove SNPs with high LD """
  shell:
    """
    module load singularity/2.x
    zcat {input} | singularity exec /apps/uibk/ngsld/1.1.1/ngsLD.sandbox /opt/ngsLD/scripts/prune_graph.pl --max_kb_dist 2000 --min_weight 0.5 --out {output} --print_excl ngsLD/angsd_LC_GL2_cutoff_maf0018_nInd55_sub.excluded_nodes.csv
    """

#rule LDpruned_SNPlist:
#  input:
#   'ngsLD/angsd_LC_GL2_cutoff_maf0018_nInd55_sub.unlinked.id'
#   #arg2 = 'angsd/angsd_LC_GL2_cutoff.nInd55.mafs.gz',
#   #arg3 = 'list/LDpruned_snps.list'
#  output:
#    touch('LDpruned_SNPlist.done')
#  log: 'log/LDpruned_SNPlist.log'
#  message:
#    """ Make a list of LDpruned SNPs in R """
#  shell:
#    """
#    Rscript scripts/getLDpruned_SNPlist.R
#    """
