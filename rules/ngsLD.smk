
rule prepare_genoFile:
  input:
    'angsd/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done'
  output:
    'ngsLD/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_sub_geno.beagle.gz'
  log: 'log/{sets}/prepare_subgenoFile_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.log'
  threads: 12
  message:
    """ Prepare beagle formatted genotype likelihood file generated from ANGSD (-doGlf 2) by removing the header row and the first three columns (i.e. positions, major allele, minor allele) """
  shell:
    """
    zcat angsd/{wildcards.sets}/angsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth}.beagle.gz | cut -f 4- | awk 'NR != 1' | awk 'NR % 50 == 0' | gzip  > {output} 2> {log}
    """


rule prepare_posFile:
  input:
    'angsd/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done'
  output:
    'ngsLD/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_sub_pos.gz'
  log: 'log/{sets}/prepare_subposFile_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.log'
  threads: 12
  message:
    """ Prepare position file generated from ANGSD (-doGlf 2) by removing the header row and the first three columns (i.e. positions, major allele, minor allele) """
  shell:
    """
    zcat angsd/{wildcards.sets}/angsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth}.mafs.gz | cut -f 1,2 |  awk 'NR != 1' | awk 'NR % 50 == 0' | sed 's/:/_/g'| gzip > {output}
    """


rule run_ngsLD:
  input:
   position = 'ngsLD/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_sub_pos.gz',
   geno = 'ngsLD/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_sub_geno.beagle.gz'
  output:
    'ngsLD/{sets}/run_ngsLD_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_sub.ld.gz'
  log: 'log/{sets}/run_ngsLD_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.log'
  threads: 36
  message:
    """ Estimate LD using ngsLD, which employs a maximum likelihood approach to account for genotype uncertainty in low-coverage whole genome sequencing data. """
  shell:
    """
    module load ngsld/1.1.1
    N_SITES=`zcat {input.position} | wc -l` &&
    ngsld --geno {input.geno} --pos {input.position} --probs --n_ind {wildcards.IND} --n_sites $N_SITES --max_kb_dist 0 --max_snp_dist 0 --rnd_sample 0.01 --n_threads {threads} | gzip --best > {output} 2> {log}
    """


rule run_LDpruning:
  input:
    'ngsLD/{sets}/run_ngsLD_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_sub.ld.gz'
  output:
   touch('ngsLD/{sets}/run_LDpruning_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done')
  log: 'log/{sets}/run_LDpruning_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.log'
  threads: 12
  message:
    """ Prune your data, remove SNPs with high LD """
  shell:
    """
    module load singularity/2.x
    zcat {input} | singularity exec /apps/uibk/ngsld/1.1.1/ngsLD.sandbox /opt/ngsLD/scripts/prune_graph.pl --max_kb_dist 2000 --min_weight 0.5 --out ngsLD/{wildcards.sets}/angsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth}_sub_unlinked.id --print_excl ngsLD/{wildcards.sets}/angsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth}_sub_excluded_nodes.csv
    """


rule LDpruned_SNPlist:
  input:
    touched = 'ngsLD/{sets}/run_LDpruning_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done'
  output:
    arg3 = 'ngsLD/{sets}/LDpruned_snps_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.list'  
  log: 'log/{sets}/LDpruned_SNPlist_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.log'
  message:
    """ Make a list of LDpruned SNPs in R and index SNP list """
  shell:
    """
    Rscript scripts/getLDpruned_SNPlist.R ngsLD/{wildcards.sets}/angsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth}_sub_unlinked.id angsd/{wildcards.sets}/angsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth}.mafs.gz {output.arg3} 2> {log}
    """


#rule LDblocks:
#  input:
#   'ngsLD/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_sub.ld.gz'
#  output:
#   touch('ngsLD/{sets}/run_LDpruning_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done')
#  log: 'log/{sets}/LD_pruned_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.log'
#  threads: 12
#  message:
#    """ Prune your data, remove SNPs with high LD """
#  shell:
#    """
#    module load singularity/2.x
#    zcat {input} | singularity exec /apps/uibk/ngsld/1.1.1/ngsLD.sandbox /opt/ngsLD/scripts/prune_graph.pl --max_kb_dist 2000 --min_weight 0.5 --out ngsLD/{wildcards.sets}/angsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth}_sub_unlinked.id --print_excl ngsLD/{wildcards.sets}/angsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth}_sub_excluded_nodes.csv
#    """

