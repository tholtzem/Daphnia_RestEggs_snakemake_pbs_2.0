
#rule newline2csv:
#  input:
#    pops = 'list/pop_list/{pops}.txt'
#  output:
#    pops = 'list/pop_list/{pops}.csv'
#  log: 'log/{pops}_newline2csv.log'
#  threads: 2
#  message: """ Convert newline character separated list to comma separated list """
#  shell:
#    """
#    cat {input.pops} | cut -f1 -d'.' | cut -f2 -d'/' | sed -z 's/\n/,/g;s/,$/\n/' > {output.pops} 2> {log}
#    """

rule get_fixedSites:
  input:
    vcf = 'ancestry/LC/data/{prefix}_hf.DP{genoDP}_vmiss20.vcf',
    pop1 = 'list/pop_list/longispina.csv',
    pop2 = 'list/pop_list/galeata.csv',
    hybrids = 'list/pop_list/hybrids_LC.csv'
  output:
    fixed_sites = 'ancestry/LC/{prefix}_hf.DP{genoDP}_fixed_sites.txt',
    report = 'ancestry/LC/{prefix}_hf.DP{genoDP}_get_fixedSites_sites.report'
  log:
    'log/get_fixedSites_{prefix}_hf.DP{genoDP}_fixed_sites.log'
  message: """ --- Run ruby script (https://github.com/mmatschiner/tutorials/blob/master/analysis_of_introgression_with_snp_data/src/get_fixed_site_gts.rhttps://github.com/mmatschiner/tutorials/blob/master/analysis_of_introgression_with_snp_data/src/get_fixed_site_gts.rbb) to get the alleles at sites that are fixed differently in two parental species (complete fixation of at least 100% in parental species) --- """
  threads: 12
  shell:
    """
    pop1=$(cat {input.pop1})
    pop2=$(cat {input.pop2})
    hy=$(cat {input.hybrids})
    ruby scripts/get_fixed_site_gts.rb {input.vcf} {output.fixed_sites} $pop1 $hy $pop2 1.0 > {output.report} 2> {log}
    """ 


rule plot_fixedSites:
  input:
    fixed_sites = 'ancestry/LC/{prefix}_hf.DP{genoDP}_fixed_sites.txt',
  output:
    svg = 'ancestry/LC/{prefix}_hf.DP{genoDP}_fixed_sites_80percent_thinned1000bp.svg',
    report = 'ancestry/LC/{prefix}_hf.DP{genoDP}_plot_fixed_sites_80percent_thinned1000bp.report'
  log:
    'log/plot_fixedSites_{prefix}_hf.DP{genoDP}_fixed_sites_80percent.log'
  message: """ --- Run ruby script (https://github.com/mmatschiner/tutorials/blob/master/analysis_of_introgression_with_snp_data/src/get_fixed_site_gts.rhttps://github.com/mmatschiner/tutorials/blob/master/analysis_of_introgression_with_snp_data/src/plot_fixed_site_gts.rbb) to get the alleles at sites that are fixed differently in two parental species (complete fixation of 80-100% in all individuals) --- """
  threads: 12
  shell:
    """
    ruby scripts/plot_fixed_site_gts.rb {input.fixed_sites} {output.svg} 0.8 1000 | head -n -3 | tail -n+2 | awk -v OFS='\t' '{{$1=$1; print}}' > {output.report} 2> {log}
    """


rule list_fixedSites:
  input:
    'ancestry/fixed_sites.txt'
  output:
    'ancestry/fixed_sites.list'
  log: 'list_fixedSites.log'
  threads: 12
  message: """ Reformat list with fixed sites """
  shell:
    """
    cat {input} | grep -e 'HiC_scaffold' | cut -f1,2 > {output}
    """


rule index_fixedSites:
  input:
    'ancestry/fixed_sites.list'
  output:
    touch('ancestry/index_fixed_sites.done')
  log: 'log/index_fixed_sites.log'
  threads: 12
  message:
    """ Index list of fixed sites """
  shell:
    """
    module load angsd/0.938
    /apps/uibk/bin/sysconfcpus -n {threads} angsd sites index {input} 2> {log}
    """



rule saf_fixedSites:
  input:
    ref = config["ref_rapid"],
    bam = 'realigned/{sample}.realigned.bam',
    fixed_sites = 'ancestry/fixed_sites.list',
    touched = 'ancestry/index_fixed_sites.done'
  output:
    touch('ancestry/{sample}.GL2.saf.idx.done')
  log:
    'log/saf_fixedSites_{sample}.GL2.saf.idx.log'
  threads: 12
  message:
    """ Compute site allele frequency likelihood (.saf.idx) for single samples using angsd (required for the estimation of heterozygosity from single samples and fixed sites)"""
  shell:
    """
    module load angsd/0.938
    /apps/uibk/bin/sysconfcpus -n 12 angsd -i {input.bam} -ref {input.ref} -anc {input.ref} -out ancestry/{wildcards.sample}.GL2 -doSaf 1 -GL 2 -sites {input.fixed_sites} 2> {log}
    """


rule realSFS_samples_fixedSites:
  input:
    touched1 = 'ancestry/{sample}.GL2.saf.idx.done',
    touched2 = 'ancestry/index_fixed_sites.done'
  output:
    'ancestry/{sample}.GL2.est.ml'
  log:
    'log/saf_fixedSites_{sample}.GL2.est.ml.log'
  threads: 12
  message:
    """ Optimize .saf.idx and estimate the site frequency spectrum (SFS) using realSFS and fold (required for the global heterozygosity estimate for single samples) """
  shell:
    """
    module load singularity/2.x
    module load angsd/0.938
    singularity exec /apps/uibk/angsd/0.938/angsd.sandbox /opt/angsd/misc/realSFS ancestry/{wildcards.sample}.GL2.saf.idx -fold 1 > {output} 2> {log}
    """
