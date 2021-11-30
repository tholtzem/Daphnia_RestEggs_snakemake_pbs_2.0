rule angsd_saf:
  input:
    ref = config["ref"],
    realigned = 'realigned/{sample}.minq20.realigned.bam'
  output:
    touch('angsd/{sample}.GL2.saf.idx.done')
  log:
    'log/{sample}.GL2.saf.log'
  threads: 24
  message:
    """ SFS Estimation for single samples - part1 """
  shell:
    """
    angsd -i {input.realigned} -out angsd/{wildcards.sample}.GL2 -anc {input.ref} -dosaf 1 -fold 1 -ref {input.ref} -C 50 -minQ 20 -minMapQ 30 -GL 2 2> {log}
    """

rule real_SFS:
  input:
    'angsd/{sample}.GL2.saf.idx.done'
  output:
    'angsd/{sample}.GL2.est.ml'
  log:
    'log/{sample}.GL2.est.ml.log'
  threads: 24
  message:
    """ SFS Estimation for single samples - part1 """
  shell:
    """
    path=(/home/uibk/c7701178/.conda/envs/eggs/bin/)
    $path/realSFS angsd/{wildcards.sample}.GL2.saf.idx > {output} 2> {log}
    """

rule plotHeterozygosity:
  input: 'angsd/{sample}.GL2.est.ml'
  output:
    'angsd/{sample}.GL2.heterozygosity.txt'
  log:
    'log/{sample}.GL2.heterozygosity.log'
  message:
    """ plot heterozygosity """
  shell:
    """
   Rscript scripts/plotHeterozygosity.R {input} {output} 2> {log}
    """




