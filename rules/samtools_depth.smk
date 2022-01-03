
rule samtools_depth:
  input:
    realigned = 'realigned/{sample}.realigned.bam'
  output:
    depth = 'depth/{sample}.realigned.bam.depth.gz'
  log: 'log/{sample}.realigned.bam.depth.log'
  threads: 12
  message:
    """ Count per position depth per sample using samtools depth """
  shell:
    """
    samtools depth -aa {input.realigned} | cut -f3 | gzip > {output.depth} 2> {log}
    """
