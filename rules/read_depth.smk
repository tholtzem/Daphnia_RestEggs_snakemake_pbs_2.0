#rule ls_realignedBam:
#  input:
#    #bam = 'bbmap/all/{sample}.bam'
#  output:
#    #touch('bbmap/all/{sample}.ls.done'),
#    'list/realignedBAM.list'
#  log: 'log/realignedBam.list.log'
#  message: """--- Creating a sample list of realigned bam files ---"""
#  shell:
#    """
#    ls realigned/*minq20.realigned.bam > {output} 2> {log}
#    """

#rule ls_depth:
#  input:
#    #'depth/{sample}.depth.gz'
#  output:
#    #touch('bbmap/all/{sample}.ls.done'),
#    'list/depth.list'
#  log: 'log/depth.list.log'
#  message: """--- Creating a sample list of samtools-depth files for Rscript ---"""
#  shell:
#    """
#    ls depth/*depth.gz | cut -f2 -d '/' > {output} 2> {log}
#    """


rule read_depth:
  #input:
    #test.list
  output:
    #pdf = "some.pdf"
    touch('genome_stats.done')
  log: 'log/genome_stats.done.log'
  threads: 12
  message:
    """ Running Rscript to plot the genome-wide distribution of coverage """
  shell:
    """
    Rscript scripts/read_depthChrom_R.R 2> {log}
    """

