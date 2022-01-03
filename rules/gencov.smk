
rule genome_coverage_bed:
  input:
    'realigned/{sample}.realigned.bam'
  output:
    'bedtools/{sample}.realigned.genomecov.bed'
  threads: 12
  message:
    """ Computes BED summaries using bedtools """
  shell:
    """
    bedtools genomecov -ibam {input} > {output}
    """


rule plot_gencov:
  input:
    bed= 'bedtools/{sample}.realigned.genomecov.bed'
  output:
    pdf = "bedtools/plots/{sample}.realigned.genomecov.pdf"
  threads: 4
  message:
    """ Running Rscript to plot the genome-wide distribution of coverage """
  shell:
    """
    Rscript scripts/plot_gene_covs.R {input.bed} {output.pdf}
    """
    

