rule fastqc:
  input:
    r1 = lambda wildcards: getFqHome(wildcards.sample)[0],
    r2 = lambda wildcards: getFqHome(wildcards.sample)[1]
  output:
    "raw/qc/fastqc/{sample}_R1_001_fastqc.html",
    "raw/qc/fastqc/{sample}_R1_001_fastqc.zip",
    "raw/qc/fastqc/{sample}_R2_001_fastqc.html",
    "raw/qc/fastqc/{sample}_R2_001_fastqc.zip"
  log: "log/{sample}.qc.log"
  threads: 1
  message: """ Quality check of raw data with FastQC before trimming. """
  shell:
    """
    fastqc -o raw/qc/fastqc/ -f fastq {input.r1} & fastqc -o raw/qc/fastqc/ -f fastq {input.r2}
    """

