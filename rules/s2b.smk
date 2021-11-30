
#rule sam2bam:
#  input:
#    sam = '/home/uibk/c7701125/scratch/EXCHANGE/RESTEGGS_LC/MAPPING/bbmap/{sample}.sam.gz'
#  output:
#    bam = 'bbmap/rapid/{sample}.bam'
#  message: """--- Converting sam files to sorted bam ---"""
#  threads: 12 
#  shell:
#    """
#    samtools view -F 4 -Shu {input.sam} | samtools sort - -o {output.bam}
#    """

#rule getRefClones:
#  input:
#    #'ref_clones.txt'
#  output:
#    touch('getRefClones.done')
#  message: """ Get a set of bam files from reference clone data set (extant samples) """
#  shell:
#    """
#    rsync -avP --files-from=ref_clones.txt /home/uibk/c7701178/scratch/DAPHNIA/daphnia_snakemake_pbs/bbmap/rapid/ bbmap/rapid/
#    """

rule bam_q20:
  input:
    bam = 'bbmap/rapid/{sample}.bam'
  output:
    bamq20 = 'bbmap/rapid/minq20/{sample}.minq20.bam'
  message: """--- Filter bam reads with a mapping quality lower than 20 and convert to sorted bam ---"""
  threads: 12 
  shell:
    """
    samtools view -Shu -q 20 {input.bam} | samtools sort - -o {output.bamq20}
    """

