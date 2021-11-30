rule remove_duplicates:
  input:
    bam = 'bbmap/rapid/minq20/{sample}.bam'
  output:
    deDup = 'deDup/{sample}.dedup.bam',
    metrics = 'deDup/{sample}.dedup.metrics.txt'
  log: 'log/{sample}.dedup.bam.log'
  threads: 12
  message: """--- Removing duplicates of merged bam files with Picard ---"""
  shell:
    """
    /apps/uibk/bin/sysconfcpus -n 12 java -Xmx110g -jar /home/uibk/c7701178/.conda/envs/eggs/share/picard-2.25.0-1/picard.jar MarkDuplicates TMP_DIR=./ MAX_RECORDS_IN_RAM=15000000 REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 CREATE_INDEX=true INPUT={input.bam} OUTPUT={output.deDup} METRICS_FILE={output.metrics} 2> {log}
    """


rule mark_duplicates:
  input:
    mrgd = 'bbmap/rapid/minq20/{sample}.bam'
  output:
    mrkDup = 'markDup/{sample}.mrkdup.bam',
    metrics = 'markDup/{sample}.mrkdup.metrics.txt'
  log: 'log/{sample}.mrkdup.bam.log'
  threads: 12
  message:
    """--- Marking duplicates of merged bam files with Picard ---"""
  shell:
    """
    /apps/uibk/bin/sysconfcpus -n 12 java -Xmx110g -jar /home/uibk/c7701178/.conda/envs/eggs/share/picard-2.25.0-1/picard.jar MarkDuplicates TMP_DIR=./ MAX_RECORDS_IN_RAM=15000000 REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 CREATE_INDEX=true INPUT={input.mrgd} OUTPUT={output.mrkDup} METRICS_FILE={output.metrics} 2> {log}
    """

rule clip_overlap:
  input:
    deDup = 'deDup/{sample}.dedup.bam'
  output:
    clip = 'deDup/{sample}.overlapclipped.bam' 
  log: 'log/{sample}.overlapclipped.bam.log'
  threads: 12
  message:
    """ Clip overlapping paired end reads """
  shell:
    """
    bam clipOverlap --in {input.deDup} --out {output.clip} --stats 2> {log}
    """

rule Index_clippedBAM:
  input:
    clip = 'deDup/{sample}.overlapclipped.bam'
  output:
    idx = 'deDup/{sample}.overlapclipped.bam.bai'
  log: 'log/{sample}.overlapclipped.bam.log'
  threads: 2
  message: """--- Indexing clipped BAM files with samtools ---"""
  shell:
    """
    samtools index {input.clip} {output.idx} 2> {log}
    """







