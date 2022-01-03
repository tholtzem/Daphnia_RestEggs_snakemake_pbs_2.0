rule ls_ClipBam:
  input:
    #bam = 'bbmap/all/{sample}.bam'
  output:
    #touch('bbmap/all/{sample}.ls.done'),
    'list/overlapclippedBAM.list'
  log: 'log/lsBam.log'
  message: """--- Creating a sample list of clipped bam files for indel-realignement ---"""
  shell:
    """
    ls deDup/*minq20.overlapclipped.bam > {output} 2> {log}
    """


rule list_indels:
  input:
    #clip = 'deDup/{sample}.overlapclipped.bam',
    clip = 'list/overlapclippedBAM.list',
    ref = config["ref"]
  output:
    #intervals = 'list/{sample}.intervals'
    intervals = 'list/all_samples.intervals'
  log: 'log/listIndels.log'
  threads: 12
  message:
    """ Create list of potential in-dels """
  shell:
    """
    GATK=(/home/uibk/c7701178/.conda/envs/eggs/opt/gatk-3.8/GenomeAnalysisTK.jar)
    /apps/uibk/bin/sysconfcpus -n 12 java -Xmx120g -jar $GATK -T RealignerTargetCreator -R {input.ref} -I {input.clip} -o {output.intervals} -drf BadMate 2> {log}
    """

rule realign_indel:
  input:
    clip = 'deDup/{sample}.overlapclipped.bam',
    ref = config["ref"],
    intervals = 'list/all_samples.intervals'
  output:
    realigned = 'realigned/{sample}.realigned.bam'
  log: 'log/{sample}.realigned.bam.log'
  threads: 12
  message:
    """ Realign in-dels """
  shell:
    """
    GATK=(/home/uibk/c7701178/.conda/envs/eggs/opt/gatk-3.8/GenomeAnalysisTK.jar)
    /apps/uibk/bin/sysconfcpus -n 12 java -Xmx120g -jar $GATK -T IndelRealigner -R {input.ref} -I {input.clip} -targetIntervals {input.intervals} -o {output.realigned} --consensusDeterminationModel USE_READS 2> {log}
    """
