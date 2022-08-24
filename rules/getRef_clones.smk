rule clip_overlap:
  input:
    deDup = 'deDup/{sample}.dedup.bam'
    #log_files = 'log/{sample}.dedup.Indexref.log'
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


rule refIndex:
  input:
    ref = config['ref_rapid']
  output:
    "ref/dgal_ra_pb-target-and-other_ill-confilter.blobfilter.rmmt_sspace-lr3_lrgc3_pg3_pilon_3.simple.fasta.fai"
  log: 'log/dgal_ra_refIndex.log'
  shell:
    """
    samtools faidx {input.ref} 2> {log}
    """


rule ref_Dict:
  input:
    ref = config['ref_rapid']
  output:
    'ref/dgal_ra_pb-target-and-other_ill-confilter.blobfilter.rmmt_sspace-lr3_lrgc3_pg3_pilon_3.simple.dict'
  log: 'log/dgal_ra_refDict.log'
  shell:
    """
    /apps/uibk/bin/sysconfcpus -n 12 java -Xmx110g -jar /home/uibk/c7701178/.conda/envs/eggs2.0/share/picard-2.25.0-1/picard.jar CreateSequenceDictionary R={input.ref} O={output} 2> {log}
    """

    
rule ls_ClipBam:
  input:
    #clip = 'deDup/{sample}.overlapclipped.bam'
  output:
    'list/overlapclippedBAM.list'
  log: 'log/lsClipBam.log'
  message: """--- Creating a sample list of clipped bam files for indel-realignement ---"""
  shell:
    """
    ls deDup/*overlapclipped.bam > {output} 2> {log}
    """


rule list_indels:
  input:
    clip = 'list/overlapclippedBAM.list',
    ref = config['ref_rapid'],
    idx_ref = "ref/dgal_ra_pb-target-and-other_ill-confilter.blobfilter.rmmt_sspace-lr3_lrgc3_pg3_pilon_3.simple.fasta.fai",
    dict_ref = 'ref/dgal_ra_pb-target-and-other_ill-confilter.blobfilter.rmmt_sspace-lr3_lrgc3_pg3_pilon_3.simple.dict'
  output:
    indels = 'list/indels.list'
  log: 'log/listIndels.log'
  threads: 36
  message:
    """ Create list of potential indels """
  shell:
    """
    GATK=(/home/uibk/c7701178/.conda/envs/eggs2.0/opt/gatk-3.8/GenomeAnalysisTK.jar)
    /apps/uibk/bin/sysconfcpus -n 36 java -Xmx330g -jar $GATK -T RealignerTargetCreator -R {input.ref} -I {input.clip} -o {output.indels} -drf BadMate 2> {log}
    """


rule realign_indel:
  input:
    clip = 'deDup/{sample}.overlapclipped.bam',
    ref = config['ref_rapid'],
    indels = 'list/indels.list'
  output:
    realigned = 'realigned/{sample}.realigned.bam'
  log: 'log/{sample}.realigned.bam.log'
  threads: 12
  message:
    """ Realign in-dels """
  shell:
    """
    GATK=(/home/uibk/c7701178/.conda/envs/eggs2.0/opt/gatk-3.8/GenomeAnalysisTK.jar)
    /apps/uibk/bin/sysconfcpus -n 12 java -Xmx110g -jar $GATK -T IndelRealigner -R {input.ref} -I {input.clip} -targetIntervals {input.indels} -o {output.realigned} --consensusDeterminationModel USE_READS 2> {log}
    """


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


rule ls_depth:
  output:
    'list/depth.list'
  log: 'log/depth.list.log'
  message: """ --- Creating a sample list of samtools-depth files for Rscript --- """
  shell:
    """
    ls depth/*.realigned.bam.depth.gz | cut -f2 -d '/' > {output} 2> {log}
    """


rule read_depth:
  input:
    'list/depth.list'
  output:
    #pdf = "some.pdf"
    touch('genome_stats.done')
  log: 'log/genome_stats.log'
  threads: 12
  message:
    """ --- Running Rscript to plot the genome-wide distribution of coverage --- """
  shell:
    """
    Rscript scripts/read_depthChrom_R.R 2> {log} 
    """

rule plot_summary:
  input:
    'genome_stats.done'
  output:
    touch('plot_summary.done')
  log: 'log/plot_summary.log'
  message:
    """ --- Running Rscript to plot depth summary per sample and output depth filters --- """
  shell:
    """
    Rscript scripts/plot_summary.R 2> {log}
    """

rule new_metadata_table:
  input:
    'plot_summary.done'
  output:
    'list/new_metadata_table.done'
  message: """ --- Create a metadata table for the samples passing the mean depth treshold ( > 3) --- """
  shell:
    """
    ./scripts/create_newmetadata_list.sh 2> {log} 
    """
