rule fastqc_raw:
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

rule multiqc_raw:
  input:
    #"KRAKEN2_RESULTS/qc/fastqc/{sample}_R1_trmdfilt.keep_fastqc.html",
    #"KRAKEN2_RESULTS/qc/fastqc/{sample}_R2_trmdfilt.keep_fastqc.html"
    #r1 = lambda wildcards: fastqcHome(wildcards.sample)[0]
    #r2 = lambda wildcards: fastqcHome(wildcards.sample)[1]
  output:
    "raw/qc/multiqc/"
  log: "log/raw_multiqc_report.log"
  message: """ Multiqc of raw reads """
  shell:
    """
    multiqc raw/qc/fastqc/ -o {output}
    """


rule bbduk_trm:
  input:
    r1 = lambda wildcards: getFqHome(wildcards.sample)[0],
    r2 = lambda wildcards: getFqHome(wildcards.sample)[1],
    adapters = config["adapters"]
  output:
    r1 = "trm/{sample}_R1.trmd.fq.gz",
    r2 = "trm/{sample}_R2.trmd.fq.gz",
  log:
    'log/bbduk_trm_{sample}.log'
  threads: 2
  message: """ --- High sensitivity adapter trimming and polyG removal and minlength --- """
  shell:
    """
    /apps/uibk/bin/sysconfcpus -n 12 bbduk.sh -Xms110g -Xmx110g in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} ref={input.adapters} trimpolygright=10 qtrim=rl trimq=10 ordered=t ktrim=r k=21 mink=10 hdist=2 hdist2=1 minlength=40 stats=trm/stats/{wildcards.sample}.stats1 refstats=trm/stats/{wildcards.sample}.refstats1 bhist=trm/hist/{wildcards.sample}.bhist1 qhist=trm/hist/{wildcards.sample}.qhist1 lhist=trm/hist/{wildcards.sample}.lhist1 tpe tbo 2> {log}
    """


rule bbduk_trmfilt:
  input:
    #r1 = lambda wildcards: getFqHome(wildcards.sample)[0],
    #r2 = lambda wildcards: getFqHome(wildcards.sample)[1],
    r1 = "trm/{sample}_R1.trmd.fq.gz",
    r2 = "trm/{sample}_R2.trmd.fq.gz",
    phiX = config["phiX"]
  output:
    r1 = "trm/{sample}_R1.trmdfilt.fq.gz",
    r2 = "trm/{sample}_R2.trmdfilt.fq.gz"
  log:
    'log/bbduk_trmfilt_{sample}.log'
  threads: 2
  message: """ --- Additional quality and PhiX filtering and removal of reads with more than 1 N --- """
  shell:
    """
    /apps/uibk/bin/sysconfcpus -n 12 bbduk.sh -Xms110g -Xmx110g in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} ref={input.phiX} maq=10 k=31 hdist=1 maxns=1 ordered=t stats=trm/stats/{wildcards.sample}.stats2 bhist=trm/hist/{wildcards.sample}.bhist2 qhist=trm/hist/{wildcards.sample}.qhist2 lhist=trm/hist/{wildcards.sample}.lhist2 tpe tbo 2> {log} """


rule fastqc_trm:
  input:
    #r1 = lambda wildcards: getTrmFiltHome(wildcards.sample)[0],
    #r2 = lambda wildcards: getTrmFiltHome(wildcards.sample)[1]
    r1 = "trm/{sample}_R1.trmdfilt.fq.gz",
    r2 = "trm/{sample}_R2.trmdfilt.fq.gz"
  output:
    "trm/qc/fastqc/{sample}_R1_trmdfilt_fastqc.html",
    "trm/qc/fastqc/{sample}_R1_trmdfilt_fastqc.zip",
    "trm/qc/fastqc/{sample}_R2_trmdfilt_fastqc.html",
    "trm/qc/fastqc/{sample}_R2_trmdfilt_fastqc.zip"
  log: "log/{sample}_trmdfilt_qc.log"
  threads: 1
  message: """ Quality check with FastQC after trimming with bbduk """
  shell:
    """
    fastqc -o trm/qc/fastqc/ -f fastq {input.r1} & fastqc -o trm/qc/fastqc/ -f fastq {input.r2} 2> {log}
    """

rule multiqc_trm:
  output:
    "trm/qc/multiqc/"
  log: "log/trm_multiqc_report.log"
  message: """ Multiqc of trimmed reads """
  shell:
    """
    multiqc trm/qc/multiqc/ -o {output}
    """


rule kraken2:
  input:
    r1 = "trm/{sample}_R1.trmdfilt.fq.gz",
    r2 = "trm/{sample}_R2.trmdfilt.fq.gz",
    DB = config["KRAKEN2_DB"]
  output:
    #r1_class = "KRAKEN2_RESULTS/{sample}_R1.trmdfilt.classified.fq",
    #r2_class = "KRAKEN2_RESULTS/{sample}_R2.trmdfilt.classified.fq",
    #r1_unclass = "KRAKEN2_RESULTS/{sample}_R1.trmdfilt.unclassified.fq",
    #r2_unclass = "KRAKEN2_RESULTS/{sample}_R2.trmdfilt.unclassified.fq"
    #report = "KRAKEN2_RESULTS/{sample}.kraken_conf_0.9.report",
    #out = "KRAKEN2_RESULTS/{sample}.kraken_conf_0.9.out"
    touch("KRAKEN2_RESULTS/{sample}_kraken2_classified.done")
  log: "log/{sample}_kraken2_classified.log"
  shell:
    """
    module load ncbi-blast/2.7.1
    /apps/uibk/bin/sysconfcpus -n 12 kraken2 --db {input.DB} --threads 12 --paired --confidence 0.9 --classified-out KRAKEN2_RESULTS/{wildcards.sample}_R#.trmdfilt.classified.fq --unclassified-out KRAKEN2_RESULTS/{wildcards.sample}_R#.trmdfilt.unclassified.fq --report KRAKEN2_RESULTS/{wildcards.sample}.kraken_conf_0.9.report --output KRAKEN2_RESULTS/{wildcards.sample}.kraken_conf_0.9.out {input.r1} {input.r2} 2> {log}
    """


rule process_kraken2:
  input:
    touched = "KRAKEN2_RESULTS/{sample}_kraken2_classified.done",
    r1 = "KRAKEN2_RESULTS/{sample}_R_1.trmdfilt.classified.fq",
    r2 = "KRAKEN2_RESULTS/{sample}_R_2.trmdfilt.classified.fq"
  output:
    r1 = "KRAKEN2_RESULTS/{sample}_R1.trmdfilt.classified_retain.fq",
    r2 = "KRAKEN2_RESULTS/{sample}_R2.trmdfilt.classified_retain.fq"
  log: "log/{sample}.classified_retained.log"
  message: """ Process kraken2 report: from the classified reads grep all reads that are classified as opisthokonta or higher tax ids as they could still contain Daphnia reads (the settings below are for the big database containing protozoa, plants, fungi and humans. Clades below opisthokonta do not need to be considered for this database (no reads until human stuff) but this needs to be adapted depending on database """
  shell:
    """
    grep -A3 --no-group-separator -e 'kraken:taxid|1$' -e 'kraken:taxid|131567$' -e 'kraken:taxid|2759$' -e 'kraken:taxid|33154$' {input.r1} > {output.r1} &&
    grep -A3 --no-group-separator -e 'kraken:taxid|1$' -e 'kraken:taxid|131567$' -e 'kraken:taxid|2759$' -e 'kraken:taxid|33154$' {input.r2} > {output.r2}
    """


rule join_kraken_reads:
  input:
    touched = "KRAKEN2_RESULTS/{sample}_kraken2_classified.done",
    unclassified_r1 = "KRAKEN2_RESULTS/{sample}_R_1.trmdfilt.unclassified.fq",
    unclassified_r2 = "KRAKEN2_RESULTS/{sample}_R_2.trmdfilt.unclassified.fq",
    retained_r1 = "KRAKEN2_RESULTS/{sample}_R1.trmdfilt.classified_retain.fq",
    retained_r2 = "KRAKEN2_RESULTS/{sample}_R2.trmdfilt.classified_retain.fq"
  output:
    r1 = "KRAKEN2_RESULTS/{sample}_R1.trmdfilt.keep.fq.gz",
    r2 = "KRAKEN2_RESULTS/{sample}_R2.trmdfilt.keep.fq.gz"
  log: "log/{sample}.join_kraken_reads.log"
  message: """ Join the unclassified reads and the retained classified reads, remove intermediate files if desired """
  shell:
    """
    cat {input.unclassified_r1} {input.retained_r1} | gzip > {output.r1}  &&
    cat {input.unclassified_r2} {input.retained_r2} | gzip > {output.r2} 2> {log} 
    """


#rule gzip_kept:
#  input:
#    "KRAKEN2_RESULTS/{sample}.keep.fq"
#  output:
#    "KRAKEN2_RESULTS/{sample}.keep.fq.gz"
#  log: "log/{sample}.gzip_keptkrakenreads.log"
#  message: """ Gzip kept reads """
#  shell:
#    """
#    gzip {input} 2> {log}
#    """

rule fastqc_KRAKEN:
  input:
    r1 = "KRAKEN2_RESULTS/{sample}_R1.trmdfilt.keep.fq.gz",
    r2 = "KRAKEN2_RESULTS/{sample}_R2.trmdfilt.keep.fq.gz"
  output:
    "KRAKEN2_RESULTS/qc/fastqc/{sample}_R1_trmdfilt.keep_fastqc.html",
    "KRAKEN2_RESULTS/qc/fastqc/{sample}_R2_trmdfilt.keep_fastqc.html"
    #"KRAKEN2_RESULTS/qc/fastqc/{sample}_R1_trmdfilt.keep_fastqc.zip",
    #"KRAKEN2_RESULTS/qc/fastqc/{sample}_R1_trmdfilt.keep_fastqc.zip"
  log: "log/{sample}.qcKRAKEN.log"
  threads: 12
  message: """ Fastqc quality check of kept reads after cleaning with KRAKEN """
  shell:
    """
    fastqc -o KRAKEN2_RESULTS/qc/fastqc/ -f fastq {input}
    """

rule multiqc_KRAKEN:
  input:
    #"KRAKEN2_RESULTS/qc/fastqc/{sample}_R1_trmdfilt.keep_fastqc.html",
    #"KRAKEN2_RESULTS/qc/fastqc/{sample}_R2_trmdfilt.keep_fastqc.html"
  output:
    "KRAKEN2_RESULTS/qc/multiqc/"
  log: "log/KRAKEN_multiqc_report.log"
  message: """ Multiqc of kept reads after cleaning with KRAKEN """
  shell:
    """
    multiqc KRAKEN2_RESULTS/qc/fastqc/ -o {output}
    """

#rule remove_intermediate:
#rm ${OUTPUT}/${id}_R1.trmdfilt.unclassified.fastq &&
#rm ${OUTPUT}/${id}_R2.trmdfilt.unclassified.fastq &&


rule bb_refIndex:
  input:
    ref_rapid = config['ref_rapid'],
    #ref_canu = config['ref_canu'],
    #ref_flye = config['ref_flye'],
    #ref_long = config['ref_long']
  output:
    'ref/genome/1/summary.txt'
  shell:
    """
    /apps/uibk/bin/sysconfcpus -n 12 bbmap.sh ref={input.ref_rapid}
    """

rule bbmap:
  input:
    kraken_r1 = lambda wildcards: getKrakenHome(wildcards.sample)[0],
    kraken_r2 = lambda wildcards: getKrakenHome(wildcards.sample)[1],
    ref = config['ref_rapid']
  output:
    "bbmap/rapid/{sample}.bam"
  log: "log/{sample}_bam.log"
  threads: 24
  message: """ --- Mapping reads to reference genome, convert 2 bam, exclude unmapped reads, only keep reads with minq => 20 --- """
  shell:
    """
    /apps/uibk/bin/sysconfcpus -n 24 bbmap.sh -Xmx200g t={threads} ref={input.ref} in1={input.kraken_r1} in2={input.kraken_r2} out=stdout.sam minid=0.76 k=13 bw=0 ordered=t rgid={wildcards.sample} rglb=Novogene rgsm={wildcards.sample} rgpl=ILLUMINA overwrite=f unpigz=t | samtools view -F 4 -Shu -q 20 | samtools sort - -o {output} 2> {log}
    """

rule bamIndex:
  input:
    "bbmap/rapid/{sample}.bam"
  output:
    "bbmap/rapid/{sample}.bam.bai"
  threads: 2
  message: """--- Indexing with samtools ---"""
  shell:
    """
    samtools index {input} {output}
    """

rule remove_duplicates:
  input:
    bam = "bbmap/rapid/{sample}.bam",
    bai = "bbmap/rapid/{sample}.bam.bai"
  output:
    deDup = 'deDup/{sample}.dedup.bam',
    metrics = 'deDup/{sample}.dedup.metrics.txt'
  log: 'log/{sample}.dedup.bam.log'
  threads: 12
  message: """--- Removing duplicates of bam files with Picard ---"""
  shell:
    """
    /apps/uibk/bin/sysconfcpus -n 12 java -Xmx110g -jar /home/uibk/c7701178/.conda/envs/eggs/share/picard-2.25.0-1/picard.jar MarkDuplicates TMP_DIR=./ MAX_RECORDS_IN_RAM=15000000 REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 CREATE_INDEX=true INPUT={input.bam} OUTPUT={output.deDup} METRICS_FILE={output.metrics} 2> {log}
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

rule refIndex:
  input:
    ref = config['ref_rapid']
  output:
    "ref/CAJNAO01.fasta.fai"
  log: 'log/CAJNAO01_refIndex.log'
  shell:
    """
    samtools faidx {input.ref} 2> {log}
    """

rule ref_Dict:
  input:
    ref = config['ref_rapid']
  output:
    'ref/CAJNAO01.fasta.dict'
  log: 'log/CAJNAO01_refDict.log'
  shell:
    """
    /apps/uibk/bin/sysconfcpus -n 12 java -Xmx110g -jar /home/uibk/c7701178/.conda/envs/eggs/share/picard-2.25.0-1/picard.jar CreateSequenceDictionary R={input.ref} O={output} 2> {log}
    """


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
    ref = config['ref_rapid']
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
    /apps/uibk/bin/sysconfcpus -n 12 java -Xmx110g -jar $GATK -T RealignerTargetCreator -R {input.ref} -I {input.clip} -o {output.intervals} -drf BadMate 2> {log}
    """

rule realign_indel:
  input:
    clip = 'deDup/{sample}.overlapclipped.bam',
    ref = config['ref_rapid'],
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
    /apps/uibk/bin/sysconfcpus -n 12 java -Xmx110g -jar $GATK -T IndelRealigner -R {input.ref} -I {input.clip} -targetIntervals {input.intervals} -o {output.realigned} --consensusDeterminationModel USE_READS 2> {log}
    """
