localrules: ls_ClipBam, ls_ClipBam_REFs, sym_link_REFs
#rule fastqc_raw:
#  input:
#    r1 = lambda wildcards: getFqHome(wildcards.sample)[0],
#    r2 = lambda wildcards: getFqHome(wildcards.sample)[1]
#  output:
#    "raw/qc/fastqc/{sample}_R1_001_fastqc.html",
#    "raw/qc/fastqc/{sample}_R1_001_fastqc.zip",
#    "raw/qc/fastqc/{sample}_R2_001_fastqc.html",
#    "raw/qc/fastqc/{sample}_R2_001_fastqc.zip"
#  log: "log/{sample}.qc.log"
#  threads: 12
#  message: """ Quality check of raw data with FastQC before trimming. """
#  shell:
#    """
#    fastqc -o raw/qc/fastqc/ -f fastq {input.r1} & fastqc -o raw/qc/fastqc/ -f fastq {input.r2}
#    """
#
#
#rule multiqc_raw:
#  input:
#  output:
#    "raw/qc/multiqc/"
#  log: "log/raw_multiqc_report.log"
#  message: """ Multiqc of raw reads """
#  shell:
#    """
#    multiqc raw/qc/fastqc/ -o {output}
#    """
#
#
#rule bbduk_trm:
#  input:
#    r1 = lambda wildcards: getFqHome(wildcards.sample)[0],
#    r2 = lambda wildcards: getFqHome(wildcards.sample)[1],
#    adapters = config['adapters']
#  output:
#    r1 = "trm/{sample}_R1.trmd.fq.gz",
#    r2 = "trm/{sample}_R2.trmd.fq.gz"
#  log:
#    'log/bbduk_trm_{sample}.log'
#  threads: 12
#  message: """ --- High sensitivity adapter trimming and polyG removal and minlength --- """
#  shell:
#    """
#    /apps/uibk/bin/sysconfcpus -n 12 bbduk.sh -Xms110g -Xmx110g in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} ref={input.adapters} trimpolygright=10 qtrim=rl trimq=10 ordered=t ktrim=r k=21 mink=10 hdist=2 hdist2=1 minlength=40 stats=trm/stats/{wildcards.sample}.stats1 refstats=trm/stats/{wildcards.sample}.refstats1 bhist=trm/hist/{wildcards.sample}.bhist1 qhist=trm/hist/{wildcards.sample}.qhist1 lhist=trm/hist/{wildcards.sample}.lhist1 tpe tbo 2> {log}
#    """
#
#
#rule bbduk_trmfilt:
#  input:
#    r1 = "trm/{sample}_R1.trmd.fq.gz",
#    r2 = "trm/{sample}_R2.trmd.fq.gz",
#    phiX = config["phiX"]
#  output:
#    r1 = "trm/{sample}_R1.trmdfilt.fq.gz",
#    r2 = "trm/{sample}_R2.trmdfilt.fq.gz"
#  log:
#    'log/bbduk_trmfilt_{sample}.log'
#  threads: 12
#  message: """ --- Additional quality and PhiX filtering and removal of reads with more than 1 N --- """
#  shell:
#    """
#    /apps/uibk/bin/sysconfcpus -n 12 bbduk.sh -Xms110g -Xmx110g in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} ref={input.phiX} maq=10 k=31 hdist=1 maxns=1 ordered=t stats=trm/stats/{wildcards.sample}.stats2 bhist=trm/hist/{wildcards.sample}.bhist2 qhist=trm/hist/{wildcards.sample}.qhist2 lhist=trm/hist/{wildcards.sample}.lhist2 tpe tbo 2> {log} """
#
#
#rule fastqc_trm:
#  input:
#    r1 = "trm/{sample}_R1.trmdfilt.fq.gz",
#    r2 = "trm/{sample}_R2.trmdfilt.fq.gz"
#  output:
#    "trm/qc/fastqc/{sample}_R1.trmdfilt_fastqc.html",
#    "trm/qc/fastqc/{sample}_R1.trmdfilt_fastqc.zip",
#    "trm/qc/fastqc/{sample}_R2.trmdfilt_fastqc.html",
#    "trm/qc/fastqc/{sample}_R2.trmdfilt_fastqc.zip"
#  log: "log/{sample}_trmdfilt_qc.log"
#  threads: 1
#  message: """ Quality check with FastQC after trimming with bbduk """
#  shell:
#    """
#    fastqc -o trm/qc/fastqc/ -f fastq {input.r1} & fastqc -o trm/qc/fastqc/ -f fastq {input.r2} 2> {log}
#    """
#
#
#rule multiqc_trm:
#  output:
#    "trm/qc/multiqc/"
#  log: "log/trm_multiqc_report.log"
#  message: """ Multiqc of trimmed reads """
#  shell:
#    """
#    multiqc trm/qc/fastqc/ -o {output}
#    """
#
#
#rule kraken2:
#  input:
#    r1 = "trm/{sample}_R1.trmdfilt.fq.gz",
#    r2 = "trm/{sample}_R2.trmdfilt.fq.gz",
#    DB = config["KRAKEN2_DB"]
#  output:
#     touch("KRAKEN2_RESULTS/{sample}_kraken2_classified.done")
#  log: "log/{sample}_kraken2_classified.log"
#  shell:
#    """
#    module load ncbi-blast/2.7.1
#    /apps/uibk/bin/sysconfcpus -n 12 kraken2 --db {input.DB} --threads 12 --paired --confidence 0.9 --classified-out KRAKEN2_RESULTS/{wildcards.sample}_R#.trmdfilt.classified.fq --unclassified-out KRAKEN2_RESULTS/{wildcards.sample}_R#.trmdfilt.unclassified.fq --report KRAKEN2_RESULTS/{wildcards.sample}.kraken_conf_0.9.report --output KRAKEN2_RESULTS/{wildcards.sample}.kraken_conf_0.9.out {input.r1} {input.r2} 2> {log}
#    """
#
#
#rule process_kraken2:
#  input:
#    touched = "KRAKEN2_RESULTS/{sample}_kraken2_classified.done",
#    r1 = "KRAKEN2_RESULTS/{sample}_R_1.trmdfilt.classified.fq",
#    r2 = "KRAKEN2_RESULTS/{sample}_R_2.trmdfilt.classified.fq"
#  output:
#    r1 = "KRAKEN2_RESULTS/{sample}_R1.trmdfilt.classified_retain.fq",
#    r2 = "KRAKEN2_RESULTS/{sample}_R2.trmdfilt.classified_retain.fq"
#  log: "log/{sample}.classified_retained.log"
#  message:
#    """
#    Process kraken2 report: from the classified reads grep all reads that are classified as opisthokonta or higher tax ids as they could still contain Daphnia reads (the settings below are for the big database containing protozoa, plants, fungi and humans. Clades below opisthokonta do not need to be considered for this database (no reads until human stuff) but this needs to be adapted depending on database
#    """
#  shell:
#    """
#    grep -A3 --no-group-separator -e 'kraken:taxid|1$' -e 'kraken:taxid|131567$' -e 'kraken:taxid|2759$' -e 'kraken:taxid|33154$' {input.r1} > {output.r1} && grep -A3 --no-group-separator -e 'kraken:taxid|1$' -e 'kraken:taxid|131567$' -e 'kraken:taxid|2759$' -e 'kraken:taxid|33154$' {input.r2} > {output.r2} 2> {log}
#    """
#
#
#rule join_kraken_reads:
#  input:
#    touched = "KRAKEN2_RESULTS/{sample}_kraken2_classified.done",
#    unclassified_r1 = "KRAKEN2_RESULTS/{sample}_R_1.trmdfilt.unclassified.fq",
#    unclassified_r2 = "KRAKEN2_RESULTS/{sample}_R_2.trmdfilt.unclassified.fq",
#    retained_r1 = "KRAKEN2_RESULTS/{sample}_R1.trmdfilt.classified_retain.fq",
#    retained_r2 = "KRAKEN2_RESULTS/{sample}_R2.trmdfilt.classified_retain.fq"
#  output:
#    r1 = "KRAKEN2_RESULTS/{sample}_R1.trmdfilt.keep.fq.gz",
#    r2 = "KRAKEN2_RESULTS/{sample}_R2.trmdfilt.keep.fq.gz"
#  log: "log/{sample}.join_kraken_reads.log"
#  message: """ Join the unclassified reads and the retained classified reads, remove intermediate files if desired """
#  shell:
#    """
#    cat {input.unclassified_r1} {input.retained_r1} | gzip > {output.r1}  &&
#    cat {input.unclassified_r2} {input.retained_r2} | gzip > {output.r2} 2> {log} 
#    """
#
#
#rule fastqc_KRAKEN:
#  input:
#    r1 = "KRAKEN2_RESULTS/{sample}_R1.trmdfilt.keep.fq.gz",
#    r2 = "KRAKEN2_RESULTS/{sample}_R2.trmdfilt.keep.fq.gz"
#  output:
#    "KRAKEN2_RESULTS/qc/fastqc/{sample}_R1.trmdfilt.keep_fastqc.html",
#    "KRAKEN2_RESULTS/qc/fastqc/{sample}_R2.trmdfilt.keep_fastqc.html",
#    "KRAKEN2_RESULTS/qc/fastqc/{sample}_R1.trmdfilt.keep_fastqc.zip",
#    "KRAKEN2_RESULTS/qc/fastqc/{sample}_R2.trmdfilt.keep_fastqc.zip"
#  log: "log/{sample}.qcKRAKEN.log"
#  threads: 12
#  message: """ Fastqc quality check of kept reads after cleaning with KRAKEN """
#  shell:
#    """
#    fastqc -o KRAKEN2_RESULTS/qc/fastqc/ -f fastq {input.r1} & fastqc -o KRAKEN2_RESULTS/qc/fastqc/ -f fastq {input.r2} 2> {log}
#    """
#
#
#rule multiqc_KRAKEN:
#  input:
#    #"KRAKEN2_RESULTS/qc/fastqc/{sample}_R1_trmdfilt.keep_fastqc.html",
#    #"KRAKEN2_RESULTS/qc/fastqc/{sample}_R2_trmdfilt.keep_fastqc.html"
#  output:
#    "KRAKEN2_RESULTS/qc/multiqc/"
#  log: "log/KRAKEN_multiqc_report.log"
#  message: """ Multiqc of kept reads after cleaning with KRAKEN """
#  shell:
#    """
#    multiqc KRAKEN2_RESULTS/qc/fastqc/ -o {output}
#    """
#

rule bb_refIndex:
  input:
    #ref_rapid = config['ref_rapid'],
    ref = config['ref_HiC']
  output:
    touch('ref/bb_indexRef.done')
  log: 'log/bb_indexRef.log'
  shell:
    """
    /apps/uibk/bin/sysconfcpus -n 12 bbmap.sh -Xmx120g ref={input.ref}
    """


#rule bbmap:
#  input:
#    kraken_r1 = lambda wildcards: getKrakenHome(wildcards.sample)[0],
#    kraken_r2 = lambda wildcards: getKrakenHome(wildcards.sample)[1],
#    ref = config['ref_HiC'],
#    idx = 'ref/bb_indexRef.done',
#    SAMPLETABLE = 'list/da_eggs_hatchlings_pelagial_metadata.csv'
#  output:
#    bam = 'bbmap/{sample}.bam',
 #   bhist = 'bbmap/stats/bhist/{sample}.bhist.txt',
 #   qhist = 'bbmap/qhist/{sample}.qhist.txt',
 ##   lhist = 'bbmap/lhist/{sample}.lhist.txt',
 #   covstats = 'bbmap/cov/{sample}.covstats.txt',
 #   covhist = 'bbmap/cov/{sample}.covhist.txt',
 #   basecov = 'bbmap/cov/{sample}.basecov.txt',
 #   bincov = 'bbmap/cov/{sample}.bincov.txt'
 # log: 'log/bbmap_HiC_{sample}_bam.log' 
 # threads: 24
 # message: """ --- Mapping reads to reference genome, convert 2 bam, exclude unmapped reads, only keep reads with minq => 20 --- """
 # shell:
 #   """
 #   id=`echo {input.kraken_r1} | sed -e "s/_R1.trmdfilt.keep.fq.gz$//" | cut -f7 -d '/'`
 #   echo $id
 #   
 #   RG_ID=`grep -P "${{id}}\t" {input.SAMPLETABLE} | cut -f4`
 #   RG_LB=LIB1of_${{RG_ID}}
 #   RG_SM=$RG_ID
 #   RG_PL=ILLUMINA
 #   RG_PU=LIB1of_${{RG_ID}}
 #                                   
 #   echo $RG_ID
 #   echo $RG_LB
 #   echo $RG_SM
 #   echo $RG_PL
 #   echo $RG_PU
 #   
 #   /apps/uibk/bin/sysconfcpus -n 24 bbmap.sh -Xmx200g t={threads} ref={input.ref} in1={input.kraken_r1} in2={input.kraken_r2} out=stdout.sam minid=0.95 k=13 bw=0 ordered=t rgid=$RG_ID rglb=$RG_LB rgsm=$RG_SM rgpl=$RG_PL rgpu=$RG_PU overwrite=f unpigz=t bhist={output.bhist} qhist={output.qhist} lhist={output.lhist} covstats={output.covstats} covhist={output.covhist} basecov={output.basecov} bincov={output.bincov} | samtools view -F 4 -Shu -q 20 | samtools sort - -o {output.bam} 2> {log}
 #   """
#
#
#rule bamIndex:
#  input:
#    'bbmap/{sample}.bam'
#  output:
#    'bbmap/{sample}.bam.bai'
#  threads: 12
#  log: "log/indexBam_HiC_{sample}.log"
#  message: """--- Indexing with samtools ---"""
#  shell:
#    """
#    samtools index {input} {output} 2> {log}
#    """
#
#
#rule remove_duplicates:
#  input:
#    bam = "bbmap/{sample}.bam",
#    bai = "bbmap/{sample}.bam.bai"
#  output:
#    deDup = 'deDup/{sample}.dedup.bam',
#    metrics = 'deDup/{sample}.dedup.metrics.txt'
#  log: 'log/{sample}.dedup.bam.log'
#  threads: 12
#  message: """--- Removing duplicates of bam files with Picard ---"""
#  shell:
#    """
#    /apps/uibk/bin/sysconfcpus -n 12 java -Xmx110g -jar /home/uibk/c7701178/.conda/envs/eggs2.0/share/picard-2.25.0-1/picard.jar MarkDuplicates TMP_DIR=./ MAX_RECORDS_IN_RAM=15000000 REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 CREATE_INDEX=true INPUT={input.bam} OUTPUT={output.deDup} METRICS_FILE={output.metrics} 2> {log}
#    """
#
#
#rule clip_overlap:
#  input:
#    deDup = 'deDup/{sample}.dedup.bam'
#  output:
#    clip = 'deDup/{sample}.overlapclipped.bam' 
#  log: 'log/{sample}.overlapclipped.bam.log'
#  threads: 12
#  message:
#    """ Clip overlapping paired end reads """
#  shell:
#    """
#    bam clipOverlap --in {input.deDup} --out {output.clip} --stats --poolSize 5000000 2> {log}
#    """

rule sym_link_REFs:
  input:
    '/home/uibk/c7701178/scratch/DAPHNIA/daphnia_snakemake_pbs/deDup/HiC/minid95/{sample}.overlapclipped.bam'
  output:
    touch('deDup/{sample}.symlink.done')
  log: 'log/{sample}.symlink.log'
  message: """--- Create a symlink for clipped REF files ---"""
  shell:
    """
    ln -s {input} deDup/{wildcards.sample}.overlapclipped.bam 2> {log}
    """

#
#rule Index_clippedBAM:
#  input:
#    clip = 'deDup/{sample}.overlapclipped.bam'
#  output:
#    idx = 'deDup/{sample}.overlapclipped.bam.bai'
#  log: 'log/{sample}.overlapclipped.bam.log'
#  threads: 2
#  message: """--- Indexing clipped BAM files with samtools ---"""
#  shell:
#    """
#    samtools index {input.clip} {output.idx} 2> {log}
#    """
#


rule refIndex:
  input:
    ref = config['ref_HiC']
  output:
    "ref/{ref_name}.fasta.fai"
  log: 'log/{ref_name}_refIndex.log'
  shell:
    """
    samtools faidx {input.ref} 2> {log}
    """


rule ref_Dict:
  input:
    ref = config['ref_HiC']
  output:
    'ref/{ref_name}.dict'
  log: 'log/{ref_name}_refDict.log'
  shell:
    """
    /apps/uibk/bin/sysconfcpus -n 12 java -Xmx110g -jar /home/uibk/c7701178/.conda/envs/eggs2.0/share/picard-2.25.0-1/picard.jar CreateSequenceDictionary R={input.ref} O={output} 2> {log}
    """

    
rule ls_ClipBam:
  input:
    clip = 'deDup/{sample}.overlapclipped.bam'
  output:
    touch('deDup/{sample}.added2ClippedList.done')
  log: 'log/{sample}.added2ClippedList.log'
  message: """--- Creating a sample list of clipped bam files for indel-realignement ---"""
  shell:
    """
    ls {input.clip} >> list/overlapclippedBAM.list 2> {log}
    """


rule ls_ClipBam_REFs:
  input:
    clip = '/home/uibk/c7701178/scratch/DAPHNIA/daphnia_snakemake_pbs/deDup/HiC/minid95/{sample}.overlapclipped.bam'
  output:
    touch('deDup/{sample}.REFsadded2ClippedList.done')
  log: 'log/{sample}.REFsadded2ClippedList.log'
  message: """--- Creating a sample list of clipped bam files for indel-realignement ---"""
  shell:
    """
    ls {input.clip} >> list/overlapclippedBAM.list 2> {log}
    """


rule list_indels:
  input:
    clip = 'list/overlapclippedBAM_curvi.list',
    ref = config['ref_HiC'],
    idx_ref = "ref/{ref_name}.fasta.fai",
    dict_ref = 'ref/{ref_name}.dict'
  output:
    'list/indels_{ref_name}_curvi.list'
  log: 'log/listIndels{ref_name}_curvi.log'
  threads: 48
  message:
    """ Create list of potential indels """
  shell:
    """
    GATK=(/home/uibk/c7701178/.conda/envs/eggs2.0/opt/gatk-3.8/GenomeAnalysisTK.jar)
    /apps/uibk/bin/sysconfcpus -n 48 java -Xmx440g -jar $GATK -T RealignerTargetCreator -R {input.ref} -I {input.clip} -o {output} -drf BadMate 2> {log}
    """



rule realign_indel:
  input:
    clip = 'deDup/{sample}.overlapclipped.bam',
    ref = config['ref_HiC'],
    indels = 'list/indels_Dgaleata_M5_PBasm.FINAL_curvi.list'
    #indels = 'list/indels_Dgaleata_M5_PBasm.FINAL.list'
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


rule samtools_coverage:
  input:
    realigned = 'realigned/{sample}.realigned.bam'
  output:
    coverage = 'depth/{sample}.realigned.bam.coverage.hist'
  log: 'log/{sample}.realigned.bam.coverage.hist.log'
  threads: 12
  message: """ Produce an ASCII-art histogram of coverage per chromosome   """
  shell:
    """
    samtools coverage -A -o {output.coverage} {input.realigned} 2> {log}
    """


rule samtools_depth:
  input:
    realigned = 'realigned/{sample}.realigned.bam'
  output:
    depth = 'depth/{sample}.realigned.bam.depth.gz'
  log: 'log/{sample}.realigned.bam.depth.log'
  threads: 12
  message: """ Count per position depth per sample using samtools depth """
  shell:
    """
    samtools depth -aa {input.realigned} | cut -f3 | gzip > {output.depth} 2> {log}
    """


rule ls_depth:
  input:
    'depth/{sample}.realigned.bam.depth.gz'
  output:
    touch('depth/{sample}.added2DepthList.list')
  log: 'log/{sample}.depth.list.log'
  message: """ --- Creating a sample list of samtools-depth files for Rscript --- """
  shell:
    """
    ls {input} | cut -f2 -d '/' >> depth/depth_curvi.list 2> {log}
    """


rule read_depth:
  input:
    'depth/depth_curvi.list'
  output:
    'depth/stats/depth_statistics_curvi.txt'
  log: 'log/depth_statistics_curvi.log'
  threads: 12
  message:
    """ --- Running Rscript to plot the genome-wide distribution of coverage --- """
  shell:
    """
    Rscript scripts/read_depth.R {input} {output} 2> {log} 
    """


rule plot_summary:
  input:
    'depth/stats/depth_statistics_{sets}.txt'
  output:
    touch('depth/plots/plot_summary_{sets}.done')
  log: 'log/plot_summary_{sets}log'
  message:
    """ --- Running Rscript to plot depth summary per sample and output depth filters --- """
  shell:
    """
    Rscript scripts/plot_summary.R {input} {wildcards.sets} 2> {log}
    """


rule rbind_depthFilter:
  output:
    touch('depth/stats/rbind_depthFilter.done')
  log: 'log/rbind_dfFilter.log'
  message:
    """ --- Running Rscript to plot depth summary per sample and output depth filters --- """
  shell:
    """
    Rscript scripts/rbind_dfFilter.R depthFilter.list 2> {log}
    """

rule genome_coverage_bed:
  input:
    'realigned/{sample}.realigned.bam'
  output:
    'bedtools/{sample}.realigned.genomecov.bed'
  log: 'log/{sample}.realigned.genomecov.log'
  threads: 12
  message:
    """ Computes BED summaries using bedtools """
  shell:
    """
    bedtools genomecov -ibam {input} > {output} 2> {log}
    """


rule plot_gencov:
  input:
    bed= 'bedtools/{sample}.realigned.genomecov.bed'
  output:
    pdf = 'bedtools/plots/{sample}.realigned.genomecov.pdf'
  log: 'log/{sample}.realigned.genomecov_plot.log'
  threads: 4
  message:
    """ Running Rscript to plot the genome-wide distribution of coverage """
  shell:
    """
    Rscript scripts/plot_gene_covs.R {input.bed} {output.pdf} 2> {log}
    """
