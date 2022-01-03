include: "rules/common.smk"



# -----------------------------------------------

rule all:
	input:
		#expand('raw/qc/fastqc/{sample}_{pair}_001_fastqc.{ext}', sample=sample_names, pair=['R1', 'R2'], ext=['html', 'zip']),
		#expand('raw/qc/multiqc/'),
		## !!! works, but some intermediate files from trimming were deleted due to space limitation on mach2 !!! ##
		#expand('trm/{sample}_{pair}.{ext}', sample=sample_names, pair=['R1','R2'], ext=['trmd.fq.gz']),
		#expand('trm/{sample}_{pair}.{ext}', sample=sample_names, pair=['R1','R2'], ext=['trmdfilt.fq.gz']),
		#expand('trm/qc/fastqc/{sample}_{pair}.trmdfilt_fastqc.{ext}', sample=sample_names, pair=['R1', 'R2'], ext=['html', 'zip']),
		#expand('trm/qc/multiqc/'),
		## !!! works, but some intermediate files from Kraken were deleted due to space limitation on mach2 !!! ##
		#expand('KRAKEN2_RESULTS/{sample}_kraken2_classified.done', sample=sample_names),
		#expand('KRAKEN2_RESULTS/{sample}_{pair}.{ext}', sample=samples183_names, pair=['R1', 'R2'], ext=['trmdfilt.classified_retain.fq']),
		#expand('KRAKEN2_RESULTS/{sample}_{pair}.{ext}', sample=samples183_names, pair=['R1','R2'], ext=['trmdfilt.keep.fq.gz']),
		#expand('KRAKEN2_RESULTS/qc/fastqc/{sample}_{pair}.trmdfilt.keep_fastqc.{ext}', sample=samples183_names, pair=['R1', 'R2'], ext=['zip', 'html']),
		#expand('KRAKEN2_RESULTS/qc/multiqc/'),
		expand('ref/genome/1/summary.txt'),
		expand('bbmap/rapid/{sample}.{ext}', sample=samples183_names, ext=['bam']),
		expand('bbmap/rapid/{sample}.{ext}', sample=samples183_names, ext=['bam.bai']),
		expand('deDup/{sample}.{ext}', sample=samples183_names, ext=['dedup.bam', 'dedup.metrics.txt']),
                expand('deDup/{sample}.{ext}', sample=refclone_names, ext=['dedup.bam'])
		#expand('deDup/{sample}.{ext}', sample=samples253_names, ext=['overlapclipped.bam'])
		#expand('deDup/{sample}.{ext}', sample=samples253_names, ext=['overlapclipped.bam.bai']),
		#expand('ref/dgal_ra_pb-target-and-other_ill-confilter.blobfilter.rmmt_sspace-lr3_lrgc3_pg3_pilon_3.simple.fasta.{ext}', ext=['fai', 'dict'])#,
#-----------------------------------------------------------------------------------------------------------
		#expand('list/overlapclippedBAM.list'),
		#expand('list/all_samples.{ext}', ext=['intervals']),
		#expand('realigned/{sample}.{ext}', sample=new_sample_names, ext=['minq20.realigned.bam']),
		#expand('depth/{sample}.{ext}', sample=new_sample_names, ext=['minq20.realigned.bam.depth.gz']),
		#expand('bedtools/{sample}.{ext}', sample=new_sample_names, ext=['minq20.realigned.genomecov.bed']),
		#expand('bedtools/plots/{sample}.{ext}', sample=new_sample_names, ext=['minq20.realigned.genomecov.pdf']),
		#expand('genome_stats.done'),
		##expand('angsd/angsd_GL2_cutoffs_maf0018_nInd55.done'),
		##expand('pcangsd/PCAngsd_GL2_cutoffs_maf0018_nInd55.done'),
		##expand('pcangsd/PCAngsd_GL2_cutoffs_maf0018_nInd55_covmat.pdf'),
		#expand('ngsLD/angsd_LC_GL2_cutoff_maf0018_nInd55_sub_geno.beagle.gz'),
		##expand('ngsLD/angsd_LC_GL2_cutoff_maf0018_nInd55_sub_pos.gz'),
		##expand('ngsLD/angsd_LC_GL2_cutoff_maf0018_nInd55_sub.ld.gz'),
                #expand('ngsLD/angsd_LC_GL2_cutoff_maf0018_nInd55_sub.unlinked.id'),
                #expand('LDpruned_SNPlist.done'),
                #expand('list/index_SNP.done'),
                #expand('list/getChrom_from_sites.done'),
                #expand('angsd/angsd_LC_GL2_maf0018_LDpruned.done'),
                #expand('pcangsd/PCAngsd_GL2_maf0018_LDpruned.done'),
                #expand('pcangsd/PCAngsd_GL2_maf0018_LDpruned.pdf'),
                #expand('pcangsd/PCAngsd_GL2_maf0018_LDpruned_admix_K3.done'),
                #expand('angsd/{sample}.GL2.saf.idx.done', sample=BAM),
                #expand('angsd/{sample}.GL2.est.ml', sample=BAM),
                #expand('test_angsd935.done')
    
	
# -----------------------------------------------


include: "rules/eggs2.0.smk"

