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
		#expand('trm/qc/multiqc/')#,
		## !!! works, but some intermediate files from Kraken were deleted due to space limitation on mach2 !!! ##
		#expand('KRAKEN2_RESULTS/{sample}_kraken2_classified.done', sample=sample_names),
		#expand('KRAKEN2_RESULTS/{sample}_{pair}.{ext}', sample=sample_names, pair=['R1', 'R2'], ext=['trmdfilt.classified_retain.fq']),
		#expand('KRAKEN2_RESULTS/{sample}_{pair}.{ext}', sample=sample_names, pair=['R1','R2'], ext=['trmdfilt.keep.fq.gz']),
		#expand('KRAKEN2_RESULTS/qc/fastqc/{sample}_{pair}.trmdfilt.keep_fastqc.{ext}', sample=sample_names, pair=['R1', 'R2'], ext=['zip', 'html']),
		#expand('KRAKEN2_RESULTS/qc/multiqc/'),
		#expand('ref/genome/1/summary.txt'),
		#expand('bbmap/rapid/{sample}.{ext}', sample=sample_names, ext=['bam']),
		#expand('bbmap/rapid/{sample}.{ext}', sample=sample_names, ext=['bam.bai']),
		#expand('deDup/{sample}.{ext}', sample=sample_names, ext=['dedup.bam', 'dedup.metrics.txt'])#,
		#expand('deDup/{sample}.{ext}', sample=refclone_names, ext=['dedup.bam']),
		#expand('deDup/{sample}.{ext}', sample=refclone_names, ext=['dedup.bam.bai']),
		#expand('deDup/{sample}.{ext}', sample=ALL_names, ext=['overlapclipped.bam']),
		#expand('deDup/{sample}.{ext}', sample=ALL_names, ext=['overlapclipped.bam.bai']),
		#expand('ref/dgal_ra_pb-target-and-other_ill-confilter.blobfilter.rmmt_sspace-lr3_lrgc3_pg3_pilon_3.simple.{ext}', ext=['fasta.fai', 'dict']),
		#expand('deDup/{sample}.added2ClippedList.done', sample=ALL_names)#,
		#expand('list/indels.list'),
		#expand('realigned/{sample}.{ext}', sample=ALL_names, ext=['realigned.bam']),
		#expand('depth/{sample}.{ext}', sample=ALL_names, ext=['realigned.bam.depth.gz']),
		#expand('depth/{sample}.added2DepthList.list', sample=ALL_names),
		#expand('depth/stats/depth_statistics_{sets}.txt', sets=['curvi']),
		#expand('depth/plots/plot_summary_{sets}.done', sets=['curvi']),
		#expand('depth/stats/depthFilter.list'),
		#expand('list/new_metadata_table_{sets}.done', sets=sets)#,
		#expand('bedtools/{sample}.realigned.{ext}', sample=ALL_names, ext=['genomecov.bed']),
		#expand('bedtools/plots/{sample}.realigned.{ext}', sample=ALL_names, ext=['genomecov.pdf']),
		#expand('angsd/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done', zip, sets=sets, GL=GL, minMaf=minMaf, IND=N, MinDepth=MinDepth, MaxDepth=MaxDepth),
		#expand('pcangsd/{sets}/PCAngsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done', zip, sets=sets, GL=GL, minMaf=minMaf, IND=N, MinDepth=MinDepth, MaxDepth=MaxDepth),
		#expand('pcangsd/{sets}/PCAngsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_covmat.pdf', zip, sets=sets, GL=GL, minMaf=minMaf, IND=N, MinDepth=MinDepth, MaxDepth=MaxDepth),
                #expand('pcangsd/{sets}/PCAngsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_unpruned_admix_K{K}.done', sets=['LC_REF'], GL=['2'], IND=['221'], minMaf=['0.05'], MinDepth=['221'], MaxDepth=['10006'], K=admix_K),
		#expand('pcangsd/{sets}/PCAngsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_unpruned_admix_K{K}.done', sets=['LC_withoutREF'], GL=['2'], IND=['177'], minMaf=['0.05'], MinDepth=['177'], MaxDepth=['8091'], K=admix_K),
		#expand('ngsLD/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_sub_geno.beagle.gz', zip, sets=sets, GL=GL, minMaf=minMaf, IND=N, MinDepth=MinDepth, MaxDepth=MaxDepth),
		#expand('ngsLD/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_sub_pos.gz', zip, sets=sets, GL=GL, minMaf=minMaf, IND=N, MinDepth=MinDepth, MaxDepth=MaxDepth),
		#expand('ngsLD/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_sub.ld.gz', zip, sets=sets, GL=GL, minMaf=minMaf, IND=N, MinDepth=MinDepth, MaxDepth=MaxDepth),
                #expand('ngsLD/{sets}/run_LDpruning_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done', zip, sets=sets, GL=GL, minMaf=minMaf, IND=N, MinDepth=MinDepth, MaxDepth=MaxDepth),
                #expand('ngsLD/{sets}/LDpruned_snps_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.list', zip, sets=sets, GL=GL, minMaf=minMaf, IND=N, MinDepth=MinDepth, MaxDepth=MaxDepth),
		#expand('ngsLD/{sets}/index_SNPs_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done', zip, sets=sets, GL=GL, minMaf=minMaf, IND=N, MinDepth=MinDepth, MaxDepth=MaxDepth),
		#expand('ngsLD/{sets}/LDpruned_snps_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.chr', zip, sets=sets, GL=GL, minMaf=minMaf, IND=N, MinDepth=MinDepth, MaxDepth=MaxDepth),
		#expand('angsd/{sets}/LDpruned_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done', zip, sets=sets, GL=GL, minMaf=minMaf, IND=N, MinDepth=MinDepth, MaxDepth=MaxDepth),
		#expand('pcangsd/{sets}/PCAngsd_LDpruned_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done', zip, sets=sets, GL=GL, minMaf=minMaf, IND=N, MinDepth=MinDepth, MaxDepth=MaxDepth),
		#expand('pcangsd/{sets}/PCAngsd_LDpruned_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_admix_K{K}.done', sets=['LC_REF'], GL=['2'], IND=['221'], minMaf=['0.05'], MinDepth=['221'], MaxDepth=['10006'], K=admix_K),
		#expand('pcangsd/{sets}/PCAngsd_LDpruned_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_admix_K{K}.done', sets=['LC_withoutREF'], GL=['2'], IND=['177'], minMaf=['0.05'], MinDepth=['177'], MaxDepth=['8091'], K=admix_K),
                #expand('saf/saf_samples/{sample}.GL2.saf.idx.done', sample=saf_prefix)#,
                #expand('saf/saf_samples/{sample}.GL2.est.ml', sample=saf_prefix),
                expand('saf/POPS/df_{POPS}.{ext}', POPS=['REF_long', 'REF_gal', 'REF_cuc', 'PRE_long', 'POST_long_LC'], ext=['GL2.saf.idx.done']),
		expand('saf/POPS/df_{POPS}.{ext}', POPS=['REF_long', 'REF_gal', 'REF_cuc', 'PRE_long', 'POST_long_LC'], ext=['GL2.sfs']),
                expand('saf/POPS/df_{POPS}.{ext}', POPS=['REF_long', 'REF_gal', 'REF_cuc', 'PRE_long', 'POST_long_LC'], ext=['sfs2theta.done']),
                expand('saf/POPS/df_{POPS}.{ext}', POPS=['REF_long', 'REF_gal', 'REF_cuc', 'PRE_long', 'POST_long_LC'], ext=['theta_stats_chrom.done']),
                expand('saf/POPS/df_{POPS}.{ext}', POPS=['REF_long', 'REF_gal', 'REF_cuc', 'PRE_long', 'POST_long_LC'], ext=['theta.thetasWindow10kb.gz']),
                expand('saf/POPS/{POP1}_vs_{POP2}_2D.sfs', zip, POP1=POP1, POP2=POP2),
                expand('saf/POPS/{POP1}_vs_{POP2}.FstIndex.done', zip, POP1=POP1, POP2=POP2),
                expand('saf/POPS/{POP1}_vs_{POP2}.Fst_Global.done', zip, POP1=POP1, POP2=POP2)
                #expand('saf/{sample}.heterozygosity.txt', sample=sample_prefix)
                #expand('angsd/{sample}.GL2.est.ml', sample=BAM),
                #expand('test_angsd935.done')
    
	
# -----------------------------------------------


include: "rules/eggs2.0.smk"
#include: "rules/get_beagle_allSNPs.smk"
#include: "rules/PCAngsd_allSNPs.smk"
#include: "rules/ngsLD.smk"
#include: "rules/get_beagle_LDprunedSNPs.smk"
#include: "rules/PCAngsd_LDprunedSNPs.smk"
include: "rules/realSFS.smk"
#include: "rules/gencov.sm"
