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
		#expand('ref/bb_indexRef.done'),
		#expand('bbmap/{sample}.{ext}', sample=sample_names, ext=['bam']),
		#expand('bbmap/{sample}.{ext}', sample=sample_names, ext=['bam.bai']),
		#expand('deDup/{sample}.{ext}', sample=sample_names, ext=['dedup.bam', 'dedup.metrics.txt']),
		#expand('deDup/{sample}.{ext}', sample=sample_names, ext=['overlapclipped.bam']),
		#expand('deDup/{sample}.{ext}', sample=sample_names, ext=['overlapclipped.bam.bai']),
		#expand('ref/{ref_name}.{ext}', ref_name=['Dgaleata_M5_PBasm.FINAL'] ,ext=['fasta.fai', 'dict']),
		#expand('deDup/{sample}.added2ClippedList.done', sample=sample_names),
		#expand('deDup/{sample}.REFsadded2ClippedList.done', sample=refclone_names),
		##Run until here
		#expand('list/indels_{ref_name}.list', ref_name=['Dgaleata_M5_PBasm.FINAL']),
                #expand('deDup/{sample}.symlink.done', sample=refclone_names),
		#expand('realigned/{sample}.{ext}', sample=Dlgc_names, ext=['realigned.bam']),
                #expand('depth/{sample}.{ext}', sample=Dlgc_names, ext=['realigned.bam.coverage.hist']),
		#expand('depth/{sample}.{ext}', sample=Dlgc_names, ext=['realigned.bam.depth.gz']),
		#expand('depth/{sample}.added2DepthList.list', sample=Dlgc_names),
		#expand('depth/stats/depth_statistics.txt'),
		#expand('depth/plots/plot_summary_{sets}.done', sets=['LC', 'LZ', 'LCwithoutREF', 'LZwithoutREF', 'longREF7', 'longREF9', 'cucREF7', 'galREF9', 'ALL_gal_LC', 'hatched_gal_LC', 'EGGS_gal_LC', 'PEL_cuc_LC', 'PEL_long_LC', 'PRE_long_LC']),
                ## run until here
		#expand('depth/stats/rbind_depthFilter.done'),
		#expand('bedtools/{sample}.realigned.{ext}', sample=ALL_names, ext=['genomecov.bed']),
		#expand('bedtools/plots/{sample}.realigned.{ext}', sample=ALL_names, ext=['genomecov.pdf']),
		expand('angsd/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done', zip, sets=sets_structure[0], GL=GL_structure[0], minMaf=minMaf[0], IND=N_structure[0], MinDepth=MinDepth_structure[0], MaxDepth=MaxDepth_structure[0]),
		expand('angsd/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_globalSNP.list', zip, sets=sets_structure[0], GL=GL_structure[0], minMaf=minMaf[0], IND=N_structure[0], MinDepth=MinDepth_structure[0], MaxDepth=MaxDepth_structure[0]),
		expand('angsd/{sets}/index_globalSNP_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done', zip, sets=sets_structure[0], GL=GL_structure[0], minMaf=minMaf[0], IND=N_structure[0], MinDepth=MinDepth_structure[0], MaxDepth=MaxDepth_structure[0]),
		expand('angsd/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_globalSNP.chr', zip, sets=sets_structure[0], GL=GL_structure[0], minMaf=minMaf[0], IND=N_structure[0], MinDepth=MinDepth_structure[0], MaxDepth=MaxDepth_structure[0]),
		expand('pcangsd/{sets}/PCAngsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_cov_admix.done', zip, sets=sets_structure[0], GL=GL_structure[0], minMaf=minMaf[0], IND=N_structure[0], MinDepth=MinDepth_structure[0], MaxDepth=MaxDepth_structure[0]),
		expand('ngsLD/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_sub_geno.beagle.gz', zip, sets=sets_structure[0], GL=GL_structure[0], minMaf=minMaf[0], IND=N_structure[0], MinDepth=MinDepth_structure[0], MaxDepth=MaxDepth_structure[0]),
		expand('ngsLD/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_sub_pos.gz', zip, sets=sets_structure[0], GL=GL_structure[0], minMaf=minMaf[0], IND=N_structure[0], MinDepth=MinDepth_structure[0], MaxDepth=MaxDepth_structure[0]),
		expand('ngsLD/{sets}/run_ngsLD_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_sub.ld.gz', zip, sets=sets_structure[0], GL=GL_structure[0], minMaf=minMaf[0], IND=N_structure[0], MinDepth=MinDepth_structure[0], MaxDepth=MaxDepth_structure[0]),
		expand('ngsLD/{sets}/run_LDpruning_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done', zip, sets=sets_structure[0], GL=GL_structure[0], minMaf=minMaf[0], IND=N_structure[0], MinDepth=MinDepth_structure[0], MaxDepth=MaxDepth_structure[0]),
		expand('ngsLD/{sets}/LDpruned_snps_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.list', zip, sets=sets_structure[0], GL=GL_structure[0], minMaf=minMaf[0], IND=N_structure[0], MinDepth=MinDepth_structure[0], MaxDepth=MaxDepth_structure[0]),
		expand('ngsLD/{sets}/index_SNPs_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done', zip, sets=sets_structure[0], GL=GL_structure[0], minMaf=minMaf[0], IND=N_structure[0], MinDepth=MinDepth_structure[0], MaxDepth=MaxDepth_structure[0]),
		expand('ngsLD/{sets}/LDpruned_snps_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.chr', zip, sets=sets_structure[0], GL=GL_structure[0], minMaf=minMaf[0], IND=N_structure[0], MinDepth=MinDepth_structure[0], MaxDepth=MaxDepth_structure[0]),
		expand('angsd/{sets}/LDpruned_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done', zip, sets=sets_structure[0], GL=GL_structure[0], minMaf=minMaf[0], IND=N_structure[0], MinDepth=MinDepth_structure[0], MaxDepth=MaxDepth_structure[0]),
		expand('pcangsd/{sets}/PCAngsd_LDpruned_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done', zip, sets=sets_structure[0], GL=GL_structure[0], minMaf=minMaf[0], IND=N_structure[0], MinDepth=MinDepth_structure[0], MaxDepth=MaxDepth_structure[0]),
                #expand('saf/saf_samples/{sample}.GL2.saf.idx.done', sample=ALL_names),
                #expand('saf/saf_samples/{sample}.GL2.est.ml', sample=ALL_names),
                #expand('saf/POPS/df_{POPS}.{ext}', POPS=['REF_long', 'REF_gal', 'REF_cuc', 'PRE_long', 'POST_long_LC'], ext=['GL2.saf.idx.done']),
		#expand('saf/POPS/df_{POPS}.{ext}', POPS=['REF_long', 'REF_gal', 'REF_cuc', 'PRE_long', 'POST_long_LC'], ext=['GL2.sfs']),
                #expand('saf/POPS/df_{POPS}.{ext}', POPS=['REF_long', 'REF_gal', 'REF_cuc', 'PRE_long', 'POST_long_LC'], ext=['sfs2theta.done']),
                #expand('saf/POPS/df_{POPS}.{ext}', POPS=['REF_long', 'REF_gal', 'REF_cuc', 'PRE_long', 'POST_long_LC'], ext=['theta_stats_chrom.done']),
                #expand('saf/POPS/df_{POPS}.{ext}', POPS=['REF_long', 'REF_gal', 'REF_cuc', 'PRE_long', 'POST_long_LC'], ext=['theta.thetasWindow10kb.gz.pestPG']),
                #expand('saf/POPS/{POP1}_vs_{POP2}_2D.sfs', zip, POP1=POP1, POP2=POP2),
                #expand('saf/POPS/{POP1}_vs_{POP2}.FstIndex.done', zip, POP1=POP1, POP2=POP2),
                #expand('saf/POPS/{POP1}_vs_{POP2}.Fst_Global.done', zip, POP1=POP1, POP2=POP2),
                #expand('list/abbababa/4pop_combinations.bam.list.done'),
                #expand('list/abbababa/4pop_combinations.sizefile.done'),
                #expand('abbababa/{combi}.abbababa.done', combi=Dstats_combi),
                #expand('list/abbababa/errorFile.txt'),
                #expand('list/abbababa/4pop_combinations.Popnames.list.done'),
                #expand('abbababa/{combi}.Dstat_results.done', combi=Dstats_combi),
                #expand('ancestry/get_bcf.done'),
                #expand('ancestry/angsd_GL2_minInd221_maf0.05_minDepth221_maxDepth10006.{ext}', ext=['vcf.gz']),
		#'list/parents_hybrids.list',
		#expand('ancestry/parents_hybrids.{ext}', ext=['vcf.gz']),
		#expand('list/pop_list/{pop}.csv', pop=['longispina', 'galeata', 'cucullata', 'hybrids_LC']),
                #expand('ancestry/fixed_sites.{ext}', ext=['txt']),
                #expand('ancestry/fixed_sites.{ext}', ext=['svg']),
                #'ancestry/index_fixed_sites.done',
		#expand('ancestry/{sample}.GL2.saf.idx.done', sample=bams),
		#expand('ancestry/{sample}.GL2.est.ml', sample=bams)
# -----------------------------------------------


#include: "rules/eggs2.0.smk"
#include: "rules/get_vcf.smk"
#include: "rules/ancestry.smk"
include: "rules/get_beagle_allSNPs.smk"
include: "rules/PCAngsd_allSNPs.smk"
include: "rules/ngsLD.smk"
include: "rules/get_beagle_LDprunedSNPs.smk"
include: "rules/PCAngsd_LDprunedSNPs.smk"
#include: "rules/realSFS.smk"
#include: "rules/abbababa2.smk"
#include: "rules/gencov.sm"
