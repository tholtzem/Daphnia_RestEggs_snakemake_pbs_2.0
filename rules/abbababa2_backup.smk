# Here, we use the global list of SNPs, we called/calculated from GLs of the entire D.lgc-complex (including reference panel, resting eggs and pelagial). This set should explain the variation in our data from the D.lgc-complex. Then, we add the outgroup to the data set, and call SNPs from GLs, but only for those sites present in the D. lgc-complex.


#rule do_fasta:
#  input:
#    'realigned/{sample}.realigned.bam'
#  output:
#    'abbababa/perfect.sample{sample}.fa.fai'
#  log:
#    'log/abbababa/perfect.sample{sample}.fa.fai.log'
#  threads: 24
#  message:
#    """ Generate a fasta file for one of the bam file in our data. We assume such a genome has very high quality and we can use it as a reference for estimating error rates in others of our data. """
#  shell:
#    """
#    module load angsd/0.938
#    /apps/uibk/bin/sysconfcpus -n 24 angsd -i {input} -doFasta 1 -doCounts 1 -out abbababa/perfectSample{sample} &&
#    gunzip abbababa/perfectSample{sample}.fa.gz &&
#    samtools faidx perfect.sample{sample}.fa 2> {log}
#    """


#rule AncError:
#  input:
#    'list/pop_list/df_PRE_long.txt' 
#  output:
#    'abbababa/perfect.sample{sample}.fa.fai'
#  log:
#    'log/abbababa/perfect.sample{sample}.fa.fai.log'
#  threads: 12
#  message:
#    """ We apply error correction to the PRE-eutrophication group (which might be affected by transition error) using "perfectSampleCEU" as high-quality reference genome. """
#  shell:
#    """
#    module load angsd/0.938
#    /apps/uibk/bin/sysconfcpus -n 12 angsd -doAncError 1 -anc chimpHg19.fa -ref perfectSampleCEU.fa -out errorFile -bam bamWithErrors.filelist 2> {log}
#    """

localrules: all, prepare_bamlist, prepare_sizefile, create_errorFile, get_popNames


rule prepare_bamlist:
  input:
    'list/pop_list/4pop_combinations.list'
  output:
    touch('list/abbababa/4pop_combinations.bam.list.done')
  log:
    'log/abbababa/4pop_combinations.bam.list.log'
  threads: 12
  message:
    """ Prepare list of bam files required for the ABBA-BABA test (D-statistic) """
  shell:
    """
    Lines=$(cat {input} | tail -n +2)
    path2in=$(echo {input} | cut -f1,2 -d'/')
    path2out=list/abbababa
    
    for Line in $Lines; do
      p1=$(echo $Line | cut -f1 -d',')
      p2=$(echo $Line | cut -f2 -d',')
      p3=$(echo $Line | cut -f3 -d',')
      p4=$(echo $Line | cut -f4 -d',')
      echo $p1 $p2 $p3 $p4
      cat $path2in/${{p1}}.txt $path2in/${{p2}}.txt $path2in/${{p3}}.txt $path2in/${{p4}}.txt > $path2out/${{p1}}_${{p2}}_${{p3}}_${{p4}}.bam.list 2> {log}
    done
    """


rule prepare_sizefile:
  input:
    'list/pop_list/4pop_combinations.list'
  output:
    touch('list/abbababa/4pop_combinations.sizefile.done')
  log:
    'log/abbababa/4pop_combinations.sizefile.log'
  threads: 12
  message:
    """ Prepare list of bam files required for the ABBA-BABA test (D-statistic) """
  shell:
    """
    Lines=$(cat {input} | tail -n +2)
    path2in=$(echo {input} | cut -f1,2 -d'/')
    path2out=list/abbababa
    
    for Line in $Lines; do
      p1=$(echo $Line | cut -f1 -d',')
      p2=$(echo $Line | cut -f2 -d',')
      p3=$(echo $Line | cut -f3 -d',')
      p4=$(echo $Line | cut -f4 -d',')
      echo $p1 $p2 $p3 $p4

      NBR_pop1=$(cat $path2in/${{p1}}.txt | wc -l)
      NBR_pop2=$(cat $path2in/${{p2}}.txt | wc -l)
      NBR_pop3=$(cat $path2in/${{p3}}.txt | wc -l)
      NBR_pop4=$(cat $path2in/${{p4}}.txt | wc -l)
      printf '%s\n' $NBR_pop1 $NBR_pop2 $NBR_pop3 $NBR_pop4 > $path2out/${{p1}}_${{p2}}_${{p3}}_${{p4}}.sizefile.list
    done 2> {log}
    """



rule ABBA_BABA:
  input:
    'list/abbababa/4pop_combinations.bam.list.done',
    'list/abbababa/4pop_combinations.sizefile.done',
    ref = config["ref_rapid"],
    #bamlist = 'list/abbababa/{combi}.bam.list',
    #popsize = 'list/abbababa/{combi}.sizefile.list',
    #sites = 'ngsLD/LC_REF/LDpruned_snps_angsd_GL2_minInd221_maf0.05_minDepth221_maxDepth10006.list',
    sites = 'angsd/LC_REF/angsd_GL2_minInd221_maf0.05_minDepth221_maxDepth10006_globalSNP.list',
    chroms = 'angsd/LC_REF/angsd_GL2_minInd221_maf0.05_minDepth221_maxDepth10006_globalSNP.chr'
  output:
   touch('abbababa/{combi}.abbababa.done')
  log:
    'log/abbababa/{combi}.abbababa.log'
  threads: 24
  message:
    """ Compute ABBA-BABA test (D-statistic) using angsd, allows for multiple individuals in each group """
  shell:
    """
    bamlist=(list/abbababa/{wildcards.combi}.bam.list)
    popsize=(list/abbababa/{wildcards.combi}.sizefile.list)

    module load angsd/0.938
    /apps/uibk/bin/sysconfcpus -n 24 angsd -doAbbababa2 1 -b $bamlist -sizeFile $popsize -useLast 1 -doCounts 1 -out abbababa/abbababa_{wildcards.combi} -blockSize 1000000 -ref {input.ref} -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minQ 20 -minMapQ 30 -rmTrans 0 -sites {input.sites} -rf {input.chroms} 2> {log}
    """


rule create_errorFile:
  output:
    'list/abbababa/errorFile.txt'
  log: 'log/abbababa/errorFile.log'
  message:
    """ Define the error file in text file (if available, if not write NA) """
  shell:
    """
    printf '%s\n' NA NA NA NA > {output} 2> {log}
    """


rule get_popNames:
  input:
    'list/pop_list/4pop_combinations.list'
  output:
    touch('list/abbababa/4pop_combinations.Popnames.list.done')
  log:
    'log/abbababa/4pop_combinations.Popnames.list.log'
  threads: 12
  message:
    """ Prepare list of population names required for the R script calculating D statistic """
  shell:
    """
    Lines=$(cat {input} | tail -n +2)
    path2out=list/abbababa
    
    for Line in $Lines; do
      p1=$(echo $Line | cut -f1 -d',')
      p2=$(echo $Line | cut -f2 -d',')
      p3=$(echo $Line | cut -f3 -d',')
      p4=$(echo $Line | cut -f4 -d',')
      printf '%s\n' ${{p1}} ${{p2}} ${{p3}} ${{p4}} > $path2out/${{p1}}_${{p2}}_${{p3}}_${{p4}}.POPnames.list 2> {log}
    done
    """


rule get_Dstats:
  input:
    touched = 'list/abbababa/4pop_combinations.Popnames.list.done',
    errFile = 'list/abbababa/errorFile.txt'
  output:
    touch('abbababa/{combi}.Dstat_results.done')
  log:
    'log/abbababa/{combi}.Dstat_results.log'
  threads: 12
  message: """ Run R script downloaded from https://github.com/ANGSD/angsd/blob/master/R/estAvgError.R """
  shell:
    """
    pop_list=(list/abbababa/{wildcards.combi}.POPnames.list)
    
    Rscript scripts/estAvgError.R angsdFile=abbababa/abbababa_{wildcards.combi} out=abbababa/abbababa_{wildcards.combi}.Dstatresults sizeFile=list/abbababa/{wildcards.combi}.sizefile.list errFile={input.errFile} nameFile=$pop_list 2> {log}
    """



