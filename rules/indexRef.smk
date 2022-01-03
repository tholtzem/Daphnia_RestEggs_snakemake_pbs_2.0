
rule refIndex:
  input:
    ref = config['ref']
  output:
    'ref/dgal_ra_pb-target-and-other_ill-confilter.blobfilter.rmmt_sspace-lr3_lrgc3_pg3_pilon_3.simple.fasta.fai'
  log: 'log/refIndex.log'
  shell:
    """
    samtools faidx {input.ref} 2> {log}
    """

rule ref_Dict:
  input:
    ref = config['ref']
  output:
    'ref/dgal_ra_pb-target-and-other_ill-confilter.blobfilter.rmmt_sspace-lr3_lrgc3_pg3_pilon_3.simple.dict'
  log: 'log/refDict.log'
  shell:
    """
    /apps/uibk/bin/sysconfcpus -n 12 java -Xmx110g -jar /home/uibk/c7701178/.conda/envs/eggs/share/picard-2.25.0-1/picard.jar CreateSequenceDictionary R={input.ref} O={output} 2> {log}
    """
