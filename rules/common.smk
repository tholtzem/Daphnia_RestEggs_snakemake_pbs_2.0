import os
import pandas as pd

###### Config file and sample sheets #####


configfile: 
  "config/config.yaml"

# load sample info from raw data
samples003 = pd.read_csv("samples003.txt", sep='\t', index_col=False)
samples190 = pd.read_csv("samples190.txt", sep='\t', index_col=False)
sample_names = list(samples190['sample'])
sample_locations = list(samples190['location'])
samples_set = zip(sample_names, sample_locations)
samples_dict = dict(zip(sample_names, sample_locations))


samples183 = pd.read_csv("samples183.txt", sep='\t', index_col=False)
samples183_names = list(samples183['sample'])


# load sample info from reference clone files
refclones = pd.read_csv("daphnia70.list", sep='\t', index_col=False)
refclone_names = list(refclones['sample'])


# load sample info from duplicate removed bam files 
#samples253 = pd.read_csv("samples253.txt", sep='\t', index_col=False)
#sample253_names = list(samples253['sample'])
#sample253_locations = list(samples253['location'])
#samples253_set = zip(sample253_names, sample253_locations)
#samples253_dict = dict(zip(sample253_names, sample253_locations))


# load number of individuals for angsd_cutoffs
#nInd = pd.read_csv("list/cutoff_nInd.txt", sep='\t')
#IND = list(nInd['nInd'])

#prefixBAMdepth10 = pd.read_csv("list/prefix_realignedBAM_depth10.list", sep='\t')
#BAM = list(prefixBAMdepth10['prefix'])


###### helper functions ######

def getFqHome(sample):
  return(list(os.path.join(samples_dict[sample],"{0}_{1}_001.fastq.gz".format(sample,pair)) for pair in ['R1','R2']))

#def fastqcHome(sample):
#  return(list(os.path.join('raw/qc/fastqc',"{0}_001.fastqc.html")))

def getKrakenHome(sample):
    return(list(os.path.join('KRAKEN2_RESULTS', "{0}_{1}.trmdfilt.keep.fq.gz".format(sample,pair)) for pair in ['R1','R2']))

