import os
import pandas as pd

###### Config file and sample sheets #####


configfile: 
  "config/config.yaml"

# load sample info from raw data
#samples003 = pd.read_csv("samples003.txt", sep='\t', index_col=False)
#samples190 = pd.read_csv("samples190.txt", sep='\t', index_col=False)
#samples_info = pd.read_csv("list/LC_clones2021_hatchlings_metadata.list", sep='\t', index_col=False)
samples_info = pd.read_csv("list/da_eggs_hatchlings_pelagial_metadata.csv", sep='\t', index_col=False)
sample_names = list(samples_info['fastq_ID'])
sample_dir = list(samples_info['fastq_dir'])
#Kraken_dir = list(samples_info['Kraken_dir'])
sample_ID = list(samples_info['sample_id'])
samples_set = zip(sample_names, sample_dir)
samples_dict = dict(zip(sample_names, sample_dir))

# load sample info from reference clone files
samples_ALL = pd.read_csv("list/Da_Resteggs_metadata_cleaned_NEW.csv", sep='\t', index_col=False)
refclones = samples_ALL[samples_ALL['period']=='REF']
refclone_names = list(refclones['prefix'])
# load sample info from duplicate removed bam files 
ALL_names = list(samples_ALL['prefix'])

# load sample info from metadata
#samples253_new = pd.read_csv("samples253_metadata.txt", sep='\t', index_col=False)
#samples253_prefix = list(samples253_new['prefix'])
#samples253_id = list(samples253_new['sample_id'])

# load Chrom Info
#chromosom_information = pd.read_csv("list/dgal_rapid_ChromInfo.csv", sep=',', index_col=False)
#chromosom_information = pd.read_csv("list/dgal_rapid_Chrom_dgal1_dgal10.csv", sep=',', index_col=False)
### chromosom_information
#chrom = list(map(str, chromosom_information['chrom']))
#chromosom_length= list(chromosom_information['length'])
#start = list(chromosom_information['start'])
#end = list(chromosom_information['end'])


# load depth filters
if os.path.isfile("depth/stats/depthFilter.list"):
  print ("Depth filter file exists")
  depth_information = pd.read_csv("depth/stats/depthFilter.list", sep='\t')

  # number of samples (individuals)
  N = list(depth_information['Number_of_samples'])
  # minimum depth over all samples: 1 x number of samples
  MinDepth = 1*N
  # maximum over all samples: (mean depth + 3 x standard deviation) x number of samples
  #MaxDepth = int(depth_information['HengLi_max']*N)
  MaxDepth = [int(max) for max in depth_information['HengLi_max']*depth_information['Number_of_samples']]
else:
  print ("Depth file does not exist")

sets=['LC_REF', 'LC_withoutREF']

GL = ['2', '2']
minMaf = ['0.05', '0.05']

#filters = list(zip(sets, GL[1], minMaf[1], N, MinDepth, MaxDepth))

## angsd parameters
# minor allele frequency for PCA and admixture plots for 1% and 5% of the data
#minMaf = ['0.01', '0.05']
# minor allele frequency for SFS, keeping singletons
# SFS
#minMaf_SFS = round(2/(2*N),3)

# load minimum number of individuals a read has to be present
nInd = pd.read_csv("list/cutoff_nInd.txt", sep='\t')
# in 50 %, 75 % and 100 % of individuals
IND = list(nInd['nInd'])


# number of Ks for admixture proportions
#admix_K = ['2', '2', '3', '3', '4', '4', '5', '5', '6', '6', '7', '7', '8', '8', '9', '9', '10', '10']
admix_K = ['2', '3', '4', '5', '6']

POP = ['longispina_March21', 'longispina_June21']

## Population pairs for 2Dsfs and fst
pop_pair = pd.read_csv("list/pop_list/pop_pairs.list", sep='\t', index_col=False)
POP1 = list(pop_pair['POP1'])
POP2 = list(pop_pair['POP2'])


## 4-population combinations for ABBA-BABA (Dstats)
pop_combi = pd.read_csv("list/pop_list/4pop_combinations.list", sep=',', index_col=False)
pop_combi['combi'] = pop_combi[['P1', 'P2', 'P3', 'outgroup']].apply(lambda x: '_'.join(x), axis=1)
Dstats_combi = list(pop_combi['combi'])



## samples for heterozygosity estimates of fixed sites
bam_list = pd.read_csv("ancestry/parents_hybrids_BAM.list", sep='\t', index_col=False)
bams = list(bam_list['sample'])
#prefixBAMdepth10 = pd.read_csv("list/prefix_realignedBAM_depth10.list", sep='\t')
#BAM = list(prefixBAMdepth10['prefix'])


###### helper functions ######

def getFqHome(sample):
  return(list(os.path.join(samples_dict[sample],"{0}_{1}_001.fastq.gz".format(sample,pair)) for pair in ['R1','R2']))

def getKrakenHome(sample):
    return(list(os.path.join("/home/uibk/c7701178/lotte-c7701178-local/KRAKEN2_RESULTS/","{0}_{1}.trmdfilt.keep.fq.gz".format(sample,pair)) for pair in ['R1','R2']))

