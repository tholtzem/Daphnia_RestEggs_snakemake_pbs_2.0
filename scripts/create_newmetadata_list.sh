#! /bin/bash



# Path to the snakemake base directory
BASEDIR=(/home/uibk/c7701178/scratch/DAPHNIA/Daphnia_RestEggs_snakemake_pbs_2.0) 

# Path to a list of realigned bam files that passed that passed a certain mean depth (here: > 3)
BAMS=$BASEDIR/depth/stats/realignedBAM_df1.list

# Path to a list of prefixes of the realigned bam files. It should be a subset of the the 1st column of the table containing the metadata.
SAMPLELIST=$BASEDIR/list/realignedBAM_df3_prefixes.list

# Path to table containing the metadata (prefix, lane_number, sample_id, data_type, species, location, region,...)
#SAMPLETABLE=$BASEDIR/samples253_metadata.txt
SAMPLETABLE=$BASEDIR/list/Da_Resteggs_metadata_cleaned_NEW.csv

# Path to temporary table of metadata (without header)
TMP=$BASEDIR/list/tmp.tsv

NEW_SAMPLETABLE=$BASEDIR/list/angsd_metadata.tsv

# Path to the header of the table
HEADER=$BASEDIR/list/header_metadata.tsv


# First, create a list of prefixes of the realigned bam files
cat $BAMS | cut -f2 -d'/' | sed 's/.realigned.bam//g' > $SAMPLELIST

# Second, count the number of realigned bam files
N=$(cat $BAMS | wc -l)
echo $N

# Third, print the header of the table containing the metadata
cat $SAMPLETABLE | head -n1 > $HEADER 

# Finally, create a new metadata file for the samples that will be further processed
for SAMPLEFILE in `cat $SAMPLELIST`; do   # Loop through each of the prefixes
	#For each prefix, extract the associated sample ID (column 4) from the table
	SAMPLE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE`
	#echo $SAMPLEFILE refers to sample $SAMPLE_ID
	echo -e "$SAMPLE_ID"
done > $TMP


cat $HEADER $TMP > $NEW_SAMPLETABLE
