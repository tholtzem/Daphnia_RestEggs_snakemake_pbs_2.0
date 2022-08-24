#! /bin/bash


SAMPLELIST=list/realignedBAM_df3_prefixes.list
 # Path to a list of prefixes of the raw fastq files. It should be a subset of the the 1st column of the fastq table.

SAMPLETABLE=list/samples184_metadata.tsv # Path to a fastq table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID. 


#for SAMPLEFILE in `cat $SAMPLELIST`; do
#	# For each prefix, extract the associated species (column 2) from the table
#	SPECIES=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 3`
#	#echo $SAMPLEFILE
#	echo $SPECIES
#        #For each species, make a list of bam files
#        #ls -v realigned/${SAMPLEFILE}*bam >> list/pop_list/${SPECIES}.list
#done


#for SAMPLEFILE in `cat list/pop_list/eggs.list | cut -f2 -d '/' | cut -f1 -d'.'`; do
#	# For each prefix, extract the associated species (column 2) from the table
#	LAKE=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 6`
#	echo $SAMPLEFILE
#	echo $LAKE
#        # For each lake, make a list of bam files
#        ls -v realigned/$SAMPLEFILE*bam >> list/pop_list/eggs_${LAKE}.list
#done


for SAMPLEFILE in `cat list/pop_list/eggs.list | cut -f2 -d '/' | cut -f1 -d'.'`; do
	# For each prefix, extract the associated species (column 2) from the table
	LAKE=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 6`
	PERIOD=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 15`
	echo $SAMPLEFILE
	#echo $LAKE
	#echo $PERIOD
        # For each period, make a list of bam files
        #ls -v realigned/$SAMPLEFILE*bam >> list/pop_list/eggs_${LAKE}_${PERIOD}.list
done


