#!/bin/sh

# Prepocessing and lanch CITE-seq
# Author: DONG Chuang CIML/PMlab 2019/07/10

#if [ $# -lt 1 ] ; then
	echo  "\nAuthor: DONG Chuang CIML/PMlab 05.31.2019"
	echo  "USAGE: $0 <command>"
	echo  "Use for Pipeline cite-seq"
	echo  "Prepocessing and lanch CITE-seq"
	echo  "For Example	: $0 toto\n"
#fi

SAMPLE_NAME="10935372"
HTO_SUFFIX="HTO_S5"

sampledir="/mnt/NAS7/Collaboration/moFluMemB/01_BIOINFO_RAWDATA/10x_190712_m_moFluMemB/AHHY3LBGX7/outs/fastq_path/HHY3LBGX7/"
outputdir="/mnt/NAS7/Collaboration/moFluMemB/03_BIOINFO_ANALYSIS/10x_190712_m_moFluMemB/05_Output/02_CITEseqCountAnalysis/"${SAMPLE_NAME}
outputdir_inter=${outputdir}"/mergefastq"

singularity_imagedir="/mnt/NAS7/Collaboration/moFluMemB/03_BIOINFO_ANALYSIS/10x_190712_m_moFluMemB/02_Container/"

mkdir -p $outputdir_inter

# ######################
# Concatenate
# ######################
echo "Concatenate... ..."
#Concatenate R1
ls $sampledir/${SAMPLE_NAME}${HTO_SUFFIX}"_L00"*"_R1_001.fastq.gz"
cat $sampledir/${SAMPLE_NAME}${HTO_SUFFIX}"_L00"*"_R1_001.fastq.gz" >$outputdir_inter/${SAMPLE_NAME}${HTO_SUFFIX}"_R1.fastq.gz"

#Concatenate R2
ls $sampledir/${SAMPLE_NAME}${HTO_SUFFIX}"_L00"*"_R2_001.fastq.gz"
cat $sampledir/${SAMPLE_NAME}${HTO_SUFFIX}"_L00"*"_R2_001.fastq.gz" >$outputdir_inter/${SAMPLE_NAME}${HTO_SUFFIX}"_S5_R2.fastq.gz"


# ##################################################################
# ##################################################################
# ##################################################################
# DANGER : TRIMMING STEP BEFORE CITE SEQ:
# ##################################################################
# ##################################################################
# ##################################################################

# NEED TO DO A TRIMMING STEP BEFORE CITE-SEQ
# #######################

#echo "trimming 1 ... ..."
#singularity run docker://comics/trimmomatic java -jar /software/applications/Trimmomatic/0.36/trimmomatic-0.36.jar PE $sampledir/$sampleID"_R1.fastq.gz" $sampledir/$sampleID"_R2.fastq.gz" $sampledir/$sampleID"_R1_TRIM_PAIRED.fastq.gz" $sampledir/#$sampleID"_R1_TRIM_UNPAIRED.fastq.gz" $sampledir/$sampleID"_R2_TRIM_PAIRED.fastq.gz" $sampledir/$sampleID"_R2_TRIM_UNPAIRED.fastq.gz" CROP:50


#echo "trimming 2 ... ..."
# AND NEED TO DO A REMOVE THE READ UNDER 26bp (for read1) and 50bp (for read2)
# Using an other docker
# #######################

#singularity run docker://dceoy/cutadapt --minimum-length 26:50 -o $sampledir/$sampleID"_R1_TRIM_PAIRED_MINCUT.fastq.gz" -p $sampledir/$sampleID"_R2_TRIM_PAIRED_MINCUT.fastq.gz" $sampledir/$sampleID"_R1_TRIM_PAIRED.fastq.gz" $sampledir/$sampleID"_R2_TRIM_PAIRED.fastq.gz"
#exit 1
# ##################################################################
# ##################################################################
# ##################################################################
# END of trimming
# ##################################################################
# ##################################################################
# ##################################################################

echo "cite-seq... ..."
# ##########################################
#     RUN THE CITEseqcount SINGULARITY IMAGE 
# ##########################################
cd $outputdir
 
singularity exec $singularity_imagedir/citeseqcount.img CITE-seq-Count -R1 $outputdir_inter/${SAMPLE_NAME}${HTO_SUFFIX}"_R1.fastq.gz" -R2 $outputdir_inter/${SAMPLE_NAME}${HTO_SUFFIX}"_R2.fastq.gz" -t TAG_LIST.csv -cbf 1 -cbl 16 -umif 17 -umil 26 --max-errors 3 -cell 5515 -wl barcodes.tsv -o $outputdir
