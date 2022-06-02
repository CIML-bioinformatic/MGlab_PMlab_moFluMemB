#!/bin/sh

# Prepocessing and launch CITE-seq-count
# Author: Lionel Spinelli, BIP, 12/11/2019

# The sample name as shown in the name of the fastq files
SAMPLE_NAME="10938436"

# The suffix of the fastq file name after the sample name
HTO_SUFFIX="HTO_S6"

# Location of the sample BCR fastq
SAMPLE_DIR="/mnt/NAS7/Collaboration/moFluMemB/01_BIOINFO_RAWDATA/10x_191105_m_moFluMemB/DIR191010_NB501865_0215_AHYGK7BGX9_Bis/AHYGK7BGX9/outs/fastq_path/HYGK7BGX9/"

# Location of the file with the HTO tag to sample name association
# This file contains two comma separated columns with no headers and no void lines
# Example:
#  TGCCATAGGTAC,M1LN
#  TGAAGAACGCCG,M1LG
HTO_TAG_LIST_FILE="//mnt/NAS7/Collaboration/moFluMemB/01_BIOINFO_RAWDATA/10x_191105_m_moFluMemB/meta-data/10x_051119_m_moFluMemB_taglist.csv"

# Location of the result of the cell ranger analysis on the RNA fastq
CELLRANGER_OUTPUT_DIR="/mnt/NAS7/Collaboration/moFluMemB/03_BIOINFO_ANALYSIS/10x_191105_m_moFluMemB/05_Output/01_CellRangerAnalysis/"${SAMPLE_NAME}

# Location of the singularity image folder
SINGULARITY_IMAGE_DIR="/mnt/NAS7/Collaboration/moFluMemB/03_BIOINFO_ANALYSIS/10x_191105_m_moFluMemB/02_Container/"

# Location of the desired CITEseq-count output folder
OUTPUT_DIR="/mnt/NAS7/Collaboration/moFluMemB/03_BIOINFO_ANALYSIS/10x_191105_m_moFluMemB/05_Output/02_CITEseqCountAnalysis/"${SAMPLE_NAME}


# ######################
# Concatenate
# ######################
echo "Concatenate fastq files..."

# Prepare intermediate folder for merging
OUTPUT_DIR_INTER=${OUTPUT_DIR}"/mergefastq"
mkdir -p $OUTPUT_DIR_INTER

#Concatenate R1
ls $SAMPLE_DIR/${SAMPLE_NAME}${HTO_SUFFIX}"_L00"*"_R1_001.fastq.gz"
cat $SAMPLE_DIR/${SAMPLE_NAME}${HTO_SUFFIX}"_L00"*"_R1_001.fastq.gz" >$OUTPUT_DIR_INTER/${SAMPLE_NAME}${HTO_SUFFIX}"_R1.fastq.gz"

#Concatenate R2
ls $SAMPLE_DIR/${SAMPLE_NAME}${HTO_SUFFIX}"_L00"*"_R2_001.fastq.gz"
cat $SAMPLE_DIR/${SAMPLE_NAME}${HTO_SUFFIX}"_L00"*"_R2_001.fastq.gz" >$OUTPUT_DIR_INTER/${SAMPLE_NAME}${HTO_SUFFIX}"_R2.fastq.gz"


# #######################################################
# Retrieve barcode file from cellranger RNA analysis
# and extract the number of lines in this file
# #######################################################
echo "Getting Cellranger count results..."

cp ${CELLRANGER_OUTPUT_DIR}"/outs/filtered_feature_bc_matrix/barcodes.tsv.gz" $OUTPUT_DIR
gunzip $OUTPUT_DIR"/barcodes.tsv.gz"
BARCODE_FILE=$OUTPUT_DIR"/barcodes.tsv"
CELL_NUMBER=$(wc -l $OUTPUT_DIR/barcodes.tsv)
CELL_NUMBER=$(echo $CELL_NUMBER| cut -d' ' -f 1)

# ##########################################
#     RUN THE CITEseqcount SINGULARITY IMAGE 
# ##########################################
echo "Run CITE-seq-count..."
cd $OUTPUT_DIR

singularity exec $SINGULARITY_IMAGE_DIR/citeseqcount.img CITE-seq-Count -R1 $OUTPUT_DIR_INTER/${SAMPLE_NAME}${HTO_SUFFIX}"_R1.fastq.gz" -R2 $OUTPUT_DIR_INTER/${SAMPLE_NAME}${HTO_SUFFIX}"_R2.fastq.gz" -t $HTO_TAG_LIST_FILE -cbf 1 -cbl 16 -umif 17 -umil 26 --max-errors 3 -cell $CELL_NUMBER -wl $BARCODE_FILE -o $OUTPUT_DIR