#!/bin/bash

# This script execute the CellRanger analysis of the 10x_191105_m_moFluMemB project data.
# To execute the script use singularity with the following command
# singularity exec /mnt/NAS7/Collaboration/moFluMemB/03_BIOINFO_ANALYSIS/10x_191105_m_moFluMemB/02_Container/cellranger3.0.1.img ./launch_cellranger_mm10_SC.sh

SAMPLE_NAME="10938436"
SINGLE_CELL_SUFFIX="SC"

# Path to the folder with the reference genome
FOLDER_TO_REFERENCE="/mnt/NAS7/PMlab/REF/V2/mm10"

# Path to the folder with the fastq files
FASTQ_FOLDER="/mnt/NAS7/Collaboration/moFluMemB/01_BIOINFO_RAWDATA/10x_191105_m_moFluMemB/DIR191010_NB501865_0215_AHYGK7BGX9_Bis/AHYGK7BGX9/outs/fastq_path/HYGK7BGX9/"${SAMPLE_NAME}${SINGLE_CELL_SUFFIX}

# Path to the desired output file
OUTPUT_FOLDER="/mnt/NAS7/Collaboration/moFluMemB/03_BIOINFO_ANALYSIS/10x_191105_m_moFluMemB/05_Output/01_CellRangerAnalysis/"

# Execution of the analysis
export PATH=/usr/local/share/cellranger/cellranger-3.0.1:$PATH
echo OUTPUT_FOLDER = $OUTPUT_FOLDER
echo FOLDER_TO_REFERENCE = $FOLDER_TO_REFERENCE

mkdir -p $OUTPUT_FOLDER
cd $OUTPUT_FOLDER
cellranger count --id=${SAMPLE_NAME} --transcriptome=$FOLDER_TO_REFERENCE --fastqs=$FASTQ_FOLDER --sample=${SAMPLE_NAME}${SINGLE_CELL_SUFFIX}


