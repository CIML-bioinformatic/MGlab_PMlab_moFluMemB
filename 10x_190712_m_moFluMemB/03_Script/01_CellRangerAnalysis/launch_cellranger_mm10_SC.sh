#!/bin/bash

# This script execute the CellRanger analysis of the 10x_190712_m_moFluMemB project data.
# To execute the script use singularity with the following command
# singularity exec /mnt/NAS7/Collaboration/moFluMemB/03_BIOINFO_ANALYSIS/10x_190712_m_moFluMemB/02_Container/cellranger3.0.1.img ./launch_cellranger_mm10_SC.sh

SAMPLE_NAME="10935372"
SINGLE_CELL_SUFFIX="SC"

FOLDER_TO_REFERENCE="/mnt/NAS7/PMlab/REF/V2/mm10"
#EXPECTED_CELLS=10000
OUTPUT_FOLDER=""/mnt/NAS7/Collaboration/moFluMemB/03_BIOINFO_ANALYSIS/10x_190712_m_moFluMemB/05_Output/01_CellRangerAnalysis/"${SAMPLE_NAME}"/mm10"
FASTQ_FOLDER="/mnt/NAS7/Collaboration/moFluMemB/01_BIOINFO_RAWDATA/10x_190712_m_moFluMemB/AHHY3LBGX7/outs/fastq_path/HHY3LBGX7/"${SAMPLE_NAME}${SINGLE_CELL_SUFFIX}

# Execution of the analysis
export PATH=/usr/local/share/cellranger/cellranger-3.0.1:$PATH
echo OUTPUT_FOLDER = $OUTPUT_FOLDER
echo FOLDER_TO_REFERENCE = $FOLDER_TO_REFERENCE

cellranger count --id=$OUTPUT_FOLDER --transcriptome=$FOLDER_TO_REFERENCE --fastqs=$FASTQ_FOLDER --sample=$SAMPLE_NAME


