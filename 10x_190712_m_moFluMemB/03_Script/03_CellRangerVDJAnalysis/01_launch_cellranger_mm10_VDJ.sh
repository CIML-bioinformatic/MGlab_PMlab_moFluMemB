#!/bin/bash

# This script execute the CellRanger analysis of the 10x_190712_m_moFluMemB project data.
# To execute the script use singularity with the following command
# singularity exec /mnt/NAS7/Collaboration/moFluMemB/03_BIOINFO_ANALYSIS/10x_190712_m_moFluMemB/02_Container/cellranger3.0.1.img /mnt/NAS7/Collaboration/moFluMemB/03_BIOINFO_ANALYSIS/10x_190712_m_moFluMemB/03_Script/03_CellRangerVDJAnalysis/01_launch_cellranger_mm10_VDJ.sh 

SAMPLE_NAME="10935372"
LIBRARY_SUFFIX="BCR"

FOLDER_TO_REFERENCE="/mnt/NAS7//Collaboration/moFluMemB/02_BIOINFO_REFERENCE/10xGenomics/VDJ/refdata-cellranger-vdj-GRCm38-alts-ensembl-4.0.0"
OUTPUT_FOLDER="/mnt/NAS7/Collaboration/moFluMemB/03_BIOINFO_ANALYSIS/10x_190712_m_moFluMemB/05_Output/03_CellRangerAnalysisVDJ/"${SAMPLE_NAME}"/mm10"
FASTQ_FOLDER="/mnt/NAS7/Collaboration/moFluMemB/01_BIOINFO_RAWDATA/10x_190712_m_moFluMemB/AHHY3LBGX7/outs/fastq_path/HHY3LBGX7/"${SAMPLE_NAME}${LIBRARY_SUFFIX}

# Execution of the analysis
export PATH=/usr/local/share/cellranger/cellranger-3.0.1:$PATH
echo FASTQ_FOLDER= $FASTQ_FOLDER
echo OUTPUT_FOLDER = $OUTPUT_FOLDER
echo FOLDER_TO_REFERENCE = $FOLDER_TO_REFERENCE

cellranger vdj --id=$SAMPLE_NAME --fastqs=$FASTQ_FOLDER --reference=$FOLDER_TO_REFERENCE --sample=${SAMPLE_NAME}${LIBRARY_SUFFIX}

