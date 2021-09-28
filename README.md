# Permissive clonal selection engenders discrete subsets of memory B cells with divergent transcriptional programs and specificities to viral antigens

## Article information

**Title:** Permissive clonal selection engenders discrete subsets of memory B cells with divergent transcriptional programs and specificities to viral antigens.

**Authors:** 
1 Aix-Marseille Univ, Centre National de la Recherche Scientifique (CNRS), Institut National de la Santé et de la Recherche Médicale (INSERM), Centre d'Immunologie de Marseille-Luminy (CIML), Marseille, France.

% Corresponding author: E-mail: 

**Summary:**
Memory B cells are key cellular components of the long-term humoral immunity that dominate recall responses by rapidly differentiating into effector cells. As the biology of memory B cells has been mainly studied in the context of synthetic antigen immunizations, the dynamics of these cells during infection remains largely unexplored. Here, we used respiratory infection models with influenza and SARS-CoV-2, and fluorescent-reporter mice to identify memory B cells regardless of antigen-specificity. Confocal imaging and scRNA-seq data revealed that three transcriptionally distinct subsets of memory B cells with divergent origin and tissue localization coexist in the lung mucosa upon resolution of infection. This comprises a minor subset of “innate-like” IgM+ cells and two large subsets of class-switched somatically mutated memory B cells that preferentially differentiate into plasma rather than germinal center cells upon activation. Concomitant analysis of antigen-specificity and B cell receptor repertoire revealed that these two larger subsets segregate “bonafide” virus-specific cells from “bystander” memory B cells with no apparent specificity for viral antigens, originated as product of permissive clonal selection. Thus, diverse transcriptional programs in memory B cells are not associated with specific effector fates but rather with divergent strategies of the immune system to simultaneously provide rapid protection from re-infection while diversifying the initial B cell repertoire.

**Keywords:**
Memory B cells - Respiratory infection - Influenza virus - SARS-CoV-2 - Lung mucosa - Permissive selection

DOI : 

---
---

## Goal of the github
This github project contains the instructions and material to reproduce the analysis reported in the article (and more).
Source code (scripts and dockerfiles) are available in the github repository. Required data and builded Docker/Singularity images are available on download. Instructions to reproduce the analysis are provided below.

To reproduce the analysis, you have to first, prepare the environments (see "Prepare the Environments" section below), then execute the analysis step by step (see "Run the analysis" section below).

---
---

## Description of the datasets

As described in the article, there is 3 datasets in this study. **TODO:describe datasets** . When downloading the code and data, you will obtains 3 sub-folders with names as below:

    moFluMemB
    ├── 10x_190712_m_moFluMemB : Single-cell RNA-seq of **TODO:describe datasets**
    ├── 10x_191105_m_moFluMemB :  : Single-cell RNA-seq of **TODO:describe datasets**
    └── custom_201216_m_moFluMemB : Single-cell RNA-seq of **TODO:describe datasets**

---
---

## Prepare the environments

In order to prepare the environment for analysis execution, it is required to:

- Clone the github repository and set the WORKING_DIR environment variable
- Download the docker image tar file and the singularity img files
- Install Docker and Singularity
- Load the docker image on your system
- Download the pre-processed data (Count table for bulk RNA-seq and CellRanger results for single-cell RNA-seq)

Below you will find detailed instruction for each of these steps.

### Clone the github repository

Use you favorite method to clone this repository in a chosen folder. This will create a folder **"moFluMemB"** with all the source code. 

Then, you must set an environment variable called **WORKING_DIR** with a value set to the path to this folder.

For instance, if you have chosen to clone the Git repository in __"/home/spinellil/workspace"__, then the WORKING_DIR variable will be set to __"/home/spinellil/workspace/moFluMemB"__

**On linux:**

    export WORKING_DIR=/home/spinellil/workspace/moFluMemB

### Download the raw data

Each sample needs its own **"00_RawData"** sub-folder containing the initial data used by the analysis. Those data can be downloaded from Zenodo and uncompressed. The Zenodo dataset DOI are **TODO: Add the 3 datasets DOI**.

To download and uncompress the data, use the following code:

**On linux:**

    cd $WORKING_DIR
    wget **TODO: Add the dataset1 URL** -O **TODO: Add the dataset1 file name**
    tar zxvf **TODO: Add the dataset1 file name**
    
    wget **TODO: Add the dataset2 URL** -O **TODO: Add the dataset1 file name**
    tar zxvf **TODO: Add the dataset2 file name**

    wget **TODO: Add the dataset3 URL** -O **TODO: Add the dataset1 file name**
    tar zxvf **TODO: Add the dataset3 file name**
    
Once done, you may obtain the following subfolder structure, each of them containing several files.

    moFluMemB
    ├── 10x_190712_m_moFluMemB
    │   └── 00_RawData
    ├── 10x_191105_m_moFluMemB
    │   └── 00_RawData
    └── custom_201216_m_moFluMemB
        └── 00_RawData

### Download the reference files

The study uses references (genome annotations) you have to download. The annotations used during the study are available on Zenodo **TODO: Add the reference DOI**. Use the following command to download the tarball file and uncompress it.

Note: Since the reference files are used for the 3 single-cell samples analysis, they must be present in all the sample folder in the same **"01_Reference"** subfolder. Instead of copying the files, we will create symbolic links:

**On linux:**

    cd $WORKING_DIR
    wget **TODO: Add the reference URL** -O **TODO: Add the reference file name**
    tar zxvf **TODO: Add the reference file name**
    ln -s 10x_190712_m_moFluMemB/01_Reference 10x_191105_m_moFluMemB/01_Reference
    ln -s 10x_190712_m_moFluMemB/01_Reference custom_201216_m_moFluMemB/01_Reference

These commands will create 4 sub-folders named 01_Reference:

    moFluMemB
    ├── 10x_190712_m_moFluMemB
    │   └── 01_Reference
    ├── 10x_191105_m_moFluMemB
    │   └── 01_Reference
    └── custom_201216_m_moFluMemB
        └── 01_Reference

### Download the Docker and Singularity images

Docker image tar file and Singularity img files are stored on Zenodo  **TODO: Add the containers DOI**. Open a shell command and change dir to the root of the cloned Git repository (WORKING_DIR). Then execute the following commands to download the tarball file and untar  it:

**On linux:**

    cd $WORKING_DIR
    wget **TODO: Add the container URL** -O **TODO: Add the container file name**
    tar zxvf **TODO: Add the container file name**

These commands will create 3 sub-folders named 02_Container:

    moFluMemB
    └── custom_201216_m_moFluMemB
        └── 02_Container

The first one contains a Docker image tar file used for the bulk RNA-seq analysis. The second one contains the Singularity images for the single-cell RNA-seq analysis. Since the singularity images are used for the 4 single-cell samples analysis, they must be present in all the sample folder in the same 02_Container subfolder. Instead of copying the image files, we will create symbolic links:

**On linux:**

    cd $WORKING_DIR
    ln -s 10x_190712_m_moFluMemB/02_Container 10x_191105_m_moFluMemB/02_Container
    ln -s 10x_190712_m_moFluMemB/02_Container custom_201216_m_moFluMemB/02_Container

### Install Docker and Singularity

You need to install Docker and Singularity v2.6 on your system.

- To install Docker, follow the instructions here : https://docs.docker.com/get-docker/

- To install Singularity v2.6, follow the instructions here : https://sylabs.io/guides/2.6/admin-guide/

### Install Snakemake

If you want to take advantage of the workflow management we used for the single-cell RNA-seq analysis, you have to install Snakemake. See the official instruction and use your prefered solution:

https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

---
---

## Run the analysis

There are two types of analysis in this study : bulk RNA-seq and single-cell RNA-seq. The bulk RNA-seq analysis uses the Docker image you loaded. The single-cell RNA-seq analysis uses the Singularity images and optionnaly Snakemake.

### Run the bulk RNA-seq analysis

The RNA-seq analysis are in two steps (step1 and step2). The first step make the QC, study the differentially expressed genes and their functionnal enrichment. The second step study the pattern of evolution of group of genes along the cell types (see article methods).

To run the step1 analysis, use the following command:

**On Linux:**

    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR splab_ilcyou_deg_gsea 'cd $WORKING_DIR/Embryo_Bulk_Stage13.5_2tissues/03_Script/step1;Rscript launch_reports_compilation.R'

To run the step2 analysis, use the following command:

**On Linux:**

     docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR splab_ilcyou_deg_gsea 'cd $WORKING_DIR/Embryo_Bulk_Stage13.5_2tissues/03_Script/step2;Rscript launch_reports_compilation.R'

Each analysis will generate a result in $WORKING_DIR/Embryo_Bulk_Stage13.5_2tissues/05_output/step1 or $WORKING_DIR/Embryo_Bulk_Stage13.5_2tissues/05_output/step2.
In the output of the analysis, you will find a HTML file that contains the report of the analysis, with all figures. Some extra file are generated to export data in plain text.


### Run the single-cell RNA-seq analysis

The study contains 4 samples of single-cell RNA-seq data. Each sample have 5 step of analysis you will find the R script files in the subfolder 03_Script. The 5 steps are:

 * 01_QC : General quality control and bad cell removal
 * 02_GlobalHeterogeneity : First study of cell heterogeneity and sample contamination by undesired cell types
 * 03_GlobalHeterogeneity_NoContamination : Study of cell heterogeniety in absence of contamination
 * 04_Dynamics_Monocle : analysis of the cellular process dynamics using pseudotime analysis by Monocle
 * 05_Dynamics_RNAVelocity : analysis of the cellular process dynamics using RNA velocity (Velocyto)

Each step of analysis generates its own HTML report file and several output files. Some output files of some steps are used by other steps, making a complete workflow of analysis.

The simpliest way to run the complete single-cell analysis of a sample is to use the Snakemake workflow dedicated to each sample. The workflow is controled by a snakefile stored in the 04_Workflow subfolder of each sample folder. This workflow uses Singularity images (see above) to control the software environment for each analysis step. So you need both Snakemake and Singularity installed on your system to use this workflow.

In order to use the snakemake workflow, please type first the following commands:

     cd $WORKING_DIR
     ln -s Embryo_Stage13.5_FetalLiver/04_Workflow/snakefile.yml Embryo_Stage13.5_FetalLiver/snakefile.yml
     ln -s Embryo_Stage13.5_Periphery_CellRangerV3/04_Workflow/snakefile.yml Embryo_Stage13.5_Periphery_CellRangerV3/snakefile.yml
     ln -s Embryo_Stage14.5_FetalLiver/04_Workflow/snakefile.yml Embryo_Stage14.5_FetalLiver/snakefile.yml
     ln -s Embryo_Stage14.5_Periphery_CellRangerV3/04_Workflow/snakefile.yml Embryo_Stage14.5_Periphery_CellRangerV3/snakefile.yml

To run the analysis for the Embryo_Stage13.5_FetalLiver (for instance), then run the following commands:

Note: you have to manually change the "$WORKING_DIR" string in the snakemake command below by the value of the environment variable (i.e the path where you clone the project) because snakemake may not interpret the variable name correctly:

     cd $WORKING_DIR/Embryo_Stage13.5_FetalLiver
     snakemake -r --snakefile snakefile.yml --use-singularity --singularity-args "-B $WORKING_DIR:$WORKING_DIR"
     
To execute the analysis of the other sample, simply change folder to the target sample and run again the same snakemake command.


