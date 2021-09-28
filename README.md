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

As described in the article, there is 3 datasets in this study. 

* 10x_190712_m_moFluMemB : 10x 5’ scRNA-Seq on single-cell suspensions from spleen, lymph nodes and lungs with enzymatic digestion of lung tissue at 37°C.
* 10x_191105_m_moFluMemB :  : 10x 5’ scRNA-Seq on single-cell suspensions from spleen, lymph nodes and lungs with mechanical dissociation of lung tissue at 4°C.
* custom_201216_m_moFluMemB : FB5P-seq protocol (Attaf et al., 2020) on single-cell suspensions from lungs with enzymatic digestion, and stained with a panel of antibodies for identifying subsets of 

When downloading the code and data (see below), you will obtains 3 sub-folders with names as below:

```
    moFluMemB
    ├── 10x_190712_m_moFluMemB
    ├── 10x_191105_m_moFluMemB
    └── custom_201216_m_moFluMemB
```

---
---

## Prepare the environments

In order to prepare the environment for analysis execution, it is required to:

- Clone the github repository and set the WORKING_DIR environment variable
- Download the pre-processed data
- Download the docker image tar file and the singularity img files
- Install Docker and Singularity
- Load the docker image on your system

Below you will find detailed instruction for each of these steps.

### Clone the github repository

Use you favorite method to clone this repository in a chosen folder. This will create a folder **moFluMemB** with all the source code. 

Then, you must set an environment variable called **WORKING_DIR** with a value set to the path to this folder.

For instance, if you have chosen to clone the Git repository in __"/home/spinellil/workspace"__, then the **WORKING_DIR** variable will be set to __"/home/spinellil/workspace/moFluMemB"__

**On linux:**

```
    export WORKING_DIR=/home/spinellil/workspace/moFluMemB
```

### Download the data

Each sample needs its own sub-folder containing the initial data used by the analysis. Those data can be downloaded from Zenodo and uncompressed. The Zenodo dataset DOI are **TODO: Add the 3 datasets DOI**. The initial data from the analysis are the pre-processed data :

* CellRangerAnalysis (for the two first datasets) : contains the result of Cell Ranger count analysis from the mRNA fastq
* CITESeqCountAnalysis (for the two first datasets) : contains the result of CITE-seq-Count analysis from the HTP fastq
* BCRAnalysis (for all datasets): contains the BCR repertoire analysis and for the FB5P-seq dataset also the mRNA count tables and the Protein count tables.

Scripts of the pre-processing steps is provided for Cell Ranger count analysis and CITE-seq-Count analysis and corresponding raw data (fastq files) can be downloaded from GEO (see article). For FB5P-seq pre-processing see pipeline description in (Attaf et al., 2020).

To download and uncompress the data, use the following code:

**On linux:**

```
    cd $WORKING_DIR
    wget **TODO: Add the dataset1 URL** -O 10x_190712_m_moFluMemB_processedData.tar.gz
    tar zxvf 10x_190712_m_moFluMemB_processedData.tar.gz
    
    wget **TODO: Add the dataset2 URL** -O 10x_191105_m_moFluMemB_processedData.tar.gz
    tar zxvf 10x_191105_m_moFluMemB_processedData.tar.gz

    wget **TODO: Add the dataset3 URL** -O custom_201216_m_moFluMemB_processedData.tar.gz
    tar zxvf custom_201216_m_moFluMemB_processedData.tar.gz
```
 
Once done, you may obtain the following subfolder structure, each of them containing several files.

```
    moFluMemB
    ├── 10x_190712_m_moFluMemB
    │   ├── 05_Output/01_CellRangerAnalysis
    │   ├── 05_Output/02_CITEseqCountAnalysis
    │   └── 05_Output/03_BCRAnalysis
    ├── 10x_191105_m_moFluMemB
    │   ├── 05_Output/01_CellRangerAnalysis
    │   ├── 05_Output/02_CITEseqCountAnalysis
    │   └── 05_Output/03_BCRAnalysis
    └── custom_201216_m_moFluMemB
        └── 05_Output/03_BCRAnalysis

```

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

The analysis uses the Singularity images and optionnaly Snakemake.

The study contains 3 samples of single-cell RNA-seq data. Each sample have several steps of analysis you will find the R script files in the subfolder **03_Script**.

Each step of analysis generates its own HTML report file and several output files. Some output files of some steps are used by other steps, making a complete workflow of analysis.

The simpliest way to run the complete single-cell analysis of a sample is to use the Snakemake workflow dedicated to each sample. The workflow is controled by a snakefile stored in the **04_Workflow** subfolder of each sample folder. This workflow uses Singularity images (see above) to control the software environment for each analysis step. So you need both Snakemake and Singularity installed on your system to use this workflow.

In order to use the snakemake workflow, please type first the following commands:

     cd $WORKING_DIR
     ln -s Embryo_Stage13.5_FetalLiver/04_Workflow/snakefile.yml Embryo_Stage13.5_FetalLiver/snakefile.yml
     ln -s Embryo_Stage13.5_Periphery_CellRangerV3/04_Workflow/snakefile.yml Embryo_Stage13.5_Periphery_CellRangerV3/snakefile.yml
     ln -s Embryo_Stage14.5_FetalLiver/04_Workflow/snakefile.yml Embryo_Stage14.5_FetalLiver/snakefile.yml
     ln -s Embryo_Stage14.5_Periphery_CellRangerV3/04_Workflow/snakefile.yml Embryo_Stage14.5_Periphery_CellRangerV3/snakefile.yml

To run the analysis for the Embryo_Stage13.5_FetalLiver (for instance), then run the following commands:

Note: you have to manually change the **$WORKING_DIR** string in the snakemake command below by the value of the environment variable (i.e the path where you clone the project) because snakemake may not interpret the variable name correctly:

     cd $WORKING_DIR/Embryo_Stage13.5_FetalLiver
     snakemake -r --snakefile snakefile.yml --use-singularity --singularity-args "-B $WORKING_DIR:$WORKING_DIR"
     
To execute the analysis of the other sample, simply change folder to the target sample and run again the same snakemake command.



