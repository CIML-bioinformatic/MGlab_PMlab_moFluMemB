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

**DOI:**

---
---

## Goal of the github
This github project contains the instructions and material to reproduce the analyses reported in the article (and more).
Source code (scripts and dockerfiles) are available in the github repository. Required data and built Docker/Singularity images are available on download. Instructions to reproduce the analyses are provided below.

To reproduce the analysis, you have to first, prepare the environments (see "Prepare the Environments" section below), then execute the analysis step by step (see "Run the analysis" section below).

## Description of the datasets

As described in the article, there are 3 datasets in this study. 

* 10x_190712_m_moFluMemB : 10x 5’ scRNA-Seq on memory B cells sorted from single-cell suspensions of spleen, lymph nodes and lungs with enzymatic digestion of lung tissue at 37°C.
* 10x_191105_m_moFluMemB :  : 10x 5’ scRNA-Seq on memory B cells sorted from single-cell suspensions of spleen, lymph nodes and lungs with mechanical dissociation of lung tissue at 4°C.
* custom_201216_m_moFluMemB : FB5P-seq protocol (Attaf et al., 2020) on memory B cells sorted from single-cell suspensions of lungs with enzymatic digestion of lung tissue at 37°C, with index sorting information for a panel of antibodies identifying subsets of memory B cells.

When downloading the code and data (see below), you will obtains 3 sub-folders with names as below:

```
    moFluMemB
    │
    ├── 10x_190712_m_moFluMemB
    │
    ├── 10x_191105_m_moFluMemB
    │
    └── custom_201216_m_moFluMemB
```

---
---

## Prepare the environments

In order to prepare the environment for analysis execution, it is required to:

- Clone the github repository and set the WORKING_DIR environment variable
- Download the pre-processed data
- Install Singularity and Docker
- Install Snakemake
- Download the Singularity images
- (Optional) Download the Docker images and load them on your system

Below you will find detailed instruction for each of these steps.

---

### Clone the github repository

Use you favorite method to clone this repository in a chosen folder. This will create a folder **moFluMemB** with all the source code. 

---

### Set the WORKING_DIR variable

Then, you must set an environment variable called **WORKING_DIR** with a value set to the path to this folder.

For instance, if you have chosen to clone the Git repository in __"/home/spinellil/workspace"__, then the **WORKING_DIR** variable will be set to __"/home/spinellil/workspace/moFluMemB"__

**On linux:**

```
    export WORKING_DIR=/home/spinellil/workspace/moFluMemB
```

---

### Add you working dir in the code

The code uses variables that are stored in different "parameters" file. One important variable is the PATH_PROJECT which indicate to the code where your project is stored.
You have to modify this variable in the code to reflect your project setup. Each dataset has a file called **globalParams.R** in the subfolder **03_Script**

```
    moFluMemB
    │
    ├── 10x_190712_m_moFluMemB
    │   │
    │   └── 03_Script/globalParams.R
    │
    ├── 10x_191105_m_moFluMemB
    │   │
    │   └── 03_Script/globalParams.R
    │
    └── custom_201216_m_moFluMemB
        │
        └── 03_Script/globalParams.R
```

Edit those files and in each of them locate the line defining the **PATH_PROJECT** variable and change its value to the same value as the **WORKING_DIR** variable you defined before. Then save the files.

```
PATH_PROJECT = "/home/spinellil/workspace/moFluMemB"
```

---

### Download the data

Each sample needs its own sub-folder containing the initial data used by the analysis. Those data can be downloaded from Zenodo and uncompressed. The Zenodo dataset DOI are [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5564625.svg)](https://doi.org/10.5281/zenodo.5564625), [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5565864.svg)](https://doi.org/10.5281/zenodo.5565864) and [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5566675.svg)](https://doi.org/10.5281/zenodo.5566675). The initial data from the analysis are the pre-processed data :

* CellRangerAnalysis (for the two first datasets) : contains the result of Cell Ranger count analysis from the mRNA fastq
* CITESeqCountAnalysis (for the two first datasets) : contains the result of CITE-seq-Count analysis from the HTO fastq
* BCRAnalysis (for all datasets): contains the BCR repertoire analysis and for the FB5P-seq dataset also the mRNA count tables and the Protein count tables (index sorting).

Scripts of the pre-processing steps are provided for Cell Ranger count analysis and CITE-seq-Count analysis and corresponding raw data (fastq files) can be downloaded from GEO (see article). For FB5P-seq pre-processing see pipeline description in (Attaf et al., 2020).

To download and uncompress the data, use the following code:

**On linux:**

```
    cd $WORKING_DIR
    wget https://zenodo.org/record/5564625/files/10x_190712_m_moFluMemB_processedData.tar.gz -O 10x_190712_m_moFluMemB_processedData.tar.gz
    tar zxvf 10x_190712_m_moFluMemB_processedData.tar.gz
    
    wget https://zenodo.org/record/5565864/files/10x_191105_m_moFluMemB_processedData.tar.gz -O 10x_191105_m_moFluMemB_processedData.tar.gz
    tar zxvf 10x_191105_m_moFluMemB_processedData.tar.gz

    wget https://zenodo.org/record/5566675/files/custom_201216_m_moFluMemB_processedData.tar.gz -O custom_201216_m_moFluMemB_processedData.tar.gz
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

---

### Install Singularity and Docker

You need to install Singularity v2.6 on your system to run the complete analysis. Follow the instructions here : https://sylabs.io/guides/2.6/admin-guide/

Optionally, you can install Docker on your system to take advantage of interactive analysis environment with Rstudio, follow the instructions here : https://docs.docker.com/get-docker/

---

### Install Snakemake

You need to install Snakemake to run the complete analysis workflow. Use your preferred solution : https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

---

### Download the Singularity images

Singularity images files are stored on Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5566675.svg)](https://doi.org/10.5281/zenodo.5566675). Open a shell command and change dir to the root of the cloned Git repository (WORKING_DIR). Then execute the following commands to download the tarball file and untar  it:

**On linux:**

```
    cd $WORKING_DIR
    wget https://zenodo.org/record/5566675/files/moFluMemB_SingularityImages.tar.gz -O moFluMemB_SingularityImages.tar.gz
    tar zxvf moFluMemB_SingularityImages.tar.gz
```

These commands will create a sub-folder named **02_Container** in the first dataset folder:

```
    moFluMemB
    └── 10x_190712_m_moFluMemB
        └── 02_Container
```

This folder contains the Singularity images for the single-cell RNA-seq analysis. Since the Singularity images are used for the 3 single-cell samples analyses, they must be present in all the sample folder in the same **02_Container** subfolder. Instead of copying the image files, we will create symbolic links to spare disk space:

**On linux:**

```
    cd $WORKING_DIR
    ln -s $WORKING_DIR/10x_190712_m_moFluMemB/02_Container/mglab_moflumemb_seuratrf.img $WORKING_DIR/10x_191105_m_moFluMemB/02_Container/mglab_moflumemb_seuratrf.img
    ln -s $WORKING_DIR/10x_190712_m_moFluMemB/02_Container/mglab_moflumemb_seuratrf.img $WORKING_DIR/custom_201216_m_moFluMemB/02_Container/mglab_moflumemb_seuratrf.img
```

### (Optional) Download the Docker images and load them on your system

Docker image tar files are stored on Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5566675.svg)](https://doi.org/10.5281/zenodo.5566675). Open a shell command and change dir to the root of the cloned Git repository (WORKING_DIR). Then execute the following commands to download the tarball file, untar it and load the docker images on your system: 

```
    cd $WORKING_DIR
    wget https://zenodo.org/record/5566675/files/moFluMemB_DockerImages.tar.gz -O moFluMemB_DockerImages.tar.gz
    tar zxvf moFluMemB_DockerImages.tar.gz
    docker load -i mglab_moflumemb_seurat_r3_6.tar
```

---
---

## Run the analysis

### Run the analysis workflow using Snakemake

The analysis workflow uses the Singularity images and Snakemake.

The study contains 3 samples of single-cell RNA-seq data. Each sample have several steps of analysis for which you will find the R script files in the subfolder **03_Script**.

Each step of analysis generates its own HTML report file and several output files. Some output files of some steps are used by other steps, making a complete workflow of analysis. The output files are stored in each datatset folder in a sub-folder named "05_Output".

The simpliest way to run the complete single-cell analysis of a sample is to use the Snakemake workflow dedicated to each sample. The workflow is controlled by a snakefile stored in the **04_Workflow** subfolder of each sample folder. This workflow uses Singularity images (see above) to control the software environment for each analysis step. So you need both Snakemake and Singularity installed on your system to use this workflow (see above).

In order to use the snakemake workflow, please type first the following commands:

```
     cd $WORKING_DIR/10x_190712_m_moFluMemB
     ln -s 04_Workflow/snakefile_AllAnalysis.yml snakefile_AllAnalysis.yml

     cd $WORKING_DIR/10x_191105_m_moFluMemB
     ln -s 04_Workflow/snakefile_AllAnalysis.yml snakefile_AllAnalysis.yml

     cd $WORKING_DIR/custom_201216_m_moFluMemB
     ln -s 04_Workflow/snakefile_AllAnalysis.yml snakefile_AllAnalysis.yml
```

To run the analysis for the three datasets, run the following commands:


* For dataset **10x_190712_m_moFluMemB**:
```
     cd $WORKING_DIR/10x_190712_m_moFluMemB
     snakemake -r --snakefile snakefile_AllAnalysis.yml --use-singularity
```

* For dataset **10x_191105_m_moFluMemB**:
```
     cd $WORKING_DIR/10x_191105_m_moFluMemB
     snakemake -r --snakefile snakefile_AllAnalysis.yml --use-singularity
```

* For dataset **custom_201216_m_moFluMemB**:
```
     cd $WORKING_DIR/custom_201216_m_moFluMemB
     snakemake -r --snakefile snakefile_AllAnalysis.yml --use-singularity
```

### (Optional) Run the analysis individually using Docker

If you have loaded the docker images (see above), you can use Rstudio in Docker to run the analysis individually.

To start a docker container, use the following command:

```
docker run -d -p 8787:8787 -v /$WORKING_DIR:/$WORKING_DIR -e PASSWORD=<PASSWORD> -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g)  <IMAGE_NAME>
```

where:

* <PASSWORD> is a simple string you will use as password to login into Rstudio
* <IMAGE_NAME> is the Docker image name to run

One started, you can open an internet browser and use the URL https://127.0.0.1:8787.

At the login prompt, enter the name of the user session you are connected with and the password you type in place of <PASSWORD>. You are now in a Rstudio environment and the container is able to connect to the **WORKING_DIR**
of your system. Inside you will find the project files. To tun the analysis, look at the scripts "launch_report_compilation.R" that will indicates you what to do.






