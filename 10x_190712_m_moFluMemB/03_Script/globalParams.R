###############################################################################
# This file defines PROJECT parameters as global variables that will be loaded
# before analysis starts. It should define common parameters shared by several
# samples. 
#
# Parameters for individual samples (e.g. SAMPLE_NAME, PATH_DATA) must be set
# in the file 'sampleParams.R' (stored in the sample folder). Any value defined 
# there will override values defined in this file.
#




#### General

PROJECT_NAME = "MGlab Mouse Influenza Memory B cells - data 2019/07/12";
PROJECT_SHORT_NAME = "moFluMemB"

PATH_PROJECT = "/mnt/NAS7/BIP/Temp/moFluMemB"

EXPERIMENT_NAME = "10x_190712_m_moFluMemB"

#### Input / Output

# Output folder name in data folder (for R session object, lists of cells/genes) 
PATH_EXPERIMENT_OUTPUT = file.path( PATH_PROJECT, EXPERIMENT_NAME, "05_Output")

#### Debug

.SHOWFLEXBORDERS = FALSE;
.VERBOSE = FALSE;



