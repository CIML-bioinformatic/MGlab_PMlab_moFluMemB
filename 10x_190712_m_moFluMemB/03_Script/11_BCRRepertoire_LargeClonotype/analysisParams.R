###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used byt the current
# analysis
#

LITERAL_TITLE = "Analysis of BCR repertoire on large clonotypes"
ANALYSIS_STEP_NAME = "11_BCRRepertoire_LargeClonotype"
PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

# Cluster markers genes per tissues
PATH_LARGE_CLONOTYPE_CELL_FILE = file.path( PATH_EXPERIMENT_OUTPUT, "08_BCRRepertoire", "")
