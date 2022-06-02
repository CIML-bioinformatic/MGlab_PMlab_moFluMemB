###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used byt the current
# analysis
#

ANALYSIS_STEP_NAME = "04_StringentQC"
LITERAL_TITLE = "Quality Control"

PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME);

# Path to the file containing the cell gene expression matrix
PATH_RNA_COUNT_TABLE = file.path( PATH_EXPERIMENT_OUTPUT, "/03_BCRAnalysis/A_Preprocessing/all_endogenes_UMI.csv")

# Path to cells metadata issued from pre-processing
PATH_CELL_PREPROCESSING_METADATA = file.path( PATH_EXPERIMENT_OUTPUT, "03_BCRAnalysis/D_MergingData/nofiltered_metadata_allFeatures.csv")

#### Filtering / Normalization

LOAD_MIN_CELLS     = 0;    # Retain cells with at least this many features (annotations)
LOAD_MIN_FEATURES  = 0;  # Retain annotations appearing in at least this many cells

# Cells with number of UMIs outside the range will be excluded
FILTER_UMI_MIN     = NULL;
FILTER_UMI_MAX     = NULL;

# Cells with number of genes outside the range will be excluded
FILTER_FEATURE_MIN = 500;
FILTER_FEATURE_MAX = NULL;

# Cells with percentage of mitochondrial gene above this value will be excluded
FILTER_MITOPCT_MAX = 3;

# Cells with percentage of ribosomal gene below this value will be excluded
FILTER_RIBOPCT_MIN = 5;

# Cells with Percentage of ERCC gene above this value will be xcluded
FILTER_ERCC_MAX = 10;

# Cells with a Pearson correlation value below this value will be excluded
FILTER_ACCURACYPEARSONERCC_MIN = 0.7

# Normalization parameters (see Seurat::NormalizeData())
DATA_NORM_METHOD      = "LogNormalize";
DATA_NORM_SCALEFACTOR = 10000;

# Scaling parameters (see Seurat::ScaleData())
DATA_CENTER       = TRUE;
DATA_SCALE        = FALSE;
DATA_VARS_REGRESS = NULL;  # c("nCount_RNA") for UMIs (NULL to ignore)

