###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used byt the current
# analysis
#

LITERAL_TITLE = "Analysis of BCR repertoire"
ANALYSIS_STEP_NAME = "08_BCRRepertoire"
PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

# ##########################################
# VARIABLES ON DATA SOURCES
# ##########################################

# path to the seurat object to use for analyse
PATH_SEURAT_FILTERED_OBJECT = file.path( PATH_EXPERIMENT_OUTPUT, "06_GlobalHeterogeneityNoBCR")

# Path to the output of migmap to get the details on sequences of BCR
PATH_BCR_MIGMAP_LIGHTCHAIN_ANALYSIS_FILE = file.path( PATH_EXPERIMENT_OUTPUT, "03_BCRAnalysis/C_BCR_retrieval/IGKLV_filtered_for_all.csv")
PATH_BCR_MIGMAP_HEAVYCHAIN_ANALYSIS_FILE = file.path( PATH_EXPERIMENT_OUTPUT, "03_BCRAnalysis/C_BCR_retrieval/IGHV_filtered_for_all.csv")

# Path to cells metadata issued from pre-processing
PATH_CELL_PREPROCESSING_METADATA = file.path( PATH_EXPERIMENT_OUTPUT, "03_BCRAnalysis/D_MergingData/nofiltered_metadata_allFeatures.csv")

# ##########################################
# VARIABLES FOR HETEROGENEITY ANALYSIS
# ##########################################

# Variable value for QC filter
MINIMUM_FEATURE_RNA = 200
MAXIMUM_PERCENT_MITO = 0.1
MAXIMUM_PERCENT_ERCC = 0.1
MINIMUM_ACCURACYPEARSONERCC = 0.7

# The length of the FR4 kept for the similiaty sequence analysis
FR4_AA_LENGTH_CUT = 36

# Path to the Seurat object to use
PATH_SEURAT_FILTERED_OBJECT = file.path( PATH_EXPERIMENT_OUTPUT, "06_GlobalHeterogeneityNoBCR")

PATH_SEURAT_FILTERED_OBJECT_BY_TISSUE = list( 
                                              LG = file.path( PATH_EXPERIMENT_OUTPUT, "06_GlobalHeterogeneityNoBCR", "custom.rna.seurat.RDS")
                                            )

# Cluster markers genes per tissues
PATH_CLUSTER_MARKER_GENES_PER_TISSUE = list(  
                                              LG = file.path( PATH_EXPERIMENT_OUTPUT, "06_GlobalHeterogeneityNoBCR", "allMarkers.tab")
                                            )


# Boolean indicating if the ClustalOmega must be done
# If FALSE, the analysis will try to get file from previous ClustalOmega run
# If TRUE, the ClustalOmega analysis will be launched and the possible previous resutls will be removed
FORCE_CLUSTALO = FALSE

# The minium cell number defining a "large cluster" when all tissues together
GLOBAL_LARGE_CLUSTER_THRESHOLD = 5

# The minium cell number defining a "large cluster" when considering tissues separatly
LARGE_CLUSTER_THRESHOLD = 5

# Parameters for identification of marker annotations for clusters (see Seurat::FindAllMarkers())
FINDMARKERS_LOGFC_THR = 0.25;     # Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is '0.25'. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
FINDMARKERS_PVAL_THR  = 0.001;    # PValue threshold for identification of significative markers
FINDMARKERS_SHOWTOP   = 10;       # Number of marker genes to show in report and tables (NULL for all)

# ##########################################
# VARIABLES FOR PHENOTYPE ANALYSIS
# ##########################################

# -- Name of the index sort data for each phenotype
PHENOTYPE_CXCR3_INDEXSORT_SOURCE = "asinh_20_CXCR3.BV421"
PHENOTYPE_CCR6_INDEXSORT_SOURCE = "asinh_10_CCR6.PEDz594"
PHENOTYPE_HA_INDEXSORT_SOURCE = "asinh_100_HA.PE"
PHENOTYPE_NP_INDEXSORT_SOURCE = "asinh_100_NP.APC"

# -- Gate threshold for each phenotype
PHENOTYPE_ASINH_CXCR3_THRESHOLD = 0.5   #Limit to declare a cell CXCR3 positive  (0.1 for scale 100, 0.2 for scale 50, 0.5 for scale 20)
PHENOTYPE_ASINH_CCR6_THRESHOLD = 1.5    #Limit to delcare a cell CCR6 positive
PHENOTYPE_ASINH_HA_THRESHOLD=0.8        #Limit to delare a cell HA positive
PHENOTYPE_ASINH_NP_THRESHOLD=0.5        #Limit to declare a cell NP positive

# -- Names of cell phenotype for each gate
PHENOTYPE_CXCR3_NAMES = c( Pos =  "CXCR3p", Neg =  "CXCR3n")
PHENOTYPE_CCR6_NAMES = c( Pos =  "CCR6p", Neg =  "CCR6n")
PHENOTYPE_CXCR3xCCR6_NAMES = c( PosPos =  "CXCR3p.CCR6p", PosNeg =  "CXCR3p.CCR6n", NegPos =  "CXCR3n.CCR6p", NegNeg =  "CXCR3n.CCR6n")
PHENOTYPE_CXCR3xCCR6_NAMES_ORDER = c( PHENOTYPE_CXCR3xCCR6_NAMES[ "PosPos"], PHENOTYPE_CXCR3xCCR6_NAMES[ "PosNeg"], 
                                      PHENOTYPE_CXCR3xCCR6_NAMES[ "NegPos"], PHENOTYPE_CXCR3xCCR6_NAMES[ "NegNeg"])

PHENOTYPE_HA_NAMES = c( Pos =  "HAp", Neg =  "HAn")
PHENOTYPE_NP_NAMES = c( Pos =  "NPp", Neg =  "NPn")
PHENOTYPE_HAxNP_NAMES = c( PosPos =  "HAp.NPp", PosNeg =  "HAp.NPn", NegPos =  "HAn.NPp", NegNeg =  "HAn.NPn")
PHENOTYPE_HAxNP_NAMES_ORDER = c( PHENOTYPE_HAxNP_NAMES[ "PosPos"], PHENOTYPE_HAxNP_NAMES[ "PosNeg"], 
                                      PHENOTYPE_HAxNP_NAMES[ "NegPos"], PHENOTYPE_HAxNP_NAMES[ "NegNeg"])

