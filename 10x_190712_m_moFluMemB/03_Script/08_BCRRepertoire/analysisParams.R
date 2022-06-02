###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used byt the current
# analysis
#

LITERAL_TITLE = "Analysis of BCR repertoire"
ANALYSIS_STEP_NAME = "08_BCRRepertoire"
PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

# Path to the output of migmap to get the details on sequences of BCR
PATH_BCR_MIGMAP_LIGHTCHAIN_ANALYSIS_FILE = file.path( PATH_EXPERIMENT_OUTPUT, "03_BCRAnalysis/4_migmap/IGKLV_filtered_for_all.csv")
PATH_BCR_MIGMAP_HEAVYCHAIN_ANALYSIS_FILE = file.path( PATH_EXPERIMENT_OUTPUT, "03_BCRAnalysis/4_migmap/IGHV_filtered_for_all.csv")

PATH_BCR_BLAST_LIGHTCHAIN_ANALYSIS_FILE = file.path( PATH_EXPERIMENT_OUTPUT, "03_BCRAnalysis/4_blast/10935372BCR_S6_blastn_csr_LightCstRegion.out")
PATH_BCR_BLAST_HEAVYCHAIN_ANALYSIS_FILE = file.path( PATH_EXPERIMENT_OUTPUT, "03_BCRAnalysis/4_blast/10935372BCR_S6_blastn_csr_HeavyCstRegion.out")

# Cluster markers genes per tissues
PATH_CLUSTER_MARKER_GENES_PER_TISSUE = list(  LG = file.path( PATH_EXPERIMENT_OUTPUT, "07_GlobalHeterogeneityNoBCRByTissue/allMarkers_LG.tab"),
                                              LN = file.path( PATH_EXPERIMENT_OUTPUT, "07_GlobalHeterogeneityNoBCRByTissue/allMarkers_LN.tab"),
                                              SP = file.path( PATH_EXPERIMENT_OUTPUT, "07_GlobalHeterogeneityNoBCRByTissue/allMarkers_SP.tab")
)

PATH_CLUSTER_PAIRWISE_MARKER_GENES_PER_TISSUE = list(  LG = file.path( PATH_EXPERIMENT_OUTPUT, "07_GlobalHeterogeneityNoBCRByTissue/pairwiseMarkers_LG.tab"),
                                              LN = file.path( PATH_EXPERIMENT_OUTPUT, "07_GlobalHeterogeneityNoBCRByTissue/pairwiseMarkers_LN.tab"),
                                              SP = file.path( PATH_EXPERIMENT_OUTPUT, "07_GlobalHeterogeneityNoBCRByTissue/pairwiseMarkers_SP.tab")
)

# The length of the FR4 kept for the similiaty sequence analysis
FR4_AA_LENGTH_CUT = 36

# Path to the Seurat object to use
PATH_SEURAT_FILTERED_OBJECT = file.path( PATH_EXPERIMENT_OUTPUT, "06_GlobalHeterogeneityNoBCR")

PATH_SEURAT_FILTERED_OBJECT_BY_TISSUE = list( 
                                              LG = file.path( PATH_EXPERIMENT_OUTPUT, "07_GlobalHeterogeneityNoBCRByTissue", "sc10x.rna.seurat_LG.RDS"),
                                              LN = file.path( PATH_EXPERIMENT_OUTPUT, "07_GlobalHeterogeneityNoBCRByTissue", "sc10x.rna.seurat_LN.RDS"),
                                              SP = file.path( PATH_EXPERIMENT_OUTPUT, "07_GlobalHeterogeneityNoBCRByTissue", "sc10x.rna.seurat_SP.RDS")
                                          )
# Boolean indicating if the ClustalOmega must be done
# If FALSE, the analysis will try to get file from previous ClustalOmega run
# If TRUE, the ClustalOmega analysis will be launched and the possible previous resutls will be removed
FORCE_CLUSTALO = FALSE

# The minium cell number defining a "large cluster" when all tissues together
GLOBAL_LARGE_CLUSTER_THRESHOLD = 20

# The minium cell number defining a "large cluster" when considering tissues separatly
LARGE_CLUSTER_THRESHOLD = 10
GLOBAL_LARGE_CLUSTER_THRESHOLD = 20


# Parameters for identification of marker annotations for clusters (see Seurat::FindAllMarkers())
FINDMARKERS_LOGFC_THR = 0.25;     # Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is '0.25'. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
FINDMARKERS_PVAL_THR  = 0.001;    # PValue threshold for identification of significative markers
FINDMARKERS_SHOWTOP   = 10;       # Number of marker genes to show in report and tables (NULL for all)
