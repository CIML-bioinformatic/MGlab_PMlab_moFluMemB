###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used byt the current
# analysis
#

LITERAL_TITLE = "Compare cells in different tissues"
ANALYSIS_STEP_NAME = "09_CompareTissues"
PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

# Path to the result of QC and demultiplexing analysis
PATH_SEURAT_FILTERED_OBJECT = file.path( PATH_EXPERIMENT_OUTPUT, "04_QC_Demux")

# Cluster of cells in the previous analysis identified as affected by dissociation process stress
DISSOCIATION_AFFECTED_CELLS_CLUSTER = list( LG = 2, LN = NA, SP = NA)

# The list of all markers genes per tissue and per cluster in tissue
PATH_ALL_MARKER_GENES_TABLE_FILE = list(
  LG = file.path( PATH_EXPERIMENT_OUTPUT, "/07_GlobalHeterogeneityNoBCRByTissue/allMarkers_LG.tab"),
  LN = file.path( PATH_EXPERIMENT_OUTPUT, "/07_GlobalHeterogeneityNoBCRByTissue/allMarkers_LN.tab"),
  SP = file.path( PATH_EXPERIMENT_OUTPUT, "/07_GlobalHeterogeneityNoBCRByTissue/allMarkers_SP.tab")
)

# The list of top markers genes per tissue and per cluster in tissue
PATH_TOP_MARKER_GENES_TABLE_FILE = list(
  LG = file.path( PATH_EXPERIMENT_OUTPUT, "/07_GlobalHeterogeneityNoBCRByTissue/topMarkersDT_LG.tab"),
  LN = file.path( PATH_EXPERIMENT_OUTPUT, "/07_GlobalHeterogeneityNoBCRByTissue/topMarkersDT_LN.tab"),
  SP = file.path( PATH_EXPERIMENT_OUTPUT, "/07_GlobalHeterogeneityNoBCRByTissue/topMarkersDT_SP.tab")
)

# The list of UMAP embedding per tissue
PATH_UMAP_EMBEDDING_FILE = list(
  LG = file.path( PATH_EXPERIMENT_OUTPUT, "/07_GlobalHeterogeneityNoBCRByTissue/UMAP_Embedding_LG.tab"),
  LN = file.path( PATH_EXPERIMENT_OUTPUT, "/07_GlobalHeterogeneityNoBCRByTissue/UMAP_Embedding_LN.tab"),
  SP = file.path( PATH_EXPERIMENT_OUTPUT, "/07_GlobalHeterogeneityNoBCRByTissue/UMAP_Embedding_SP.tab")
)

# The list of cluster mapping per tissue
PATH_CLUSTER_MAPPING_FILE = list(
  LG = file.path( PATH_EXPERIMENT_OUTPUT, "/07_GlobalHeterogeneityNoBCRByTissue/clusters_LG.tab"),
  LN = file.path( PATH_EXPERIMENT_OUTPUT, "/07_GlobalHeterogeneityNoBCRByTissue/clusters_LN.tab"),
  SP = file.path( PATH_EXPERIMENT_OUTPUT, "/07_GlobalHeterogeneityNoBCRByTissue/clusters_SP.tab")
)

# The list of variable genes expression file per tissue
PATH_VARIABLE_GENES_EXPRESSION_FILE = list(
  LG = file.path( PATH_EXPERIMENT_OUTPUT, "/07_GlobalHeterogeneityNoBCRByTissue/normalizedExpressions_LG.tsv"),
  LN = file.path( PATH_EXPERIMENT_OUTPUT, "/07_GlobalHeterogeneityNoBCRByTissue/normalizedExpressions_LN.tsv"),
  SP = file.path( PATH_EXPERIMENT_OUTPUT, "/07_GlobalHeterogeneityNoBCRByTissue/normalizedExpressions_SP.tsv")
)

# List of marker genes used for the study
# MARKER_GENES_LIST = c( "Fcer2a", "Sell", "Serpinb1a", "Dusp2",
# "Cnn3", "Lmo2", "Cd22", "Plac8", "Itga4", "Cxcr3", "Cst3", "Slpi", "Fscn1", "Ms4a6c", "Lck", "Fcmr", "Dtx1",
# "Nfatc1", "Klf2", "Marcksl1", "Rgs1", "Ahr", "Lmnb1", "S100a6", "S1pr1", "Zbtb32", "Cd55", "Ccr6", "Lta",
# "Ffar2", "Vim", "Tuba1a", "Apoe", "Cd24a", "Ccrl2", "S1pr3", "Il9r", "Tbx21", "Adgre1", "Itgax"
# )
MARKER_GENES_LIST = c( "Ccr6", "Cnn3", "Lmo2",
"Fcer2a", "Serpinb1a", "Cxcr3", "Fscn1", "Plac8", "Cst3", "Itga4", "Ahr", "Rgs1", "Lmnb1", "Lta", "Ffar2",
"S100a6", "Cd55", "Klf2", "Vim", "Tuba1a", "Marcksl1", "Apoe", "Cd24a", "Fcmr", "Ms4a6c", "Lck", "Nfatc1", "Zbtb32",
"Cr2", "Il9r", "Dtx1", "S1pr1", "S1pr3", "Sell", "Cd22", "Dusp2", "Slpi", "Adgre1", "Tbx21", "Itgax"
)

# Parameters for identification of marker annotations for clusters (see Seurat::FindAllMarkers())
FINDMARKERS_METHOD    = "wilcox"  # Method used to identify markers
FINDMARKERS_ONLYPOS   = TRUE;     # Only consider overexpressed annotations for markers ? (if FALSE downregulated genes can also be markers)
FINDMARKERS_MINPCT    = 0.1;      # Only test genes that are detected in a minimum fraction of cells in either of the two populations. Speed up the function by not testing genes that are very infrequently expressed. Default is '0.1'.
FINDMARKERS_LOGFC_THR = 0.25;     # Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is '0.25'. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
FINDMARKERS_PVAL_THR  = 0.001;    # PValue threshold for identification of significative markers
FINDMARKERS_SHOWTOP   = 10;       # Number of marker genes to show in report and tables (NULL for all)